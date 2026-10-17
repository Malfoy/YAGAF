//! Project exact seed matches to gap boundaries through known flank sequence.
use crate::{index::Sequence, is_base};

pub(crate) struct ReadView<'a> {
    pub bytes: &'a [u8],
    pub reverse: bool,
}
impl ReadView<'_> {
    pub fn len(&self) -> usize {
        self.bytes.len()
    }
    pub fn base(&self, position: usize) -> u8 {
        if !self.reverse {
            return self.bytes[position].to_ascii_uppercase();
        }
        match self.bytes[self.len() - 1 - position].to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        }
    }
    pub fn sequence(&self, start: usize, end: usize) -> Option<Sequence> {
        let raw = if self.reverse {
            &self.bytes[self.len() - end..self.len() - start]
        } else {
            &self.bytes[start..end]
        };
        if !raw.iter().all(|&b| is_base(b)) {
            return None;
        }
        Some((start..end).map(|i| self.base(i)).collect())
    }
}
#[derive(Default)]
pub(crate) struct AlignmentScratch {
    previous: Vec<u32>,
    current: Vec<u32>,
}

/// Precomputed Myers bit masks keep the usual <=100-base flank alignments in
/// registers. Longer explicitly requested windows use the banded DP fallback.
#[derive(Debug, PartialEq, Eq)]
struct Extension {
    consumed: usize,
    edits: usize,
}

#[derive(Clone)]
pub(crate) struct Pattern {
    length: usize,
    masks: [u128; 4],
    long: Vec<u8>,
    edit_limit: Option<usize>,
}
fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}
impl Pattern {
    pub fn new(reference: &[u8], reverse: bool) -> Self {
        let mut p = Self {
            length: reference.len(),
            masks: [0; 4],
            long: Vec::new(),
            edit_limit: None,
        };
        for i in 0..reference.len() {
            let b =
                reference[if reverse { reference.len() - 1 - i } else { i }].to_ascii_uppercase();
            p.long.push(b);
            if reference.len() <= 128
                && let Some(base) = base_index(b)
            {
                p.masks[base] |= 1u128 << i;
            }
        }
        p
    }

    pub fn limit_edits(&mut self, limit: usize) {
        self.edit_limit = Some(limit);
    }

    /// Anchored edit alignment to a free read endpoint. Projection requires a
    /// unique best endpoint; validation only requires a score within tolerance.
    pub fn extend(
        &self,
        read: &ReadView<'_>,
        start: usize,
        forward: bool,
        unique_endpoint: bool,
        scratch: &mut AlignmentScratch,
    ) -> Option<usize> {
        self.align(read, start, forward, unique_endpoint, scratch)
            .map(|a| a.consumed)
    }

    fn align(
        &self,
        read: &ReadView<'_>,
        start: usize,
        forward: bool,
        unique_endpoint: bool,
        scratch: &mut AlignmentScratch,
    ) -> Option<Extension> {
        if self.length == 0 {
            return Some(Extension {
                consumed: 0,
                edits: 0,
            });
        }
        let available = if forward {
            read.len().checked_sub(start)?
        } else {
            start
        };
        let target = |j| read.base(if forward { start + j } else { start - 1 - j });
        let edits = self
            .edit_limit
            .map_or(self.length / 5 + 2, |limit| limit.min(self.length / 5 + 2));
        if self.length > available.saturating_add(edits) {
            return None;
        }
        let width = available.min(self.length.saturating_add(edits));
        if available >= self.length
            && self
                .long
                .iter()
                .enumerate()
                .all(|(i, &b)| is_base(b) && b == target(i))
        {
            return Some(Extension {
                consumed: self.length,
                edits: 0,
            });
        }
        if edits == 0 {
            return None;
        }
        if edits == 1 {
            let first = (0..available.min(self.length))
                .take_while(|&i| is_base(self.long[i]) && self.long[i] == target(i))
                .count();
            let same = |query: usize, read_start: usize| {
                (query..self.length).enumerate().all(|(offset, i)| {
                    is_base(self.long[i]) && self.long[i] == target(read_start + offset)
                })
            };
            let mut endpoint: Option<usize> = None;
            let mut ties = 0;
            for (valid, consumed) in [
                (
                    available >= self.length && same(first + 1, first + 1),
                    self.length,
                ),
                (
                    available > self.length && same(first, first + 1),
                    self.length + 1,
                ),
                (
                    available + 1 >= self.length && same(first + 1, first),
                    self.length - 1,
                ),
            ] {
                if valid {
                    endpoint = Some(endpoint.map_or(consumed, |old| old.min(consumed)));
                    ties += 1;
                }
            }
            if unique_endpoint && ties != 1 {
                return None;
            }
            return endpoint.map(|consumed| Extension { consumed, edits: 1 });
        }
        // With <=d edits, at least one of d+1 disjoint reference blocks
        // survives exactly, within d bases of its original offset.
        if (2..=4).contains(&edits) && self.length >= 30 {
            let parts = edits + 1;
            let possible = (0..parts).any(|part| {
                let begin = part * self.length / parts;
                let end = (part + 1) * self.length / parts;
                let block = &self.long[begin..end];
                (begin.saturating_sub(edits)
                    ..=begin
                        .saturating_add(edits)
                        .min(width.saturating_sub(block.len())))
                    .any(|pos| {
                        pos + block.len() <= width
                            && block
                                .iter()
                                .enumerate()
                                .all(|(j, &b)| is_base(b) && b == target(pos + j))
                    })
            });
            if !possible {
                return None;
            }
        }
        if self.length > 128 {
            return self.dynamic(&target, width, edits, unique_endpoint, scratch);
        }
        let high = 1u128 << (self.length - 1);
        let mut positive = !0u128;
        let mut negative = 0u128;
        let mut score = self.length;
        let mut best = score;
        let mut endpoint = 0;
        let mut ties = 1;
        for j in 0..width {
            let equal = base_index(target(j)).map_or(0, |b| self.masks[b]);
            let vertical = equal | negative;
            let horizontal = ((equal & positive).wrapping_add(positive) ^ positive) | equal;
            let mut ph = negative | !(horizontal | positive);
            let mut mh = positive & horizontal;
            score += usize::from(ph & high != 0);
            score -= usize::from(mh & high != 0);
            ph = (ph << 1) | 1;
            mh <<= 1;
            positive = mh | !(vertical | ph);
            negative = ph & vertical;
            if score < best {
                best = score;
                endpoint = j + 1;
                ties = 1;
            } else if score == best {
                ties += 1;
            }
        }
        (best <= edits && (!unique_endpoint || ties == 1)).then_some(Extension {
            consumed: endpoint,
            edits: best,
        })
    }

    fn dynamic(
        &self,
        target: &impl Fn(usize) -> u8,
        width: usize,
        edits: usize,
        unique_endpoint: bool,
        scratch: &mut AlignmentScratch,
    ) -> Option<Extension> {
        const INF: u32 = u32::MAX / 4;
        scratch.previous.resize(width + 1, INF);
        scratch.current.resize(width + 1, INF);
        scratch.previous.fill(INF);
        for j in 0..=width.min(edits) {
            scratch.previous[j] = j as u32;
        }
        for i in 1..=self.length {
            scratch.current.fill(INF);
            if i <= edits {
                scratch.current[0] = i as u32;
            }
            let first = i.saturating_sub(edits).max(1);
            let last = width.min(i.saturating_add(edits));
            for j in first..=last {
                let substitution =
                    scratch.previous[j - 1] + u32::from(self.long[i - 1] != target(j - 1));
                scratch.current[j] = substitution
                    .min(scratch.previous[j] + 1)
                    .min(scratch.current[j - 1] + 1);
            }
            std::mem::swap(&mut scratch.previous, &mut scratch.current);
        }
        let best = *scratch.previous.iter().min()?;
        if best > edits as u32 {
            return None;
        }
        let mut endpoints = scratch
            .previous
            .iter()
            .enumerate()
            .filter(|&(_, score)| *score == best)
            .map(|(j, _)| j);
        let endpoint = endpoints.next()?;
        if unique_endpoint && endpoints.next().is_some() {
            return None;
        }
        Some(Extension {
            consumed: endpoint,
            edits: best as usize,
        })
    }
}

pub(crate) fn project(
    pattern: &Pattern,
    side: usize,
    seed_position: usize,
    k: usize,
    read: &ReadView<'_>,
    scratch: &mut AlignmentScratch,
) -> Option<usize> {
    if side == 0 {
        let start = seed_position.checked_add(k)?;
        start.checked_add(pattern.extend(read, start, true, true, scratch)?)
    } else {
        seed_position.checked_sub(pattern.extend(read, seed_position, false, true, scratch)?)
    }
}

pub(crate) fn context_score(
    context: &[u8],
    pattern: &Pattern,
    side: usize,
    boundary: usize,
    k: usize,
    read: &ReadView<'_>,
    scratch: &mut AlignmentScratch,
) -> Option<usize> {
    let available = if side == 0 {
        boundary
    } else {
        read.len().saturating_sub(boundary)
    };
    let length = context.len().min(available);
    if length < k {
        return None;
    }
    if length == context.len() {
        return pattern
            .align(read, boundary, side == 1, false, scratch)
            .map(|a| a.edits);
    }
    let reference = if side == 0 {
        &context[context.len() - length..]
    } else {
        &context[..length]
    };
    let mut cropped = Pattern::new(reference, side == 0);
    cropped.edit_limit = pattern.edit_limit;
    cropped
        .align(read, boundary, side == 1, false, scratch)
        .map(|a| a.edits + context.len() - length)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tied_projection_endpoints_and_uncovered_flanks_are_rejected() {
        let view = ReadView {
            bytes: b"ATC",
            reverse: false,
        };
        let mut scratch = AlignmentScratch::default();
        let pattern = Pattern::new(b"AC", false);
        assert!(pattern.extend(&view, 0, true, true, &mut scratch).is_none());
        assert!(
            pattern
                .extend(&view, 0, true, false, &mut scratch)
                .is_some()
        );
        assert!(
            context_score(
                b"ACGTACGT",
                &Pattern::new(b"ACGTACGT", true),
                0,
                3,
                5,
                &view,
                &mut scratch
            )
            .is_none()
        );
    }
    #[test]
    fn bit_parallel_alignment_matches_dynamic_programming() {
        let mut state = 719423437812u64;
        let mut dna = |len: usize| {
            (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGTN"[state as usize % 5]
                })
                .collect::<Vec<_>>()
        };
        for length in [1, 2, 15, 31, 63, 64, 65, 100, 127, 128, 129, 200] {
            for _ in 0..80 {
                let reference = dna(length);
                let mut read = reference.clone();
                for pos in (0..read.len()).step_by(11) {
                    read[pos] = b"ACTG"[pos % 4];
                }
                read.splice(length / 2..length / 2, b"AT".iter().copied());
                read.extend(dna(30));
                let mut pattern = Pattern::new(&reference, false);
                // Unknown reference symbols cannot match in either implementation.
                pattern.long = reference
                    .iter()
                    .map(|&b| if b == b'N' { b'X' } else { b })
                    .collect();
                let view = ReadView {
                    bytes: &read,
                    reverse: false,
                };
                for limit in [None, Some(0), Some(1), Some(2), Some(4)] {
                    pattern.edit_limit = limit;
                    let edits = limit.map_or(length / 5 + 2, |n| n.min(length / 5 + 2));
                    let width = read.len().min(length + edits);
                    for unique in [false, true] {
                        let actual =
                            pattern.align(&view, 0, true, unique, &mut AlignmentScratch::default());
                        let expected = pattern.dynamic(
                            &|j| read[j],
                            width,
                            edits,
                            unique,
                            &mut AlignmentScratch::default(),
                        );
                        assert_eq!(
                            actual, expected,
                            "length={length} unique={unique} limit={limit:?}"
                        );
                    }
                }
            }
        }
    }
}
