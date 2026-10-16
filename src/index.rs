// Adapted from Malfoy/K2Rmini src/bin/k2rminimulti.rs (GPL-3.0).
// Pinned upstream revision and adaptation details are in THIRD_PARTY.md.
use crate::{Gap, Scaffold, is_base, revcomp};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use simd_minimizers::minimizers;
use simd_minimizers::packed_seq::{PackedSeqVec, Seq, SeqVec};
use simd_minimizers::seq_hash::{KmerHasher, NtHasher};

struct Anchor {
    gap: usize,
    side: usize,
    reverse: bool,
    bases: Vec<u8>,
    minimizer_offset: usize,
}

pub struct Index {
    k: usize,
    m: usize,
    anchors: Vec<Anchor>,
    minis: FxHashMap<u32, Vec<usize>>,
    kmers: FxHashMap<u32, Vec<usize>>,
    use_offsets: bool,
}

#[derive(Default)]
pub struct Buffers {
    packed: PackedSeqVec,
    mini_positions: Vec<u32>,
    hashes: Vec<u32>,
}

#[derive(Clone, Copy, Default)]
struct Hit {
    count: u8,
    position: usize,
}

#[derive(Default)]
pub struct Scanner {
    buffers: Buffers,
    hits: FxHashMap<(usize, bool), [Hit; 2]>,
}

#[derive(Debug)]
pub enum Span {
    Sequence(Vec<u8>),
    AmbiguousPair,
    Rejected(&'static str),
}

impl Index {
    pub fn new(gaps: &[Gap], k: usize, m: usize) -> Self {
        let mut index = Self {
            k,
            m,
            anchors: Vec::new(),
            minis: FxHashMap::default(),
            kmers: FxHashMap::default(),
            use_offsets: false,
        };
        let builder = minimizers(m, k - m + 1);
        let hasher = NtHasher::<false>::new(k);
        let mut jobs = Vec::new();
        for (gap, item) in gaps
            .iter()
            .enumerate()
            .filter(|(_, g)| g.rejection.is_none())
        {
            if let Some(anchors) = &item.anchors {
                for (side, bases) in anchors.iter().enumerate() {
                    for reverse in [false, true] {
                        jobs.push((gap, side, reverse, bases));
                    }
                }
            }
        }
        let prepared: Vec<_> = jobs
            .into_par_iter()
            .with_min_len(64)
            .map_init(Buffers::default, |buffers, (gap, side, reverse, bases)| {
                let bases = if reverse {
                    revcomp(bases)
                } else {
                    bases.clone()
                };
                buffers.packed.clear();
                buffers.packed.push_ascii(&bases);
                buffers.mini_positions.clear();
                builder.run(buffers.packed.as_slice(), &mut buffers.mini_positions);
                let minimizer_offset = buffers.mini_positions[0] as usize;
                let minimizer_word = buffers
                    .packed
                    .slice(minimizer_offset..minimizer_offset + m)
                    .as_u64() as u32;
                buffers.hashes.clear();
                hasher
                    .hash_kmers_simd(buffers.packed.as_slice(), 1)
                    .collect_into(&mut buffers.hashes);
                (
                    Anchor {
                        gap,
                        side,
                        reverse,
                        bases,
                        minimizer_offset,
                    },
                    minimizer_word,
                    buffers.hashes[0],
                )
            })
            .collect();
        index.anchors.reserve(prepared.len());
        index.minis.reserve(prepared.len());
        index.kmers.reserve(prepared.len());
        for (anchor, minimizer, hash) in prepared {
            let id = index.anchors.len();
            index.minis.entry(minimizer).or_default().push(id);
            index.kmers.entry(hash).or_default().push(id);
            index.anchors.push(anchor);
        }
        index.use_offsets = m >= 9 && index.minis.values().all(|postings| postings.len() <= 64);
        index
    }

    /// Preserve anchor IDs and remove rejected postings without recomputing
    /// minimizers or k-mer hashes for the surviving anchors.
    pub fn retain_eligible(&mut self, gaps: &[Gap]) {
        for postings in self.minis.values_mut().chain(self.kmers.values_mut()) {
            postings.retain(|&id| gaps[self.anchors[id].gap].rejection.is_none());
        }
        self.minis.retain(|_, postings| !postings.is_empty());
        self.kmers.retain(|_, postings| !postings.is_empty());
        self.use_offsets = self.m >= 9 && self.minis.values().all(|postings| postings.len() <= 64);
    }

    // Each ACTG run is independent. Never join across an ambiguous base: that
    // would introduce nonexistent k-mers and corrupt read coordinates.
    fn visit(&self, sequence: &[u8], buffers: &mut Buffers, mut found: impl FnMut(usize, usize)) {
        if self.minis.is_empty() {
            return;
        }
        let builder = minimizers(self.m, self.k - self.m + 1);
        let hasher = NtHasher::<false>::new(self.k);
        let mut offset = 0;
        for run in sequence.split(|&base| !is_base(base)) {
            if run.len() >= self.k && !self.kmers.is_empty() {
                buffers.packed.clear();
                buffers.packed.push_ascii(run);
                buffers.mini_positions.clear();
                builder.run(buffers.packed.as_slice(), &mut buffers.mini_positions);
                if self.use_offsets {
                    // An exact anchor has the same minimizer at the same offset.
                    // Its postings therefore identify candidate starts directly.
                    // Every candidate still requires an exact base comparison.
                    for &pos in &buffers.mini_positions {
                        let key = buffers
                            .packed
                            .slice(pos as usize..pos as usize + self.m)
                            .as_u64() as u32;
                        if let Some(ids) = self.minis.get(&key) {
                            for &id in ids {
                                let anchor = &self.anchors[id];
                                if let Some(start) =
                                    (pos as usize).checked_sub(anchor.minimizer_offset)
                                    && start + self.k <= run.len()
                                    && run[start..start + self.k]
                                        .eq_ignore_ascii_case(&anchor.bases)
                                {
                                    found(id, offset + start);
                                }
                            }
                        }
                    }
                    offset += run.len() + 1;
                    continue;
                }
                let candidate = buffers.mini_positions.iter().any(|&pos| {
                    let key = buffers
                        .packed
                        .slice(pos as usize..pos as usize + self.m)
                        .as_u64() as u32;
                    self.minis.contains_key(&key)
                });
                if candidate {
                    buffers.hashes.clear();
                    hasher
                        .hash_kmers_simd(buffers.packed.as_slice(), 1)
                        .collect_into(&mut buffers.hashes);
                    for (pos, hash) in buffers.hashes.iter().enumerate() {
                        if let Some(ids) = self.kmers.get(hash) {
                            for &id in ids {
                                if run[pos..pos + self.k]
                                    .eq_ignore_ascii_case(&self.anchors[id].bases)
                                {
                                    found(id, offset + pos);
                                }
                            }
                        }
                    }
                }
            }
            offset += run.len() + 1;
        }
    }

    /// Strictly require each flank to occur once in the assembly, considering
    /// both strands. This costs an in-memory assembly scan, never a read rescan.
    pub fn reject_repeated_anchors(&self, scaffolds: &[Scaffold], gaps: &mut [Gap]) {
        const CHUNK: usize = 262_144;
        let jobs: Vec<_> = scaffolds
            .iter()
            .flat_map(|scaffold| {
                (0..scaffold.sequence.len()).step_by(CHUNK).map(|start| {
                    let end = (start + CHUNK).min(scaffold.sequence.len());
                    let overlap_end = (end + self.k - 1).min(scaffold.sequence.len());
                    (&scaffold.sequence[start..overlap_end], end - start)
                })
            })
            .collect();
        let (_, counts) = jobs
            .into_par_iter()
            .fold(
                || (Buffers::default(), vec![0u8; self.anchors.len()]),
                |(mut buffers, mut counts), (sequence, owned)| {
                    self.visit(sequence, &mut buffers, |id, pos| {
                        if pos < owned {
                            counts[id] = counts[id].saturating_add(1).min(2);
                        }
                    });
                    (buffers, counts)
                },
            )
            .reduce(
                || (Buffers::default(), vec![0u8; self.anchors.len()]),
                |(buffers, mut a), (_, b)| {
                    for (left, right) in a.iter_mut().zip(b) {
                        *left = (*left + right).min(2);
                    }
                    (buffers, a)
                },
            );
        for ids in (0..self.anchors.len()).step_by(2) {
            let anchor = &self.anchors[ids];
            let palindrome = anchor.bases == self.anchors[ids + 1].bases;
            let occurrences = counts[ids] + if palindrome { 0 } else { counts[ids + 1] };
            if occurrences != 1 {
                gaps[anchor.gap].rejection = Some("nonunique_scaffold_anchor");
            }
        }
    }

    pub fn extract(
        &self,
        sequence: &[u8],
        scanner: &mut Scanner,
        max_span: usize,
    ) -> Vec<(usize, Span)> {
        scanner.hits.clear();
        let hits = &mut scanner.hits;
        self.visit(sequence, &mut scanner.buffers, |id, position| {
            let anchor = &self.anchors[id];
            let hit = &mut hits.entry((anchor.gap, anchor.reverse)).or_default()[anchor.side];
            hit.position = position;
            hit.count = hit.count.saturating_add(1).min(2);
        });
        let mut result = Vec::with_capacity(scanner.hits.len());
        for (&(gap, reverse), pair) in &scanner.hits {
            if pair[0].count == 0 || pair[1].count == 0 {
                continue;
            }
            let (left, right) = if reverse {
                (pair[1], pair[0])
            } else {
                (pair[0], pair[1])
            };
            if left.count > 1 || right.count > 1 {
                result.push((gap, Span::AmbiguousPair));
                continue;
            }
            let start = left.position + self.k;
            let end = right.position;
            if start > end {
                continue;
            }
            let span = if end - start > max_span {
                Span::Rejected("span_limit")
            } else if start == end {
                Span::Rejected("empty_read_span")
            } else if !sequence[start..end].iter().all(|&b| is_base(b)) {
                Span::Rejected("ambiguous_read_span")
            } else {
                Span::Sequence(if reverse {
                    revcomp(&sequence[start..end])
                } else {
                    sequence[start..end].to_ascii_uppercase()
                })
            };
            result.push((gap, span));
        }
        result.sort_unstable_by_key(|(gap, _)| *gap);
        result.dedup_by(|later, first| {
            if later.0 == first.0 {
                first.1 = Span::AmbiguousPair;
                true
            } else {
                false
            }
        });
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizer_postings_match_exhaustive_search() {
        let mut seed = 0x123456789abcdefu64;
        let mut dna = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    b"ACGT"[(seed >> 32) as usize & 3]
                })
                .collect()
        };
        for k in [3, 5, 15, 17, 31, 33, 63] {
            let gaps: Vec<_> = (0..20)
                .map(|_| Gap {
                    scaffold: 0,
                    start: k,
                    end: k + 1,
                    anchors: Some([dna(k), dna(k)]),
                    rejection: None,
                })
                .collect();
            for m in [1, k.min(9), k.min(15), k.min(29)] {
                let mut index = Index::new(&gaps, k, m);
                let mut read = dna(2048);
                for (id, anchor) in index.anchors.iter().enumerate() {
                    let start = id * (k + 3) % (read.len() - k);
                    read[start..start + k].copy_from_slice(&anchor.bases);
                }
                for p in (127..read.len()).step_by(251) {
                    read[p] = b'N';
                }
                for p in (29..read.len()).step_by(97) {
                    read[p] = read[p].to_ascii_lowercase();
                }
                let expected = |index: &Index| {
                    let mut all = Vec::new();
                    for (id, anchor) in index.anchors.iter().enumerate() {
                        for (pos, word) in read.windows(k).enumerate() {
                            if word.eq_ignore_ascii_case(&anchor.bases) {
                                all.push((id, pos));
                            }
                        }
                    }
                    all.sort_unstable();
                    all
                };
                for offsets in [false, m >= 9] {
                    index.use_offsets = offsets;
                    let mut actual = Vec::new();
                    index.visit(&read, &mut Buffers::default(), |id, pos| {
                        actual.push((id, pos))
                    });
                    actual.sort_unstable();
                    assert_eq!(actual, expected(&index), "k={k} m={m} offsets={offsets}");
                }
            }
        }
    }

    #[test]
    fn hash_collision_cannot_create_an_anchor() {
        let gaps = vec![Gap {
            scaffold: 0,
            start: 5,
            end: 7,
            anchors: Some([b"AACGT".to_vec(), b"TTGCA".to_vec()]),
            rejection: None,
        }];
        let mut index = Index::new(&gaps, 5, 3);
        // Deliberately put every anchor under the query's hash and bypass the
        // prefilter. Exact verification must still reject these false positives.
        let mut buf = Buffers::default();
        buf.packed.push_ascii(b"CCCCC");
        NtHasher::<false>::new(5)
            .hash_kmers_simd(buf.packed.as_slice(), 1)
            .collect_into(&mut buf.hashes);
        index
            .kmers
            .insert(buf.hashes[0], (0..index.anchors.len()).collect());
        let builder = minimizers(3, 3);
        builder.run(buf.packed.as_slice(), &mut buf.mini_positions);
        for &p in &buf.mini_positions {
            index
                .minis
                .entry(buf.packed.slice(p as usize..p as usize + 3).as_u64() as u32)
                .or_default()
                .extend(0..index.anchors.len());
        }
        assert!(
            index
                .extract(b"CCCCCCCCCCCC", &mut Scanner::default(), 100)
                .is_empty()
        );
    }
}
