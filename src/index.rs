// Adapted from Malfoy/K2Rmini src/bin/k2rminimulti.rs (GPL-3.0).
// Pinned upstream revision and adaptation details are in THIRD_PARTY.md.
use crate::{Gap, Scaffold, is_base, revcomp};
use rustc_hash::{FxHashMap, FxHashSet};
use simd_minimizers::minimizers;
use simd_minimizers::packed_seq::{PackedSeqVec, Seq, SeqVec};
use simd_minimizers::seq_hash::{KmerHasher, NtHasher};

struct Anchor {
    gap: usize,
    side: usize,
    reverse: bool,
    bases: Vec<u8>,
}

pub struct Index {
    k: usize,
    m: usize,
    anchors: Vec<Anchor>,
    minis: FxHashSet<u32>,
    kmers: FxHashMap<u32, Vec<usize>>,
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
            minis: FxHashSet::default(),
            kmers: FxHashMap::default(),
        };
        let builder = minimizers(m, k - m + 1);
        let hasher = NtHasher::<false>::new(k);
        let mut buffers = Buffers::default();
        for (gap, item) in gaps
            .iter()
            .enumerate()
            .filter(|(_, g)| g.rejection.is_none())
        {
            if let Some(anchors) = &item.anchors {
                for (side, bases) in anchors.iter().enumerate() {
                    for reverse in [false, true] {
                        let bases = if reverse {
                            revcomp(bases)
                        } else {
                            bases.clone()
                        };
                        buffers.packed.clear();
                        buffers.packed.push_ascii(&bases);
                        buffers.mini_positions.clear();
                        builder.run(buffers.packed.as_slice(), &mut buffers.mini_positions);
                        for &pos in &buffers.mini_positions {
                            index.minis.insert(
                                buffers
                                    .packed
                                    .slice(pos as usize..pos as usize + m)
                                    .as_u64() as u32,
                            );
                        }
                        buffers.hashes.clear();
                        hasher
                            .hash_kmers_simd(buffers.packed.as_slice(), 1)
                            .collect_into(&mut buffers.hashes);
                        let id = index.anchors.len();
                        index.kmers.entry(buffers.hashes[0]).or_default().push(id);
                        index.anchors.push(Anchor {
                            gap,
                            side,
                            reverse,
                            bases,
                        });
                    }
                }
            }
        }
        index
    }

    // Each ACTG run is independent. Never join across an ambiguous base: that
    // would introduce nonexistent k-mers and corrupt read coordinates.
    fn visit(&self, sequence: &[u8], buffers: &mut Buffers, mut found: impl FnMut(usize, usize)) {
        let builder = minimizers(self.m, self.k - self.m + 1);
        let hasher = NtHasher::<false>::new(self.k);
        let mut offset = 0;
        for run in sequence.split(|&base| !is_base(base)) {
            if run.len() >= self.k && !self.kmers.is_empty() {
                buffers.packed.clear();
                buffers.packed.push_ascii(run);
                buffers.mini_positions.clear();
                builder.run(buffers.packed.as_slice(), &mut buffers.mini_positions);
                let candidate = buffers.mini_positions.iter().any(|&pos| {
                    let key = buffers
                        .packed
                        .slice(pos as usize..pos as usize + self.m)
                        .as_u64() as u32;
                    self.minis.contains(&key)
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
        let mut counts = vec![0u8; self.anchors.len()];
        let mut buffers = Buffers::default();
        for scaffold in scaffolds {
            self.visit(&scaffold.sequence, &mut buffers, |id, _| {
                counts[id] = counts[id].saturating_add(1).min(2)
            });
        }
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
        let mut result: FxHashMap<usize, Span> = FxHashMap::default();
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
                result.insert(gap, Span::AmbiguousPair);
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
            result
                .entry(gap)
                .and_modify(|old| *old = Span::AmbiguousPair)
                .or_insert(span);
        }
        let mut result: Vec<_> = result.into_iter().collect();
        result.sort_unstable_by_key(|(gap, _)| *gap);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
                .insert(buf.packed.slice(p as usize..p as usize + 3).as_u64() as u32);
        }
        assert!(
            index
                .extract(b"CCCCCCCCCCCC", &mut Scanner::default(), 100)
                .is_empty()
        );
    }
}
