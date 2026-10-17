// Adapted from Malfoy/K2Rmini src/bin/k2rminimulti.rs (GPL-3.0).
// Pinned upstream revision and adaptation details are in THIRD_PARTY.md.
use crate::flank::{self, AlignmentScratch, Pattern, ReadView};
use crate::{Gap, Scaffold, is_base, revcomp};
use clap::ValueEnum;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use simd_minimizers::minimizers;
use simd_minimizers::packed_seq::{PackedSeqVec, Seq, SeqVec};
use simd_minimizers::seq_hash::{KmerHasher, NtHasher};
use smallvec::SmallVec;

pub type Sequence = SmallVec<[u8; 16]>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Recruitment {
    Strict,
    BestContext,
}

struct Anchor {
    gap: usize,
    side: usize,
    reverse: bool,
    bases: Vec<u8>,
    minimizer_offset: usize,
    heavy_minimizer: bool,
    minimizer_key: u32,
    projection: Pattern,
}

pub struct Index {
    k: usize,
    m: usize,
    anchors: Vec<Anchor>,
    minis: FxHashMap<u32, Vec<usize>>,
    kmers: FxHashMap<u32, Vec<usize>>,
    use_offsets: bool,
    flanks: Vec<[Vec<u8>; 2]>,
    validation: Vec<[Pattern; 2]>,
    recruitment: Recruitment,
    recruitment_margin: usize,
    fallback: Option<(usize, Recruitment, usize)>,
}

#[derive(Default)]
pub struct Buffers {
    packed: PackedSeqVec,
    mini_positions: Vec<u32>,
    hashes: Vec<u32>,
}

#[derive(Default)]
pub struct Scanner {
    buffers: Buffers,
    masks: Vec<u8>,
    active: Vec<usize>,
    placements: Vec<Placement>,
    alignment: AlignmentScratch,
}

struct Placement {
    seen: [FxHashSet<usize>; 2],
    boundaries: [Vec<(usize, usize)>; 2],
    interval: Option<(usize, usize)>,
    best_cost: usize,
    runner_up: usize,
    terminal: bool,
}
impl Default for Placement {
    fn default() -> Self {
        Self {
            seen: Default::default(),
            boundaries: Default::default(),
            interval: None,
            best_cost: usize::MAX,
            runner_up: usize::MAX,
            terminal: false,
        }
    }
}
impl Placement {
    fn clear(&mut self) {
        for seen in &mut self.seen {
            seen.clear();
        }
        for boundaries in &mut self.boundaries {
            boundaries.clear();
        }
        self.interval = None;
        self.best_cost = usize::MAX;
        self.runner_up = usize::MAX;
        self.terminal = false;
    }
    fn add(
        &mut self,
        side: usize,
        boundary: usize,
        score: usize,
        recruitment: Recruitment,
        margin: usize,
    ) {
        for &(other, other_score) in &self.boundaries[1 - side] {
            let pair = if side == 0 {
                (boundary, other)
            } else {
                (other, boundary)
            };
            if pair.0 > pair.1 {
                continue;
            }
            if recruitment == Recruitment::Strict {
                if self.interval.is_some_and(|old| old != pair) {
                    self.terminal = true;
                    return;
                }
                self.interval = Some(pair);
            } else {
                let cost = score + other_score;
                if cost < self.best_cost {
                    self.runner_up = self.best_cost;
                    self.best_cost = cost;
                    self.interval = Some(pair);
                } else {
                    self.runner_up = self.runner_up.min(cost);
                }
                // Even a future zero-cost placement cannot establish this margin
                // against the already observed second-best candidate.
                if self.runner_up < margin {
                    self.terminal = true;
                    return;
                }
            }
        }
        self.boundaries[side].push((boundary, score));
    }
}

#[derive(Debug)]
pub enum Span {
    Sequence(Sequence),
    AmbiguousPair,
    /// Both seed sides were observed but no ordered, validated interval survived.
    Unplaced,
    Rejected(&'static str),
}

impl Index {
    /// Combine two query tiers without rescanning the read input. Virtual gap
    /// identifiers of the fallback index follow those of the primary index.
    pub fn with_fallback(mut self, mut other: Self) -> Self {
        assert!(
            self.k == other.k
                && self.m == other.m
                && self.flanks.len() == other.flanks.len()
                && self.fallback.is_none()
                && other.fallback.is_none()
        );
        let gap_offset = self.flanks.len();
        let anchor_offset = self.anchors.len();
        self.fallback = Some((gap_offset, other.recruitment, other.recruitment_margin));
        for anchor in &mut other.anchors {
            anchor.gap += gap_offset;
        }
        self.anchors.append(&mut other.anchors);
        self.flanks.append(&mut other.flanks);
        self.validation.append(&mut other.validation);
        for (key, ids) in other.minis {
            self.minis
                .entry(key)
                .or_default()
                .extend(ids.into_iter().map(|id| id + anchor_offset));
        }
        for (key, ids) in other.kmers {
            self.kmers
                .entry(key)
                .or_default()
                .extend(ids.into_iter().map(|id| id + anchor_offset));
        }
        for anchor in &mut self.anchors {
            anchor.heavy_minimizer = self.minis[&anchor.minimizer_key].len() > 64;
        }
        self
    }
    fn pairing(&self, gap: usize) -> (Recruitment, usize) {
        if let Some((split, policy, margin)) = self.fallback
            && gap >= split
        {
            (policy, margin)
        } else {
            (self.recruitment, self.recruitment_margin)
        }
    }
    pub fn with_flank_edit_limit(mut self, limit: Option<usize>) -> Self {
        if let Some(limit) = limit {
            for anchor in &mut self.anchors {
                anchor.projection.limit_edits(limit);
            }
            for pair in &mut self.validation {
                for pattern in pair {
                    pattern.limit_edits(limit);
                }
            }
        }
        self
    }
    pub fn with_recruitment(mut self, recruitment: Recruitment, margin: usize) -> Self {
        assert!(margin > 0);
        self.recruitment = recruitment;
        self.recruitment_margin = margin;
        self
    }
    pub fn new(gaps: &[Gap], k: usize, m: usize) -> Self {
        let flanks = gaps
            .iter()
            .map(|gap| gap.anchors.clone().unwrap_or_default())
            .collect();
        Self::build(gaps, flanks, k, m, 1)
    }

    pub fn with_windows(
        scaffolds: &[Scaffold],
        gaps: &[Gap],
        k: usize,
        m: usize,
        window: usize,
        step: usize,
    ) -> Self {
        let flanks = gaps
            .iter()
            .map(|gap| {
                let seq = &scaffolds[gap.scaffold].sequence;
                let mut left = gap.start;
                while left > 0 && gap.start - left < window && is_base(seq[left - 1]) {
                    left -= 1;
                }
                let mut right = gap.end;
                while right < seq.len() && right - gap.end < window && is_base(seq[right]) {
                    right += 1;
                }
                [
                    seq[left..gap.start].to_ascii_uppercase(),
                    seq[gap.end..right].to_ascii_uppercase(),
                ]
            })
            .collect();
        Self::build(gaps, flanks, k, m, step)
    }

    fn build(gaps: &[Gap], flanks: Vec<[Vec<u8>; 2]>, k: usize, m: usize, step: usize) -> Self {
        let mut index = Self {
            k,
            m,
            recruitment: Recruitment::Strict,
            recruitment_margin: 2,
            fallback: None,
            anchors: Vec::new(),
            minis: FxHashMap::default(),
            kmers: FxHashMap::default(),
            use_offsets: m >= 9,
            validation: flanks
                .iter()
                .map(|s| [Pattern::new(&s[0], true), Pattern::new(&s[1], false)])
                .collect(),
            flanks,
        };
        let builder = minimizers(m, k - m + 1);
        let hasher = NtHasher::<false>::new(k);
        let mut jobs = Vec::new();
        for (gap, _) in gaps
            .iter()
            .enumerate()
            .filter(|(_, g)| g.rejection.is_none())
        {
            for side in 0..2 {
                let context = &index.flanks[gap][side];
                if context.len() < k {
                    continue;
                }
                let farthest = context.len() - k;
                let mut distances: Vec<_> = (0..=farthest).step_by(step).collect();
                if distances.last() != Some(&farthest) {
                    distances.push(farthest);
                }
                for distance in distances {
                    let start = if side == 0 {
                        context.len() - k - distance
                    } else {
                        distance
                    };
                    for reverse in [false, true] {
                        let projection = Pattern::new(
                            if side == 0 {
                                &context[context.len() - distance..]
                            } else {
                                &context[..distance]
                            },
                            side == 1,
                        );
                        jobs.push((gap, side, reverse, &context[start..start + k], projection));
                    }
                }
            }
        }
        let prepared: Vec<_> = jobs
            .into_par_iter()
            .with_min_len(64)
            .map_init(
                Buffers::default,
                |buffers, (gap, side, reverse, bases, projection)| {
                    let bases = if reverse {
                        revcomp(bases)
                    } else {
                        bases.to_vec()
                    };
                    buffers.packed.clear();
                    buffers.packed.push_ascii(&bases);
                    buffers.mini_positions.clear();
                    builder.run(buffers.packed.as_slice(), &mut buffers.mini_positions);
                    let minimizer_offset = buffers.mini_positions[0] as usize;
                    let key = buffers
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
                            heavy_minimizer: false,
                            minimizer_key: key,
                            projection,
                        },
                        key,
                        buffers.hashes[0],
                    )
                },
            )
            .collect();
        index.anchors.reserve(prepared.len());
        index.minis.reserve(prepared.len());
        index.kmers.reserve(prepared.len());
        for (anchor, key, hash) in prepared {
            let id = index.anchors.len();
            index.minis.entry(key).or_default().push(id);
            index.kmers.entry(hash).or_default().push(id);
            index.anchors.push(anchor);
        }
        for anchor in &mut index.anchors {
            anchor.heavy_minimizer = index.minis[&anchor.minimizer_key].len() > 64;
        }
        index
    }

    // Each ACTG run is independent. Never join across an ambiguous base: that
    // would introduce nonexistent k-mers and corrupt read coordinates.
    fn candidates(
        &self,
        sequence: &[u8],
        buffers: &mut Buffers,
        mut found: impl FnMut(usize, usize),
    ) {
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
                    let mut needs_hashing = false;
                    for &pos in &buffers.mini_positions {
                        let key = buffers
                            .packed
                            .slice(pos as usize..pos as usize + self.m)
                            .as_u64() as u32;
                        if let Some(ids) = self.minis.get(&key) {
                            if ids.len() > 64 {
                                needs_hashing = true;
                                continue;
                            }
                            for &id in ids {
                                let anchor = &self.anchors[id];
                                if let Some(start) =
                                    (pos as usize).checked_sub(anchor.minimizer_offset)
                                    && start + self.k <= run.len()
                                {
                                    found(id, offset + start);
                                }
                            }
                        }
                    }
                    if needs_hashing {
                        buffers.hashes.clear();
                        hasher
                            .hash_kmers_simd(buffers.packed.as_slice(), 1)
                            .collect_into(&mut buffers.hashes);
                        for (pos, hash) in buffers.hashes.iter().enumerate() {
                            if let Some(ids) = self.kmers.get(hash) {
                                for &id in ids {
                                    let anchor = &self.anchors[id];
                                    if anchor.heavy_minimizer {
                                        found(id, offset + pos);
                                    }
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
                                found(id, offset + pos);
                            }
                        }
                    }
                }
            }
            offset += run.len() + 1;
        }
    }

    #[cfg(test)]
    fn visit(&self, sequence: &[u8], buffers: &mut Buffers, mut found: impl FnMut(usize, usize)) {
        self.candidates(sequence, buffers, |id, pos| {
            if sequence[pos..pos + self.k].eq_ignore_ascii_case(&self.anchors[id].bases) {
                found(id, pos);
            }
        });
    }

    pub fn extract(
        &self,
        sequence: &[u8],
        scanner: &mut Scanner,
        max_span: usize,
    ) -> Vec<(usize, Span)> {
        let mut spans = self.extract_audited(sequence, scanner, max_span);
        spans.retain(|(_, span)| !matches!(span, Span::Unplaced));
        spans
    }

    /// Include paired-seed candidates rejected by context or ordering checks.
    pub fn extract_audited(
        &self,
        sequence: &[u8],
        scanner: &mut Scanner,
        max_span: usize,
    ) -> Vec<(usize, Span)> {
        for key in scanner.active.drain(..) {
            scanner.masks[key] = 0;
            scanner.placements[key].clear();
        }
        scanner.masks.resize(self.flanks.len() * 2, 0);
        scanner
            .placements
            .resize_with(self.flanks.len() * 2, Placement::default);
        // The first search records only exact side presence. A second in-memory
        // search is needed only for reads with a potential paired gap. No seed-hit
        // lists or per-read sorting are needed, even for low-complexity repeats.
        let masks = &mut scanner.masks;
        let active = &mut scanner.active;
        self.candidates(sequence, &mut scanner.buffers, |id, pos| {
            let anchor = &self.anchors[id];
            let key = anchor.gap * 2 + usize::from(anchor.reverse);
            let side = 1 << anchor.side;
            if masks[key] & side != 0 {
                return;
            }
            if !sequence[pos..pos + self.k].eq_ignore_ascii_case(&anchor.bases) {
                return;
            }
            if masks[key] == 0 {
                active.push(key);
            }
            masks[key] |= side;
        });
        if !scanner.active.iter().any(|&key| scanner.masks[key] == 3) {
            return Vec::new();
        }
        let placements = &mut scanner.placements;
        let alignment = &mut scanner.alignment;
        self.candidates(sequence, &mut scanner.buffers, |id, pos| {
            let anchor = &self.anchors[id];
            let key = anchor.gap * 2 + usize::from(anchor.reverse);
            if scanner.masks[key] != 3 || placements[key].terminal {
                return;
            }
            if let Some((split, _, _)) = self.fallback
                && anchor.gap >= split
            {
                let primary = (anchor.gap - split) * 2;
                // Any primary placement will ultimately replace or block this
                // read's weaker evidence, regardless of later competing hits.
                if (primary..primary + 2)
                    .any(|id| placements[id].terminal || placements[id].interval.is_some())
                {
                    return;
                }
            }
            if !sequence[pos..pos + self.k].eq_ignore_ascii_case(&anchor.bases) {
                return;
            }
            let view = ReadView {
                bytes: sequence,
                reverse: anchor.reverse,
            };
            let pos = if anchor.reverse {
                sequence.len() - pos - self.k
            } else {
                pos
            };
            let Some(boundary) = flank::project(
                &anchor.projection,
                anchor.side,
                pos,
                self.k,
                &view,
                alignment,
            ) else {
                return;
            };
            let placement = &mut placements[key];
            if !placement.seen[anchor.side].insert(boundary) {
                return;
            }
            let Some(score) = flank::context_score(
                &self.flanks[anchor.gap][anchor.side],
                &self.validation[anchor.gap][anchor.side],
                anchor.side,
                boundary,
                self.k,
                &view,
                alignment,
            ) else {
                return;
            };
            let (policy, margin) = self.pairing(anchor.gap);
            placement.add(anchor.side, boundary, score, policy, margin);
        });
        let mut result = Vec::new();
        for &key in &scanner.active {
            let placement = &scanner.placements[key];
            let gap = key / 2;
            let (policy, margin) = self.pairing(gap);
            if placement.terminal
                || (policy == Recruitment::BestContext
                    && placement.interval.is_some()
                    && placement.runner_up.saturating_sub(placement.best_cost) < margin)
            {
                result.push((gap, Span::AmbiguousPair));
                continue;
            }
            let Some((start, end)) = placement.interval else {
                if scanner.masks[key] == 3 {
                    result.push((gap, Span::Unplaced));
                }
                continue;
            };
            let view = ReadView {
                bytes: sequence,
                reverse: key % 2 != 0,
            };
            let span = if end - start > max_span {
                Span::Rejected("span_limit")
            } else if start == end {
                Span::Rejected("empty_read_span")
            } else {
                view.sequence(start, end)
                    .map_or(Span::Rejected("ambiguous_read_span"), Span::Sequence)
            };
            result.push((gap, span));
        }
        result.sort_unstable_by_key(|(gap, _)| *gap);
        result.dedup_by(|later, first| {
            if later.0 == first.0 {
                if matches!(first.1, Span::Unplaced) {
                    std::mem::swap(&mut first.1, &mut later.1);
                } else if !matches!(later.1, Span::Unplaced) {
                    first.1 = Span::AmbiguousPair;
                }
                true
            } else {
                false
            }
        });
        if let Some((split, _, _)) = self.fallback {
            let pivot = result.partition_point(|(gap, _)| *gap < split);
            let (primary, rescue) = result.split_at_mut(pivot);
            for (gap, span) in rescue {
                if let Ok(id) = primary.binary_search_by_key(&(*gap - split), |(gap, _)| *gap) {
                    match &primary[id].1 {
                        // Primary reads are transferred only if this gap needs
                        // rescue; never count one read in both tiers.
                        Span::Sequence(_) => *span = Span::Unplaced,
                        Span::AmbiguousPair => *span = Span::AmbiguousPair,
                        Span::Rejected(reason) => *span = Span::Rejected(reason),
                        Span::Unplaced => (),
                    }
                }
            }
        }
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
    #[test]
    fn heavy_minimizer_hybrid_matches_exhaustive_hits() {
        let left = b"ACGTCAGCTAGGATCGATGCATCCGATAGCA".to_vec();
        let right = b"TGATGCCAGTCATGACCGTAGCTACGGTACC".to_vec();
        let mut gaps: Vec<_> = (0..80)
            .map(|_| Gap {
                scaffold: 0,
                start: 31,
                end: 32,
                anchors: Some([left.clone(), right.clone()]),
                rejection: None,
            })
            .collect();
        gaps.push(Gap {
            scaffold: 0,
            start: 31,
            end: 32,
            anchors: Some([b"GCTACAGAGTCTCGCAATAGGTTAGCCATGT".to_vec(), right.clone()]),
            rejection: None,
        });
        let index = Index::new(&gaps, 31, 15);
        assert!(index.anchors.iter().any(|a| a.heavy_minimizer));
        assert!(index.anchors.iter().any(|a| !a.heavy_minimizer));
        let read = [
            left.as_slice(),
            b"ACTG",
            right.as_slice(),
            b"N",
            index.anchors.last().unwrap().bases.as_slice(),
        ]
        .concat();
        let mut expected = Vec::new();
        for (id, anchor) in index.anchors.iter().enumerate() {
            for (pos, word) in read.windows(31).enumerate() {
                if word.eq_ignore_ascii_case(&anchor.bases) {
                    expected.push((id, pos));
                }
            }
        }
        let mut actual = Vec::new();
        index.visit(&read, &mut Buffers::default(), |id, pos| {
            actual.push((id, pos))
        });
        expected.sort_unstable();
        actual.sort_unstable();
        assert_eq!(actual, expected);
    }

    #[test]
    fn windows_stop_at_ambiguity_and_neighboring_gaps() {
        let scaffolds = vec![Scaffold {
            header: b"s".to_vec(),
            sequence: b"AAAAARACGTACGTNNACGTACGTNNACGTACGTRTTTT".to_vec(),
        }];
        let gaps = crate::discover(&scaffolds, 5);
        let index = Index::with_windows(&scaffolds, &gaps, 5, 3, 100, 2);
        assert_eq!(
            index.flanks[0],
            [b"ACGTACGT".to_vec(), b"ACGTACGT".to_vec()]
        );
        assert_eq!(
            index.flanks[1],
            [b"ACGTACGT".to_vec(), b"ACGTACGT".to_vec()]
        );
        assert_eq!(
            index
                .anchors
                .iter()
                .filter(|a| a.gap == 0 && a.side == 0 && !a.reverse)
                .map(|a| a.bases.clone())
                .collect::<Vec<_>>(),
            vec![b"TACGT".to_vec(), b"CGTAC".to_vec(), b"ACGTA".to_vec()]
        );
    }
}
