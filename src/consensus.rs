use crate::is_base;
use anyhow::{Result, ensure};
use clap::ValueEnum;
use poasta::aligner::PoastaAligner;
use poasta::aligner::config::AffineMinGapCost;
use poasta::aligner::scoring::{AlignmentType, GapAffine};
use poasta::graphs::poa::POAGraph;
use poasta::io::fasta::poa_graph_to_msa;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MsaBackend {
    Auto,
    Poasta,
    Spoars,
    SpoarsScalar,
}

impl MsaBackend {
    pub fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Poasta => "poasta",
            Self::Spoars => "spoars",
            Self::SpoarsScalar => "spoars-scalar",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MsaOrder {
    Lexical,
    LongestFirst,
    MedianFirst,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Preset {
    Hifi,
    Ont,
}

impl Preset {
    /// Initial scoring presets, not calibrated acceptance thresholds.
    pub fn costs(self) -> GapAffine {
        match self {
            Self::Hifi => GapAffine::new(4, 2, 6),
            Self::Ont => GapAffine::new(4, 1, 4),
        }
    }
}

#[derive(Debug)]
pub struct Consensus {
    pub sequence: Vec<u8>,
    pub unanimous: usize,
    pub columns: usize,
    pub quality: f64,
    pub rejection: Option<&'static str>,
    pub backend: Option<MsaBackend>,
}

pub fn from_rows(rows: &[Vec<u8>]) -> Result<Consensus> {
    ensure!(!rows.is_empty(), "MSA has no sequences");
    let width = rows[0].len();
    ensure!(
        rows.iter().all(|r| r.len() == width),
        "Unequal MSA row lengths"
    );
    let mut consensus = Vec::new();
    let mut unanimous = 0;
    let mut columns = 0;
    let mut rejection = None;
    for col in 0..width {
        let mut counts = [0usize; 6]; // A C G T gap ambiguous
        for row in rows {
            let index = match row[col].to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'-' => 4,
                _ => 5,
            };
            counts[index] += 1;
        }
        if counts[4] == rows.len() {
            continue;
        }
        columns += 1;
        if counts[..4].contains(&rows.len()) {
            unanimous += 1;
        }
        let best = *counts.iter().max().unwrap();
        if counts.iter().filter(|&&count| count == best).count() != 1 {
            rejection = Some("consensus_tie");
            continue;
        }
        let winner = counts.iter().position(|&count| count == best).unwrap();
        match winner {
            0..=3 => consensus.push(b"ACGT"[winner]),
            4 => (),
            _ => rejection = Some("ambiguous_consensus"),
        }
    }
    if rejection.is_none() && (columns == 0 || consensus.is_empty()) {
        rejection = Some("empty_consensus");
    }
    Ok(Consensus {
        sequence: consensus,
        unanimous,
        columns,
        quality: if columns == 0 {
            0.0
        } else {
            100.0 * unanimous as f64 / columns as f64
        },
        rejection,
        backend: None,
    })
}

pub fn build(sequences: &mut [Vec<u8>], preset: Preset) -> Result<Consensus> {
    build_with_backend(sequences, preset, MsaBackend::Poasta)
}

pub fn build_with_backend<S: AsRef<[u8]>>(
    sequences: &mut [S],
    preset: Preset,
    backend: MsaBackend,
) -> Result<Consensus> {
    build_with_strategy(sequences, preset, backend, MsaOrder::Lexical)
}

pub fn build_with_strategy<S: AsRef<[u8]>>(
    sequences: &mut [S],
    preset: Preset,
    backend: MsaBackend,
    order: MsaOrder,
) -> Result<Consensus> {
    ensure!(!sequences.is_empty(), "No consensus evidence");
    ensure!(
        sequences
            .iter()
            .all(|s| !s.as_ref().is_empty() && s.as_ref().iter().all(|&b| is_base(b))),
        "Invalid consensus evidence"
    );
    if sequences
        .iter()
        .all(|s| s.as_ref() == sequences[0].as_ref())
    {
        let len = sequences[0].as_ref().len();
        return Ok(Consensus {
            sequence: sequences[0].as_ref().to_vec(),
            unanimous: len,
            columns: len,
            quality: 100.0,
            rejection: None,
            backend: None,
        });
    }
    let backend = if backend == MsaBackend::Auto {
        automatic_backend(sequences)
    } else {
        backend
    };
    // Stable lexical order makes progressive POA independent of worker scheduling.
    sequences.sort_unstable_by(|a, b| match order {
        MsaOrder::LongestFirst => b
            .as_ref()
            .len()
            .cmp(&a.as_ref().len())
            .then_with(|| a.as_ref().cmp(b.as_ref())),
        _ => a.as_ref().cmp(b.as_ref()),
    });
    if order == MsaOrder::MedianFirst {
        let mut lengths: Vec<_> = sequences.iter().map(|s| s.as_ref().len()).collect();
        let middle = lengths.len() / 2;
        let median = *lengths.select_nth_unstable(middle).1;
        let first = sequences
            .iter()
            .position(|s| s.as_ref().len() == median)
            .unwrap();
        sequences[..=first].rotate_right(1);
    }
    let rows = match backend {
        MsaBackend::Auto => unreachable!("Automatic backend must be resolved"),
        MsaBackend::Poasta => poasta_rows(sequences, preset)?,
        MsaBackend::Spoars => spoars_rows(sequences, preset, true)?,
        MsaBackend::SpoarsScalar => spoars_rows(sequences, preset, false)?,
    };
    validate_rows(&rows, sequences)?;
    let mut consensus = from_rows(&rows)?;
    consensus.backend = Some(backend);
    Ok(consensus)
}

/// Route by exact short-word sharing, never by a preferred consensus or score.
/// Every sequence remains in the eventual alignment with its original weight.
fn automatic_backend<S: AsRef<[u8]>>(sequences: &[S]) -> MsaBackend {
    use rustc_hash::{FxHashMap, FxHashSet};
    let mut frequencies: FxHashMap<&[u8], usize> = FxHashMap::default();
    let mut lengths: Vec<_> = sequences.iter().map(|s| s.as_ref().len()).collect();
    let middle = lengths.len() / 2;
    let median = *lengths.select_nth_unstable(middle).1;
    for s in sequences {
        *frequencies.entry(s.as_ref()).or_default() += 1;
    }
    let reference = frequencies
        .into_iter()
        .min_by(|(a, ac), (b, bc)| {
            bc.cmp(ac)
                .then_with(|| a.len().abs_diff(median).cmp(&b.len().abs_diff(median)))
                .then_with(|| a.cmp(b))
        })
        .unwrap()
        .0;
    let k = reference.len().min(15).min(*lengths.iter().min().unwrap());
    let word = |s: &[u8]| {
        s.iter().fold(0u32, |key, &b| {
            (key << 2)
                | match b.to_ascii_uppercase() {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => unreachable!(),
                }
        })
    };
    let words: FxHashSet<_> = reference.windows(k).map(word).collect();
    let (mut total, mut unmatched) = (0usize, 0usize);
    for sequence in sequences {
        let sequence = sequence.as_ref();
        let count = sequence.len() - k + 1;
        for start in (0..count).step_by(count.div_ceil(64)) {
            total += 1;
            unmatched += usize::from(!words.contains(&word(&sequence[start..start + k])));
        }
    }
    if unmatched > total / 10 {
        MsaBackend::Spoars
    } else {
        MsaBackend::Poasta
    }
}

fn poasta_rows<S: AsRef<[u8]>>(sequences: &[S], preset: Preset) -> Result<Vec<Vec<u8>>> {
    let aligner = PoastaAligner::new(AffineMinGapCost(preset.costs()), AlignmentType::Global);
    let mut graph = POAGraph::<u32>::new();
    for (id, sequence) in sequences.iter().enumerate() {
        let sequence = sequence.as_ref();
        let weights = vec![1; sequence.len()];
        if id == 0 {
            graph.add_alignment_with_weights(&id.to_string(), sequence, None, &weights)?;
        } else {
            let result = aligner.align::<u32, _>(&graph, sequence);
            graph.add_alignment_with_weights(
                &id.to_string(),
                sequence,
                Some(&result.alignment),
                &weights,
            )?;
        }
    }
    Ok(poa_graph_to_msa(&graph)?)
}

fn spoars_rows<S: AsRef<[u8]>>(
    sequences: &[S],
    preset: Preset,
    simd: bool,
) -> Result<Vec<Vec<u8>>> {
    use poasta::aligner::scoring::AlignmentCosts;
    use spoars::align::{AlignmentEngine, AlignmentType as Mode, Scoring, SimdEngine, SisdEngine};
    let costs = preset.costs();
    let extend = costs.gap_extend() as i8;
    // POASTA charges open + length*extend; SPOA charges g + (length-1)*e.
    let open = costs.gap_open() as i8 + extend;
    let scoring = Scoring::new(0, -(costs.mismatch() as i8), -open, -extend, -open, -extend)?;
    let mut graph = spoars::graph::Graph::new();
    let mut add = |engine: &mut dyn AlignmentEngine| -> Result<()> {
        for sequence in sequences {
            let sequence = sequence.as_ref();
            let (alignment, _) = engine.align(sequence, &graph);
            graph.add_alignment_weight(&alignment, sequence, 1)?;
        }
        Ok(())
    };
    if simd {
        add(&mut SimdEngine::new(Mode::Global, scoring))?;
    } else {
        add(&mut SisdEngine::new(Mode::Global, scoring))?;
    }
    Ok(graph
        .generate_msa(false)
        .into_iter()
        .map(String::into_bytes)
        .collect())
}

fn validate_rows<S: AsRef<[u8]>>(rows: &[Vec<u8>], sequences: &[S]) -> Result<()> {
    ensure!(rows.len() == sequences.len(), "MSA is missing sequences");
    // Each backend must preserve every contributing sequence and its gaps.
    for (row, original) in rows.iter().zip(sequences.iter()) {
        ensure!(
            row.iter()
                .filter(|&&b| b != b'-')
                .copied()
                .eq(original.as_ref().iter().copied()),
            "MSA did not preserve input sequence"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn automatic_engine_selection_uses_all_evidence_and_ignores_input_order() {
        let mut state = 818725439263u64;
        let truth: Vec<_> = (0..1000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[state as usize & 3]
            })
            .collect();
        for (spacing, expected) in [(1000, MsaBackend::Poasta), (50, MsaBackend::Spoars)] {
            let mut sequences = Vec::new();
            for read in 0..12 {
                let mut s = truth.clone();
                for pos in (read..s.len()).step_by(spacing) {
                    s[pos] = if s[pos] == b'A' { b'C' } else { b'A' };
                }
                sequences.push(s);
            }
            assert_eq!(automatic_backend(&sequences), expected);
            sequences.reverse();
            assert_eq!(automatic_backend(&sequences), expected);
            let mut forced = sequences.clone();
            let actual =
                build_with_backend(&mut sequences, Preset::Hifi, MsaBackend::Auto).unwrap();
            let reference = build_with_backend(&mut forced, Preset::Hifi, expected).unwrap();
            assert_eq!(actual.backend, Some(expected));
            assert_eq!(actual.sequence, truth);
            assert_eq!(actual.sequence, reference.sequence);
            assert_eq!(
                (actual.unanimous, actual.columns),
                (reference.unanimous, reference.columns)
            );
        }
    }

    #[test]
    fn inline_and_heap_evidence_produce_the_same_msa_metrics() {
        for length in [16, 17, 49] {
            let truth: Vec<_> = (0..length).map(|i| b"ACGT"[i % 4]).collect();
            let mut insertion = truth.clone();
            insertion.insert(length / 2, b'A');
            for backend in [
                MsaBackend::Poasta,
                MsaBackend::Spoars,
                MsaBackend::SpoarsScalar,
            ] {
                for preset in [Preset::Hifi, Preset::Ont] {
                    let mut heap = vec![truth.clone(), insertion.clone(), truth.clone()];
                    let mut inline: Vec<_> = heap
                        .iter()
                        .map(|s| crate::index::Sequence::from_slice(s))
                        .collect();
                    let expected = build_with_backend(&mut heap, preset, backend).unwrap();
                    let actual = build_with_backend(&mut inline, preset, backend).unwrap();
                    assert_eq!(actual.sequence, expected.sequence);
                    assert_eq!(
                        (
                            actual.unanimous,
                            actual.columns,
                            actual.quality,
                            actual.rejection
                        ),
                        (
                            expected.unanimous,
                            expected.columns,
                            expected.quality,
                            expected.rejection
                        )
                    );
                }
            }
        }
    }

    #[test]
    fn alternate_backends_preserve_spans_and_the_quality_definition() {
        for backend in [
            MsaBackend::Poasta,
            MsaBackend::Spoars,
            MsaBackend::SpoarsScalar,
        ] {
            for preset in [Preset::Hifi, Preset::Ont] {
                let mut reads = vec![b"ACGTAC".to_vec(), b"ACTGTAC".to_vec(), b"ACTGTAC".to_vec()];
                let c = build_with_backend(&mut reads, preset, backend).unwrap();
                assert_eq!(c.sequence, b"ACTGTAC");
                assert_eq!((c.unanimous, c.columns), (6, 7));
                let mut reads = vec![b"A".to_vec(), b"C".to_vec()];
                assert_eq!(
                    build_with_backend(&mut reads, preset, backend)
                        .unwrap()
                        .rejection,
                    Some("consensus_tie")
                );
                let mut reads = vec![b"A".to_vec(), b"CA".to_vec(), b"CA".to_vec()];
                let c = build_with_backend(&mut reads, preset, backend).unwrap();
                assert_eq!(c.sequence, b"CA");
                assert_eq!((c.unanimous, c.columns), (1, 2));
            }
        }
    }
    #[test]
    fn unanimity_includes_indels_and_ambiguity() {
        let c = from_rows(&[b"AC-G-N".to_vec(), b"ACTG-A".to_vec(), b"ACTG-A".to_vec()]).unwrap();
        assert_eq!((c.unanimous, c.columns, c.quality), (3, 5, 60.0));
        assert_eq!(c.sequence, b"ACTGA");
    }
    #[test]
    fn poasta_aligns_insertions_and_substitutions() {
        let mut evidence = vec![b"ACGTAC".to_vec(), b"ACTGTAC".to_vec(), b"ACTGTAC".to_vec()];
        let c = build(&mut evidence, Preset::Hifi).unwrap();
        assert_eq!(c.sequence, b"ACTGTAC");
        assert_eq!((c.unanimous, c.columns), (6, 7));
        let mut evidence = vec![b"ACGTAC".to_vec(), b"ACGTTC".to_vec(), b"ACGTTC".to_vec()];
        let c = build(&mut evidence, Preset::Ont).unwrap();
        assert_eq!(c.sequence, b"ACGTTC");
        assert_eq!((c.unanimous, c.columns), (5, 6));
    }
    #[test]
    fn tied_and_empty_consensuses_are_rejected() {
        assert_eq!(
            from_rows(&[b"A".to_vec(), b"C".to_vec()])
                .unwrap()
                .rejection,
            Some("consensus_tie")
        );
        assert_eq!(
            from_rows(&[b"--".to_vec(), b"--".to_vec()])
                .unwrap()
                .rejection,
            Some("empty_consensus")
        );
    }

    #[test]
    fn leading_and_trailing_indels_keep_msa_rows_rectangular() {
        for reads in [
            vec![b"A".to_vec(), b"CA".to_vec(), b"CA".to_vec()],
            vec![b"A".to_vec(), b"AC".to_vec(), b"AC".to_vec()],
            vec![b"A".to_vec(), b"C".to_vec(), b"A".to_vec()],
        ] {
            let c = build(&mut reads.clone(), Preset::Hifi).unwrap();
            assert!(c.columns > 0);
        }
    }
}
