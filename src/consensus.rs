use crate::{input, is_base};
use anyhow::{Result, ensure};
use clap::ValueEnum;
use poasta::aligner::PoastaAligner;
use poasta::aligner::config::AffineMinGapCost;
use poasta::aligner::scoring::{AlignmentType, GapAffine};
use poasta::graphs::poa::POAGraph;
use poasta::io::fasta::poa_graph_to_fasta;

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
    })
}

pub fn build(sequences: &mut [Vec<u8>], preset: Preset) -> Result<Consensus> {
    ensure!(!sequences.is_empty(), "No consensus evidence");
    ensure!(
        sequences
            .iter()
            .all(|s| !s.is_empty() && s.iter().all(|&b| is_base(b))),
        "Invalid consensus evidence"
    );
    // Stable lexical order makes progressive POA independent of worker scheduling.
    sequences.sort_unstable();
    if sequences.iter().all(|s| s == &sequences[0]) {
        let len = sequences[0].len();
        return Ok(Consensus {
            sequence: sequences[0].clone(),
            unanimous: len,
            columns: len,
            quality: 100.0,
            rejection: None,
        });
    }
    let aligner = PoastaAligner::new(AffineMinGapCost(preset.costs()), AlignmentType::Global);
    let mut graph = POAGraph::<u32>::new();
    for (id, sequence) in sequences.iter().enumerate() {
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
    let mut fasta = Vec::new();
    poa_graph_to_fasta(&graph, &mut fasta)?;
    let rows = input::fasta_bytes(&fasta)?;
    ensure!(
        rows.len() == sequences.len(),
        "POASTA MSA is missing sequences"
    );
    // Check that serialization/parsing retained every base and MSA gap.
    for (row, original) in rows.iter().zip(sequences.iter()) {
        ensure!(
            row.sequence
                .iter()
                .filter(|&&b| b != b'-')
                .copied()
                .eq(original.iter().copied()),
            "MSA did not preserve input sequence"
        );
    }
    from_rows(&rows.into_iter().map(|r| r.sequence).collect::<Vec<_>>())
}

#[cfg(test)]
mod tests {
    use super::*;
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
