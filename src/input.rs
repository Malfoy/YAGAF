use crate::Scaffold;
use anyhow::{Context, Result, ensure};
use helicase::input::*;
use helicase::parser::Format;
use helicase::*;
use std::path::Path;

// Keep every base, including Ns and MSA '-' characters, without changing coordinates.
pub const CONFIG: Config = ParserOptions::default()
    .dna_string()
    .keep_non_actg()
    .compute_quality()
    .config();

pub fn open(path: &Path) -> Result<FastxParser<'static, CONFIG>> {
    FastxParser::<CONFIG>::from_file(path)
        .with_context(|| format!("Cannot open FASTA/FASTQ {}", path.display()))
}

pub fn scaffolds(path: &Path) -> Result<Vec<Scaffold>> {
    let mut parser = open(path)?;
    ensure!(
        matches!(parser.format(), Format::Fasta),
        "Scaffolds must be FASTA"
    );
    let mut records = Vec::new();
    while parser.next().is_some() {
        ensure!(!parser.get_header().is_empty(), "Empty scaffold header");
        ensure!(
            parser.get_dna_string().len() <= u32::MAX as usize / 2,
            "Scaffold exceeds the current 2^31-1 base indexing limit"
        );
        records.push(Scaffold {
            header: parser.get_header().to_vec(),
            sequence: parser.get_dna_string().to_vec(),
        });
    }
    ensure!(!records.is_empty(), "No scaffold records");
    Ok(records)
}

/// Also used for POASTA's in-memory MSA export. All FASTA parsing stays in helicase.
pub fn fasta_bytes(data: &[u8]) -> Result<Vec<Scaffold>> {
    ensure!(!data.is_empty(), "Empty FASTA buffer");
    let mut parser = FastxParser::<CONFIG>::from_slice(data)?;
    ensure!(matches!(parser.format(), Format::Fasta), "Expected FASTA");
    let mut records = Vec::new();
    while parser.next().is_some() {
        records.push(Scaffold {
            header: parser.get_header().to_vec(),
            sequence: parser.get_dna_string().to_vec(),
        });
    }
    Ok(records)
}
