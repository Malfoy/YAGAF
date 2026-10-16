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

/// Map regular, uncompressed files; codecs remain on the streamed reader path.
pub fn mmap_uncompressed(path: &Path) -> Result<Option<MmapInput<'static>>> {
    let metadata = std::fs::metadata(path)?;
    if !metadata.is_file() || metadata.len() == 0 {
        return Ok(None);
    }
    let mut input = FileInput::new(path)?;
    if input.is_compressed()? {
        Ok(None)
    } else {
        Ok(Some(MmapInput::new(path)?))
    }
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
            sequence: parser.get_dna_string_owned(),
        });
    }
    ensure!(!records.is_empty(), "No scaffold records");
    Ok(records)
}

/// Parse in-memory FASTA through the same helicase configuration as file input.
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
