//! Mask contiguous reference intervals for evaluating gap recruitment and MSA.
//! Reference FASTA is read exclusively through helicase.
use anyhow::{Result, ensure};
use clap::Parser;
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::PathBuf,
};
use yagaf::{input, is_base};
#[derive(Parser)]
struct Options {
    #[arg(long)]
    reference: PathBuf,
    #[arg(long)]
    directory: PathBuf,
    #[arg(long, default_value_t = 200)]
    gaps: usize,
    #[arg(long, value_delimiter = ',', default_value = "50,100,250,500,1000")]
    lengths: Vec<usize>,
    #[arg(long, default_value_t = 42)]
    seed: u64,
}
fn main() -> Result<()> {
    let a = Options::parse();
    ensure!(
        a.seed > 0 && a.gaps > 0 && !a.lengths.is_empty() && a.lengths.iter().all(|&n| n > 0),
        "Positive seed, gap count, and lengths required"
    );
    let records = input::scaffolds(&a.reference)?;
    ensure!(records.len() == 1, "Expected one reference chromosome");
    let reference = &records[0];
    ensure!(
        reference.sequence.iter().all(|&b| is_base(b)),
        "Reference contains ambiguous bases"
    );
    let total: usize = (0..a.gaps).map(|i| a.lengths[i % a.lengths.len()]).sum();
    ensure!(
        total < reference.sequence.len() / 4,
        "Mask less than one quarter of the chromosome"
    );
    fs::create_dir(&a.directory)?;
    let mut sequence = reference.sequence.clone();
    let mut state = a.seed;
    let mut intervals = Vec::new();
    for i in 0..a.gaps {
        let length = a.lengths[i % a.lengths.len()];
        ensure!(length < sequence.len(), "Mask longer than chromosome");
        let mut placed = false;
        for _ in 0..1_000_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let start = state as usize % (sequence.len() - length + 1);
            let end = start + length;
            if sequence[start.saturating_sub(1)..(end + 1).min(sequence.len())]
                .iter()
                .any(|b| b.eq_ignore_ascii_case(&b'N'))
            {
                continue;
            }
            sequence[start..end].fill(b'N');
            intervals.push((start, end));
            placed = true;
            break;
        }
        ensure!(placed, "Could not place non-overlapping masks");
    }
    intervals.sort_unstable();
    let mut fasta = BufWriter::new(File::create(a.directory.join("masked.fa"))?);
    fasta.write_all(b">")?;
    fasta.write_all(&reference.header)?;
    fasta.write_all(b"\n")?;
    for line in sequence.chunks(80) {
        fasta.write_all(line)?;
        fasta.write_all(b"\n")?;
    }
    fasta.flush()?;
    let mut masks = BufWriter::new(File::create(a.directory.join("masks.tsv"))?);
    writeln!(masks, "start_0based\tend_exclusive\tlength")?;
    for (start, end) in intervals {
        writeln!(masks, "{start}\t{end}\t{}", end - start)?;
    }
    masks.flush()?;
    fs::write(
        a.directory.join("parameters.txt"),
        format!(
            "reference={}\nseed={}\ngaps={}\nmasked_bases={}\nlengths={:?}\n",
            a.reference.display(),
            a.seed,
            a.gaps,
            total,
            a.lengths
        ),
    )?;
    // Round-trip through helicase verifies that masking preserved all coordinates.
    let parsed = input::scaffolds(&a.directory.join("masked.fa"))?;
    ensure!(
        parsed.len() == 1 && parsed[0].header == reference.header && parsed[0].sequence == sequence,
        "Mask output round-trip failed"
    );
    Ok(())
}
