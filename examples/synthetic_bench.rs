//! Deterministic synthetic data only. Uses helicase for every FASTA read.
use anyhow::{Result, ensure};
use clap::Parser;
use gapfil::{input, revcomp};
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::PathBuf,
    process::Command,
    time::Instant,
};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    directory: PathBuf,
    #[arg(long, default_value_t = 10_000)]
    read_count: usize,
    #[arg(long, default_value_t = 10_000)]
    read_length: usize,
    #[arg(long, default_value_t = 128)]
    gaps: usize,
}

fn dna(seed: &mut u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|_| {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            b"ACGT"[(*seed >> 32) as usize & 3]
        })
        .collect()
}
fn record(writer: &mut impl Write, header: &str, sequence: &[u8]) -> Result<()> {
    writeln!(writer, ">{header}")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn main() -> Result<()> {
    let options = Options::parse();
    ensure!(
        options.read_count >= options.gaps * 5 && options.read_length >= 512 && options.gaps > 0,
        "Need at least 5 reads per gap, positive gap count, and read length >=512"
    );
    fs::create_dir(&options.directory)?;
    let mut seed = 0x123456789abcdef;
    let mut scaffolds = BufWriter::new(File::create(options.directory.join("scaffolds.fa"))?);
    let mut truth = BufWriter::new(File::create(options.directory.join("truth.fa"))?);
    let mut reads = BufWriter::new(File::create(options.directory.join("reads.fa"))?);
    let mut patterns = Vec::new();
    for gap in 0..options.gaps {
        let left = dna(&mut seed, 100);
        let right = dna(&mut seed, 100);
        let middle = dna(&mut seed, 200);
        record(
            &mut scaffolds,
            &format!("s{gap}"),
            &[left.as_slice(), &[b'N'; 100], right.as_slice()].concat(),
        )?;
        record(
            &mut truth,
            &format!("s{gap}"),
            &[left.as_slice(), middle.as_slice(), right.as_slice()].concat(),
        )?;
        patterns.push((left, right, middle));
    }
    for id in 0..options.read_count {
        let mut sequence = dna(&mut seed, options.read_length);
        if id < options.gaps * 5 {
            let (left, right, middle) = &patterns[id / 5];
            let mut span = [left.as_slice(), middle.as_slice(), right.as_slice()].concat();
            if id % 5 == 4 {
                let pos = left.len() + middle.len() / 2;
                span[pos] = if span[pos] == b'A' { b'C' } else { b'A' };
            }
            let offset = (sequence.len() - span.len()) / 2;
            sequence[offset..offset + span.len()].copy_from_slice(&span);
            if id % 2 == 1 {
                sequence = revcomp(&sequence);
            }
        }
        record(&mut reads, &format!("r{id}"), &sequence)?;
    }
    scaffolds.flush()?;
    truth.flush()?;
    reads.flush()?;
    let expected = input::scaffolds(&options.directory.join("truth.fa"))?;
    let executable = std::env::current_exe()?
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("gapfil");
    ensure!(
        executable.exists(),
        "Build gapfil with cargo build --release first"
    );
    let mut summary = File::create(options.directory.join("summary.tsv"))?;
    writeln!(
        summary,
        "threads\tread_bases\tgaps\twall_seconds\tmegabases_per_second\tpeak_rss_kib\texact_fills"
    )?;
    for threads in [1, 4] {
        let out = options.directory.join(format!("out{threads}.fa"));
        let rss = options.directory.join(format!("rss{threads}.txt"));
        let start = Instant::now();
        let status = Command::new("/usr/bin/time")
            .args(["-f", "%M", "-o"])
            .arg(&rss)
            .arg(&executable)
            .arg("--scaffolds")
            .arg(options.directory.join("scaffolds.fa"))
            .arg("--reads")
            .arg(options.directory.join("reads.fa"))
            .arg("--output")
            .arg(&out)
            .arg("--report")
            .arg(options.directory.join(format!("gaps{threads}.tsv")))
            .args([
                "--preset",
                "hifi",
                "--k",
                "31",
                "--min-reads",
                "5",
                "--min-quality",
                "99",
                "--threads",
                &threads.to_string(),
            ])
            .status()?;
        let seconds = start.elapsed().as_secs_f64();
        ensure!(status.success(), "Gapfil benchmark process failed");
        let actual = input::scaffolds(&out)?;
        ensure!(
            actual.len() == expected.len()
                && actual
                    .iter()
                    .zip(&expected)
                    .all(|(a, b)| a.header == b.header && a.sequence == b.sequence),
            "Output differs from synthetic truth"
        );
        let peak = fs::read_to_string(rss)?;
        let bases = options.read_count * options.read_length;
        let row = format!(
            "{threads}\t{bases}\t{}\t{seconds:.6}\t{:.3}\t{}\t{}",
            options.gaps,
            bases as f64 / seconds / 1e6,
            peak.trim(),
            options.gaps
        );
        writeln!(summary, "{row}")?;
        println!("{row}");
    }
    Ok(())
}
