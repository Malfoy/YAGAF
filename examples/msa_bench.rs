use anyhow::{Result, ensure};
use clap::Parser;
use rayon::prelude::*;
use std::{fs::File, io::Write, path::PathBuf, time::Instant};
use yagaf::{
    consensus::{self, MsaBackend, Preset},
    input,
};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    reference: PathBuf,
    #[arg(long, value_enum)]
    backend: MsaBackend,
    #[arg(long, value_enum, default_value = "hifi")]
    preset: Preset,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = 500)]
    length: usize,
    #[arg(long, default_value_t = 32)]
    families: usize,
    #[arg(long, default_value_t = 50)]
    depth: usize,
    #[arg(long, default_value_t = 32)]
    threads: usize,
    #[arg(long, default_value_t = 3)]
    repeats: usize,
    #[arg(long, default_value_t = 0.001)]
    error_rate: f64,
}
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn uniform(&mut self) -> f64 {
        ((self.next() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
}
fn main() -> Result<()> {
    let options = Options::parse();
    ensure!(
        options.length > 0
            && options.depth >= 3
            && options.families > 0
            && options.threads > 0
            && options.repeats > 0,
        "Invalid dimensions"
    );
    ensure!(
        options.error_rate.is_finite() && (0.0..=0.1).contains(&options.error_rate),
        "Invalid error rate"
    );
    let reference = input::scaffolds(&options.reference)?;
    let sequence = &reference[0].sequence;
    let mut rng = Rng(42);
    let mut truths = Vec::new();
    let mut families = Vec::new();
    for _ in 0..options.families {
        let start = rng.next() as usize % sequence.len();
        let truth: Vec<_> = (0..options.length)
            .map(|i| sequence[(start + i) % sequence.len()].to_ascii_uppercase())
            .collect();
        let mut reads = Vec::new();
        for _ in 0..options.depth {
            let mut read = Vec::new();
            for &base in &truth {
                if rng.uniform() < options.error_rate {
                    match rng.next() % 3 {
                        0 => {
                            let i = b"ACGT".iter().position(|&b| b == base).unwrap();
                            read.push(b"ACGT"[(i + 1 + (rng.next() % 3) as usize) % 4]);
                        }
                        1 => {
                            read.push(b"ACGT"[(rng.next() % 4) as usize]);
                            read.push(base);
                        }
                        _ => (),
                    }
                } else {
                    read.push(base);
                }
            }
            ensure!(!read.is_empty(), "Empty generated read");
            reads.push(read);
        }
        if reads.iter().all(|r| r == &truth) {
            reads[0][0] = if truth[0] == b'A' { b'C' } else { b'A' };
        }
        truths.push(truth);
        families.push(reads);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    let mut output = File::create_new(&options.output)?;
    writeln!(
        output,
        "backend\tpreset\terror_rate\tlength\tfamilies\tdepth\tthreads\trepetition\tseconds\texact_consensuses\trejected_consensuses\tunanimous_columns\tmsa_columns"
    )?;
    for repetition in 0..options.repeats {
        let mut reads = families.clone();
        let start = Instant::now();
        let results = pool.install(|| {
            reads
                .par_iter_mut()
                .map(|r| consensus::build_with_backend(r, options.preset, options.backend))
                .collect::<Result<Vec<_>>>()
        })?;
        let elapsed = start.elapsed().as_secs_f64();
        let exact = results
            .iter()
            .zip(&truths)
            .filter(|(r, t)| r.sequence == **t)
            .count();
        let rejected = results.iter().filter(|r| r.rejection.is_some()).count();
        let unanimous: usize = results.iter().map(|r| r.unanimous).sum();
        let columns: usize = results.iter().map(|r| r.columns).sum();
        let row = format!(
            "{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{elapsed:.9}\t{exact}\t{rejected}\t{unanimous}\t{columns}",
            options.backend.name(),
            options.preset,
            options.error_rate,
            options.length,
            options.families,
            options.depth,
            options.threads,
            repetition + 1
        );
        writeln!(output, "{row}")?;
        println!("{row}");
    }
    Ok(())
}
