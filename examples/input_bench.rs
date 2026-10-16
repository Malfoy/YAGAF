//! Compare streamed parsing with helicase record framing and parallel parsing.
use anyhow::{Result, ensure};
use clap::Parser;
use helicase::{FastxParser, HelicaseParser, ParserOptions, input::*, parser::Event};
use rayon::prelude::*;
use std::{hint::black_box, path::PathBuf, time::Instant};
use yagaf::input;

#[derive(Parser)]
struct Options {
    #[arg(long)]
    reads: PathBuf,
    #[arg(long, default_value_t = 32)]
    threads: usize,
    #[arg(long, default_value_t = 5)]
    repeats: usize,
}
const FRAME_CONFIG: u64 = ParserOptions::default()
    .ignore_dna()
    .ignore_headers()
    .config();
fn main() -> Result<()> {
    let options = Options::parse();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    for _ in 0..options.repeats {
        let start = Instant::now();
        let mut parser = input::open(&options.reads)?;
        let mut streamed = (0usize, 0usize);
        while parser.next().is_some() {
            let sequence = parser.get_dna_string_owned();
            streamed.0 += 1;
            streamed.1 += black_box(sequence).len();
        }
        let stream_time = start.elapsed().as_secs_f64();
        let start = Instant::now();
        let mapped = MmapInput::new(&options.reads)?;
        let data = mapped.data();
        let mut parser = FastxParser::<FRAME_CONFIG>::from_slice(data)?;
        let mut ranges = Vec::new();
        let mut offset = 0;
        for event in &mut parser {
            if let Event::Record(range) = event {
                ensure!(range.end <= data.len(), "Record exceeds input");
                if range.end - offset >= 1_048_576 {
                    ranges.push(offset..range.end);
                    offset = range.end;
                }
            }
        }
        if offset < data.len() {
            ranges.push(offset..data.len());
        }
        let frame_time = start.elapsed().as_secs_f64();
        let counts = pool.install(|| {
            ranges
                .par_iter()
                .map(|range| -> Result<(usize, usize)> {
                    let mut parser =
                        FastxParser::<{ input::CONFIG }>::from_slice(&data[range.clone()])?;
                    let mut counts = (0, 0);
                    while parser.next().is_some() {
                        counts.0 += 1;
                        counts.1 += black_box(parser.get_dna_string()).len();
                    }
                    Ok(counts)
                })
                .collect::<Result<Vec<_>>>()
        })?;
        let parallel = counts
            .into_iter()
            .fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
        ensure!(parallel == streamed, "Parser counts differ");
        println!(
            "stream={stream_time:.6} framing={frame_time:.6} parallel_total={:.6} records={} bases={}",
            start.elapsed().as_secs_f64(),
            parallel.0,
            parallel.1
        );
    }
    Ok(())
}
