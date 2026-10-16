use crate::{
    Gap, Scaffold,
    consensus::{self, Consensus, Preset},
    discover,
    index::{Index, Scanner, Span},
    input,
};
use anyhow::{Context, Result, ensure};
use clap::Parser;
use crossbeam_channel::bounded;
use helicase::HelicaseParser;
use helicase::parser::Format;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use tempfile::NamedTempFile;

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Conservative parallel scaffold gap filling using exact flank anchors and POASTA"
)]
pub struct Args {
    /// Scaffold FASTA (plain, gzip, zstd, bzip2, or xz; detected by contents)
    #[arg(short, long)]
    pub scaffolds: PathBuf,
    /// Long reads, FASTA or FASTQ (plain or compressed)
    #[arg(short, long)]
    pub reads: PathBuf,
    /// Corrected FASTA; must not already exist
    #[arg(short, long)]
    pub output: PathBuf,
    /// Per-gap TSV report; must not already exist
    #[arg(long)]
    pub report: PathBuf,
    /// Optional per-phase timing TSV; must not already exist
    #[arg(long)]
    pub timings: Option<PathBuf>,
    /// Read technology; selects initial MSA scoring, not acceptance thresholds
    #[arg(long, value_enum)]
    pub preset: Preset,
    /// One adjacent anchor of this length on each side (3..=63)
    #[arg(short = 'k', long)]
    pub k: usize,
    /// Minimum distinct spanning read IDs required (at least 2)
    #[arg(long)]
    pub min_reads: usize,
    /// Minimum unanimous MSA columns, as a percentage in (0,100]
    #[arg(long)]
    pub min_quality: f64,
    /// Minimizer size (default: min(k,15); at most 29)
    #[arg(short = 'm', long)]
    pub minimizer: Option<usize>,
    /// Parallel read workers (default: available CPUs)
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
    /// Concurrent gap MSAs (default: min(threads,2))
    #[arg(long)]
    pub msa_threads: Option<usize>,
    /// Reject a gap if any unambiguous spanning interval exceeds this many bases
    #[arg(long, default_value_t = 50_000)]
    pub max_span: usize,
    /// Reject a gap if more distinct spanning reads would be retained; never subsample
    #[arg(long, default_value_t = 256)]
    pub max_reads_per_gap: usize,
    /// Maximum retained sequence plus read-ID bytes for one gap
    #[arg(long, default_value_t = 4_194_304)]
    pub max_gap_bytes: usize,
    /// Maximum total retained evidence bytes (sequence plus IDs)
    #[arg(long, default_value_t = 536_870_912)]
    pub max_evidence_bytes: usize,
    /// Conservative bound sum(span lengths)*max(span length) before a nontrivial MSA
    #[arg(long, default_value_t = 20_000_000)]
    pub max_msa_cells: usize,
    /// Target bases per input batch; two batches may be queued
    #[arg(long, default_value_t = 4_194_304)]
    pub batch_bases: usize,
    /// Fail on longer input reads (checked after helicase parses the record)
    #[arg(long, default_value_t = 10_000_000)]
    pub max_read_bases: usize,
}

impl Args {
    pub fn validate(&self) -> Result<()> {
        ensure!((3..=63).contains(&self.k), "k must be in 3..=63");
        let m = self.minimizer.unwrap_or(self.k.min(15));
        ensure!(
            m > 0 && m <= self.k && m <= 29,
            "minimizer must be in 1..=min(k,29)"
        );
        ensure!(self.min_reads >= 2, "min-reads must be at least 2");
        ensure!(
            self.min_quality.is_finite() && self.min_quality > 0.0 && self.min_quality <= 100.0,
            "min-quality must be a finite percentage in (0,100]"
        );
        ensure!(
            self.threads != Some(0) && self.msa_threads != Some(0),
            "Thread counts must be positive"
        );
        ensure!(
            self.max_reads_per_gap >= self.min_reads,
            "max-reads-per-gap must be at least min-reads"
        );
        ensure!(
            [
                self.max_span,
                self.max_gap_bytes,
                self.max_evidence_bytes,
                self.max_msa_cells,
                self.batch_bases,
                self.max_read_bases
            ]
            .iter()
            .all(|&n| n > 0),
            "Resource limits must be positive"
        );
        ensure!(
            self.max_span <= self.max_read_bases && self.max_read_bases <= u32::MAX as usize / 2,
            "max-span must be <= max-read-bases; max-read-bases must be <= 2^31-1"
        );
        for path in [&self.output, &self.report] {
            ensure!(
                std::fs::symlink_metadata(path).is_err(),
                "Output already exists: {}",
                path.display()
            );
            normalized_destination(path)?;
        }
        ensure!(
            normalized_destination(&self.output)? != normalized_destination(&self.report)?,
            "FASTA and report must have different paths"
        );
        if let Some(path) = &self.timings {
            ensure!(
                std::fs::symlink_metadata(path).is_err(),
                "Timing output already exists: {}",
                path.display()
            );
            let destination = normalized_destination(path)?;
            ensure!(
                destination != normalized_destination(&self.output)?
                    && destination != normalized_destination(&self.report)?,
                "Timing output must have a separate path"
            );
        }
        Ok(())
    }
}

fn normalized_destination(path: &Path) -> Result<PathBuf> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = path.file_name().context("Output requires a filename")?;
    Ok(parent
        .canonicalize()
        .with_context(|| format!("Output directory does not exist: {}", parent.display()))?
        .join(name))
}

struct Read {
    id: Vec<u8>,
    sequence: Vec<u8>,
}

#[derive(Default)]
struct ReadTimings {
    records: usize,
    parser: f64,
    producer_wait: f64,
    input_wait: f64,
    extraction: f64,
    merge: f64,
}

#[derive(Default)]
struct Evidence {
    candidates: usize,
    spanning_records: usize,
    ambiguous_pairs: usize,
    duplicate_ids: usize,
    ids: Vec<Vec<u8>>,
    sequences: Vec<Vec<u8>>,
    bytes: usize,
    rejection: Option<&'static str>,
}

impl Evidence {
    fn reject(&mut self, reason: &'static str, total: &mut usize) {
        self.rejection = Some(reason);
        *total -= self.bytes;
        self.bytes = 0;
        self.ids.clear();
        self.sequences.clear();
    }

    fn add(&mut self, span: Span, id: &[u8], args: &Args, total: &mut usize) {
        self.candidates += 1;
        if matches!(span, Span::AmbiguousPair) {
            self.ambiguous_pairs += 1;
            return;
        }
        self.spanning_records += 1;
        if self.rejection.is_some() {
            return;
        }
        let sequence = match span {
            Span::Sequence(seq) => seq,
            Span::Rejected(reason) => {
                self.reject(reason, total);
                return;
            }
            Span::AmbiguousPair => unreachable!(),
        };
        if let Some(previous) = self.ids.iter().position(|old| old == id) {
            self.duplicate_ids += 1;
            if self.sequences[previous] != sequence {
                self.reject("conflicting_duplicate_read_id", total);
            }
            return;
        }
        let bytes = sequence.len() + id.len();
        let reason = if self.sequences.len() == args.max_reads_per_gap {
            Some("read_support_limit")
        } else if bytes > args.max_gap_bytes.saturating_sub(self.bytes) {
            Some("gap_evidence_limit")
        } else if bytes > args.max_evidence_bytes.saturating_sub(*total) {
            Some("global_evidence_limit")
        } else {
            None
        };
        if let Some(reason) = reason {
            self.reject(reason, total);
            return;
        }
        self.bytes += bytes;
        *total += bytes;
        self.ids.push(id.to_vec());
        self.sequences.push(sequence);
    }
}

fn produce_reads(
    args: &Args,
    tx: &crossbeam_channel::Sender<Result<Vec<Read>>>,
) -> Result<ReadTimings> {
    let mut clock = Instant::now();
    let mut timings = ReadTimings::default();
    let mut parser = input::open(&args.reads)?;
    let mut batch = Vec::new();
    let mut bases = 0;
    let mut count = 0;
    while parser.next().is_some() {
        let sequence = parser.get_dna_string();
        ensure!(
            sequence.len() <= args.max_read_bases,
            "Read exceeds --max-read-bases"
        );
        let id = parser
            .get_header()
            .split(|b| b.is_ascii_whitespace())
            .next()
            .unwrap_or_default();
        ensure!(!id.is_empty(), "Read has an empty ID");
        if matches!(parser.format(), Format::Fastq) {
            ensure!(
                parser
                    .get_quality()
                    .is_some_and(|q| q.len() == sequence.len()),
                "FASTQ sequence and quality lengths differ"
            );
        }
        bases += sequence.len() + id.len();
        batch.push(Read {
            id: id.to_vec(),
            sequence: sequence.to_vec(),
        });
        count += 1;
        if bases >= args.batch_bases {
            timings.parser += clock.elapsed().as_secs_f64();
            let wait = Instant::now();
            if tx.send(Ok(std::mem::take(&mut batch))).is_err() {
                timings.records = count;
                return Ok(timings);
            }
            timings.producer_wait += wait.elapsed().as_secs_f64();
            clock = Instant::now();
            bases = 0;
        }
    }
    timings.parser += clock.elapsed().as_secs_f64();
    if !batch.is_empty() {
        let wait = Instant::now();
        let _ = tx.send(Ok(batch));
        timings.producer_wait += wait.elapsed().as_secs_f64();
    }
    ensure!(count > 0, "Read input contains no records");
    timings.records = count;
    Ok(timings)
}

fn collect(
    args: &Args,
    index: &Index,
    gaps: &[Gap],
    pool: &rayon::ThreadPool,
) -> Result<(Vec<Evidence>, ReadTimings)> {
    let mut evidence: Vec<_> = gaps
        .iter()
        .map(|gap| Evidence {
            rejection: gap.rejection,
            ..Evidence::default()
        })
        .collect();
    let mut total = 0;
    let timings = std::thread::scope(|scope| -> Result<ReadTimings> {
        let (tx, rx) = bounded(2);
        let producer = scope.spawn(move || {
            let result = produce_reads(args, &tx);
            match result {
                Ok(timings) => Ok(timings),
                Err(error) => {
                    let message = error.to_string();
                    let _ = tx.send(Err(error));
                    Err(anyhow::anyhow!(message))
                }
            }
        });
        let mut main_timings = ReadTimings::default();
        let result = (|| -> Result<()> {
            loop {
                let wait = Instant::now();
                let batch = rx.recv();
                main_timings.input_wait += wait.elapsed().as_secs_f64();
                let Ok(batch) = batch else {
                    break;
                };
                let batch = batch?;
                let clock = Instant::now();
                let spans: Vec<_> = pool.install(|| {
                    batch
                        .par_iter()
                        .map_init(Scanner::default, |scanner, read| {
                            index.extract(&read.sequence, scanner, args.max_span)
                        })
                        .collect()
                });
                main_timings.extraction += clock.elapsed().as_secs_f64();
                let clock = Instant::now();
                // Indexed parallel collect and ordered merge ensure that capacity
                // decisions and read support do not depend on scheduling.
                for (read, found) in batch.into_iter().zip(spans) {
                    for (gap, span) in found {
                        evidence[gap].add(span, &read.id, args, &mut total);
                    }
                }
                main_timings.merge += clock.elapsed().as_secs_f64();
            }
            Ok(())
        })();
        drop(rx); // unblock producer on any processing error
        let mut timings = producer
            .join()
            .map_err(|_| anyhow::anyhow!("Read parser panicked; no outputs were written"))??;
        result?;
        timings.input_wait = main_timings.input_wait;
        timings.extraction = main_timings.extraction;
        timings.merge = main_timings.merge;
        Ok(timings)
    })?;
    Ok((evidence, timings))
}

struct Decision {
    support: usize,
    consensus: Option<Consensus>,
    rejection: Option<&'static str>,
}

fn decide(evidence: &mut Evidence, args: &Args) -> Result<Decision> {
    let support = evidence.sequences.len();
    let mut rejection = evidence.rejection;
    if rejection.is_none() && support < args.min_reads {
        rejection = Some("insufficient_support");
    }
    if rejection.is_some() {
        return Ok(Decision {
            support,
            consensus: None,
            rejection,
        });
    }
    let identical = evidence
        .sequences
        .iter()
        .all(|s| s == &evidence.sequences[0]);
    let max_len = evidence.sequences.iter().map(Vec::len).max().unwrap();
    let sum_len: usize = evidence.sequences.iter().map(Vec::len).sum();
    if !identical
        && sum_len
            .checked_mul(max_len)
            .is_none_or(|cells| cells > args.max_msa_cells)
    {
        return Ok(Decision {
            support,
            consensus: None,
            rejection: Some("msa_size_limit"),
        });
    }
    let consensus = consensus::build(&mut evidence.sequences, args.preset)?;
    rejection = consensus.rejection;
    if rejection.is_none() && consensus.quality < args.min_quality {
        rejection = Some("low_msa_quality");
    }
    Ok(Decision {
        support,
        consensus: Some(consensus),
        rejection,
    })
}

fn tsv_id(bytes: &[u8]) -> String {
    let mut out = String::new();
    for &b in bytes {
        if (b'!'..=b'~').contains(&b) && b != b'%' {
            out.push(b as char);
        } else {
            use std::fmt::Write;
            write!(out, "%{b:02X}").unwrap();
        }
    }
    out
}

fn write_outputs(
    args: &Args,
    scaffolds: &[Scaffold],
    gaps: &[Gap],
    evidence: &[Evidence],
    decisions: &[Decision],
) -> Result<()> {
    let output = normalized_destination(&args.output)?;
    let report = normalized_destination(&args.report)?;
    let mut fasta_tmp = NamedTempFile::new_in(output.parent().unwrap())?;
    let mut report_tmp = NamedTempFile::new_in(report.parent().unwrap())?;
    {
        let mut writer = BufWriter::new(fasta_tmp.as_file_mut());
        let mut cursor = 0;
        for (scaffold_id, scaffold) in scaffolds.iter().enumerate() {
            writer.write_all(b">")?;
            writer.write_all(&scaffold.header)?;
            writer.write_all(b"\n")?;
            let mut pos = 0;
            while cursor < gaps.len() && gaps[cursor].scaffold == scaffold_id {
                let gap = &gaps[cursor];
                writer.write_all(&scaffold.sequence[pos..gap.start])?;
                if decisions[cursor].rejection.is_none() {
                    writer.write_all(&decisions[cursor].consensus.as_ref().unwrap().sequence)?;
                } else {
                    writer.write_all(&scaffold.sequence[gap.start..gap.end])?;
                }
                pos = gap.end;
                cursor += 1;
            }
            writer.write_all(&scaffold.sequence[pos..])?;
            writer.write_all(b"\n")?;
        }
        writer.flush()?;
    }
    {
        let mut writer = BufWriter::new(report_tmp.as_file_mut());
        writeln!(
            writer,
            "# gapfil={} preset={:?} k={} minimizer={} min_reads={} min_quality_percent={} max_span={} max_reads_per_gap={} max_gap_bytes={} max_evidence_bytes={} max_msa_cells={}",
            env!("CARGO_PKG_VERSION"),
            args.preset,
            args.k,
            args.minimizer.unwrap_or(args.k.min(15)),
            args.min_reads,
            args.min_quality,
            args.max_span,
            args.max_reads_per_gap,
            args.max_gap_bytes,
            args.max_evidence_bytes,
            args.max_msa_cells
        )?;
        writeln!(
            writer,
            "gap_id\tscaffold_index\tscaffold_id\tstart_0based\tend_exclusive\toriginal_length\tcandidate_records\tspanning_records\tambiguous_pairs\tduplicate_read_ids\tretained_reads\tconsensus_length\tunanimous_columns\tmsa_columns\tmsa_quality_percent\tstatus\treason"
        )?;
        for (id, ((gap, e), d)) in gaps.iter().zip(evidence).zip(decisions).enumerate() {
            let header = scaffolds[gap.scaffold]
                .header
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            let values = if let Some(c) = &d.consensus {
                format!(
                    "{}\t{}\t{}\t{:.8}",
                    c.sequence.len(),
                    c.unanimous,
                    c.columns,
                    c.quality
                )
            } else {
                "NA\tNA\tNA\tNA".into()
            };
            let (status, reason) = d
                .rejection
                .map_or(("filled", "accepted"), |r| ("unchanged", r));
            writeln!(
                writer,
                "{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{values}\t{status}\t{reason}",
                gap.scaffold,
                tsv_id(header),
                gap.start,
                gap.end,
                gap.end - gap.start,
                e.candidates,
                e.spanning_records,
                e.ambiguous_pairs,
                e.duplicate_ids,
                d.support
            )?;
        }
        writer.flush()?;
    }
    // Each file is published complete, using no-clobber creation. Two separate
    // filesystem paths cannot be published as one atomic transaction.
    report_tmp
        .persist_noclobber(&report)
        .with_context(|| format!("Cannot publish {}", report.display()))?;
    fasta_tmp.persist_noclobber(&output).with_context(|| {
        format!(
            "Cannot publish {}; report was already published at {}",
            output.display(),
            report.display()
        )
    })?;
    Ok(())
}

pub fn run(args: &Args) -> Result<()> {
    let start = Instant::now();
    args.validate()?;
    let mut phases = Vec::new();
    let mut phase_clock = start;
    let threads = args
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    let msa_threads = args.msa_threads.unwrap_or(threads.min(2)).min(threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    mark_phase(&mut phases, &mut phase_clock, "setup");
    let scaffolds = input::scaffolds(&args.scaffolds)?;
    mark_phase(&mut phases, &mut phase_clock, "scaffold_read");
    let mut gaps = discover(&scaffolds, args.k);
    mark_phase(&mut phases, &mut phase_clock, "gap_discovery");
    let m = args.minimizer.unwrap_or(args.k.min(15));
    let index = Index::new(&gaps, args.k, m);
    mark_phase(&mut phases, &mut phase_clock, "index_build");
    index.reject_repeated_anchors(&scaffolds, &mut gaps);
    mark_phase(&mut phases, &mut phase_clock, "anchor_uniqueness");
    drop(index);
    let index = Index::new(&gaps, args.k, m);
    mark_phase(&mut phases, &mut phase_clock, "index_filter");
    eprintln!(
        "Indexed {} gaps ({} eligible) in {:.3}s; {} read workers, {} MSA workers",
        gaps.len(),
        gaps.iter().filter(|g| g.rejection.is_none()).count(),
        start.elapsed().as_secs_f64(),
        threads,
        msa_threads
    );
    let scan_start = Instant::now();
    let (mut evidence, reads) = if gaps.iter().any(|g| g.rejection.is_none()) {
        collect(args, &index, &gaps, &pool)?
    } else {
        (
            gaps.iter()
                .map(|g| Evidence {
                    rejection: g.rejection,
                    ..Evidence::default()
                })
                .collect(),
            ReadTimings::default(),
        )
    };
    drop(index);
    drop(pool);
    eprintln!(
        "Scanned {} reads in {:.3}s",
        reads.records,
        scan_start.elapsed().as_secs_f64()
    );
    mark_phase(&mut phases, &mut phase_clock, "read_scan");
    let msa_start = Instant::now();
    let msa_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(msa_threads)
        .build()?;
    let decisions: Vec<Decision> = msa_pool.install(|| {
        evidence
            .par_iter_mut()
            .enumerate()
            .map(|(id, e)| {
                let decision = decide(e, args).with_context(|| {
                    format!(
                        "Consensus failed for gap {id}, scaffold {}, interval {}..{}",
                        gaps[id].scaffold, gaps[id].start, gaps[id].end
                    )
                });
                e.ids.clear();
                e.sequences.clear();
                decision
            })
            .collect::<Result<Vec<_>>>()
    })?;
    let filled = decisions.iter().filter(|d| d.rejection.is_none()).count();
    eprintln!(
        "Consensus: filled {filled}/{} gaps in {:.3}s",
        gaps.len(),
        msa_start.elapsed().as_secs_f64()
    );
    mark_phase(&mut phases, &mut phase_clock, "consensus");
    write_outputs(args, &scaffolds, &gaps, &evidence, &decisions)?;
    mark_phase(&mut phases, &mut phase_clock, "output");
    let total = start.elapsed().as_secs_f64();
    eprintln!(
        "Phases(s): {}",
        phases
            .iter()
            .map(|(name, seconds)| format!("{name}={seconds:.6}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
    if let Some(path) = &args.timings {
        let destination = normalized_destination(path)?;
        let mut file = NamedTempFile::new_in(destination.parent().unwrap())?;
        writeln!(file, "phase\tseconds\tscope\tread_workers\tmsa_workers")?;
        for (name, seconds) in &phases {
            writeln!(file, "{name}\t{seconds:.9}\twall\t{threads}\t{msa_threads}")?;
        }
        for (name, seconds) in [
            ("read_parse_active", reads.parser),
            ("read_producer_wait", reads.producer_wait),
            ("read_input_wait", reads.input_wait),
            ("read_extraction", reads.extraction),
            ("read_merge", reads.merge),
        ] {
            writeln!(
                file,
                "{name}\t{seconds:.9}\toverlapping\t{threads}\t{msa_threads}"
            )?;
        }
        writeln!(file, "total\t{total:.9}\ttotal\t{threads}\t{msa_threads}")?;
        file.persist_noclobber(destination)?;
    }
    eprintln!("Finished in {:.3}s", start.elapsed().as_secs_f64());
    Ok(())
}

fn mark_phase(phases: &mut Vec<(&'static str, f64)>, clock: &mut Instant, name: &'static str) {
    phases.push((name, clock.elapsed().as_secs_f64()));
    *clock = Instant::now();
}
