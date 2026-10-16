use crate::{
    Gap, Scaffold,
    consensus::{self, Consensus, MsaBackend, Preset},
    discover,
    index::{Index, Scanner, Sequence, Span},
    input,
};
use anyhow::{Context, Result, ensure};
use clap::Parser;
use crossbeam_channel::bounded;
use helicase::input::{FromSlice, InputData};
use helicase::parser::{Event, Format};
use helicase::{FastxParser, HelicaseParser, ParserOptions};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHasher};
use std::borrow::Cow;
use std::cell::RefCell;
use std::hash::Hasher;
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use tempfile::NamedTempFile;

#[derive(Debug, Parser)]
#[command(
    name = "yagaf",
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
    /// Per-gap TSV report (default: OUTPUT.gaps.tsv); must not already exist
    #[arg(long)]
    pub report: Option<PathBuf>,
    /// Optional per-phase timing TSV; must not already exist
    #[arg(long, hide_short_help = true)]
    pub timings: Option<PathBuf>,
    /// Read technology; selects initial MSA scoring, not acceptance thresholds
    #[arg(long, value_enum, default_value = "hifi", hide_short_help = true)]
    pub preset: Preset,
    /// MSA implementation; alternative engines are available for comparison
    #[arg(long, value_enum, default_value = "poasta", hide_short_help = true)]
    pub msa_backend: MsaBackend,
    /// One adjacent anchor of this length on each side (3..=63)
    #[arg(short = 'k', long, default_value = "31", hide_short_help = true)]
    pub k: usize,
    /// Minimum distinct spanning read IDs required (at least 2)
    #[arg(long, default_value_t = 3, hide_short_help = true)]
    pub min_reads: usize,
    /// Minimum unanimous MSA columns, as a percentage in (0,100]
    #[arg(long, default_value_t = 50.0, hide_short_help = true)]
    pub min_quality: f64,
    /// Minimizer size (default: min(k,15); at most 29)
    #[arg(short = 'm', long, hide_short_help = true)]
    pub minimizer: Option<usize>,
    /// Parallel read workers (default: available CPUs)
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
    /// Concurrent gap MSAs (default: threads)
    #[arg(long, hide_short_help = true)]
    pub msa_threads: Option<usize>,
    /// Reject a gap if any unambiguous spanning interval exceeds this many bases
    #[arg(long, default_value_t = 50_000, hide_short_help = true)]
    pub max_span: usize,
    /// Reject a gap if more distinct spanning reads would be retained; never subsample
    #[arg(long, default_value_t = 256, hide_short_help = true)]
    pub max_reads_per_gap: usize,
    /// Maximum retained sequence plus read-ID bytes for one gap
    #[arg(long, default_value_t = 4_194_304, hide_short_help = true)]
    pub max_gap_bytes: usize,
    /// Maximum total retained evidence bytes (sequence plus IDs)
    #[arg(long, default_value_t = 536_870_912, hide_short_help = true)]
    pub max_evidence_bytes: usize,
    /// Conservative bound sum(span lengths)*max(span length) before a nontrivial MSA
    #[arg(long, default_value_t = 20_000_000, hide_short_help = true)]
    pub max_msa_cells: usize,
    /// Target bases per input batch; two batches may be queued
    #[arg(long, default_value_t = 16_777_216, hide_short_help = true)]
    pub batch_bases: usize,
    /// Fail on longer input reads (checked after helicase parses the record)
    #[arg(long, default_value_t = 10_000_000, hide_short_help = true)]
    pub max_read_bases: usize,
    /// Stream uncompressed reads instead of memory mapping and parallel parsing
    #[arg(long, hide_short_help = true)]
    pub stream_reads: bool,
}

impl Args {
    pub fn report_path(&self) -> Cow<'_, Path> {
        self.report.as_deref().map_or_else(
            || {
                let mut path = self.output.as_os_str().to_os_string();
                path.push(".gaps.tsv");
                Cow::Owned(PathBuf::from(path))
            },
            Cow::Borrowed,
        )
    }

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
        let report = self.report_path();
        for path in [self.output.as_path(), report.as_ref()] {
            ensure!(
                std::fs::symlink_metadata(path).is_err(),
                "Output already exists: {}",
                path.display()
            );
            normalized_destination(path)?;
        }
        ensure!(
            normalized_destination(&self.output)? != normalized_destination(&report)?,
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
                    && destination != normalized_destination(&report)?,
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
    id: Arc<[u8]>,
    id_hash: u64,
    sequence: Vec<u8>,
}

enum ReadBatch {
    Records(Vec<Read>),
    Encoded(Vec<Range<usize>>),
}

fn read_id(parser: &impl HelicaseParser, args: &Args) -> Result<(Arc<[u8]>, u64)> {
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
    let mut hasher = FxHasher::default();
    hasher.write(id);
    Ok((Arc::from(id), hasher.finish()))
}

thread_local! {
    static READ_SCANNER: RefCell<Scanner> = RefCell::new(Scanner::default());
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
    ids: Vec<Arc<[u8]>>,
    id_index: FxHashMap<u64, usize>,
    sequences: Vec<Sequence>,
    bytes: usize,
    rejection: Option<&'static str>,
}

impl Evidence {
    fn reject(&mut self, reason: &'static str, total: &mut usize) {
        self.rejection = Some(reason);
        *total -= self.bytes;
        self.bytes = 0;
        self.ids.clear();
        self.id_index.clear();
        self.sequences.clear();
    }

    fn add(&mut self, span: Span, read: &Read, args: &Args, total: &mut usize) {
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
        let id = &read.id;
        let previous = self.id_index.get(&read.id_hash).and_then(|&first| {
            if self.ids[first].as_ref() == id.as_ref() {
                Some(first)
            } else {
                self.ids.iter().position(|old| old.as_ref() == id.as_ref())
            }
        });
        if let Some(previous) = previous {
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
        self.id_index.entry(read.id_hash).or_insert(self.ids.len());
        self.ids.push(Arc::clone(id));
        self.sequences.push(sequence);
    }
}

fn produce_reads(
    args: &Args,
    tx: &crossbeam_channel::Sender<Result<ReadBatch>>,
) -> Result<ReadTimings> {
    let mut clock = Instant::now();
    let mut timings = ReadTimings::default();
    let mut parser = input::open(&args.reads)?;
    let mut batch = Vec::new();
    let mut bases = 0;
    let mut count = 0;
    while parser.next().is_some() {
        let (id, id_hash) = read_id(&parser, args)?;
        bases += parser.get_dna_string().len() + id.len();
        batch.push(Read {
            id,
            id_hash,
            sequence: parser.get_dna_string_owned(),
        });
        count += 1;
        if bases >= args.batch_bases {
            timings.parser += clock.elapsed().as_secs_f64();
            let wait = Instant::now();
            if tx
                .send(Ok(ReadBatch::Records(std::mem::take(&mut batch))))
                .is_err()
            {
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
        let _ = tx.send(Ok(ReadBatch::Records(batch)));
        timings.producer_wait += wait.elapsed().as_secs_f64();
    }
    ensure!(count > 0, "Read input contains no records");
    timings.records = count;
    Ok(timings)
}

fn frame_reads(
    args: &Args,
    data: &[u8],
    tx: &crossbeam_channel::Sender<Result<ReadBatch>>,
    threads: usize,
) -> Result<ReadTimings> {
    const CONFIG: u64 = ParserOptions::default()
        .ignore_dna()
        .ignore_headers()
        .config();
    let mut parser = FastxParser::<CONFIG>::from_slice(data)?;
    let chunk_size = args
        .batch_bases
        .div_ceil(threads.saturating_mul(2))
        .clamp(65_536, 1_048_576);
    let mut timings = ReadTimings::default();
    let mut clock = Instant::now();
    let mut batch = Vec::new();
    let mut offset = 0;
    let mut batch_start = 0;
    for event in &mut parser {
        if let Event::Record(range) = event {
            ensure!(range.end <= data.len(), "Record exceeds input bounds");
            timings.records += 1;
            if range.end - offset >= chunk_size {
                batch.push(offset..range.end);
                offset = range.end;
                if offset - batch_start >= args.batch_bases {
                    timings.parser += clock.elapsed().as_secs_f64();
                    let wait = Instant::now();
                    if tx
                        .send(Ok(ReadBatch::Encoded(std::mem::take(&mut batch))))
                        .is_err()
                    {
                        return Ok(timings);
                    }
                    timings.producer_wait += wait.elapsed().as_secs_f64();
                    clock = Instant::now();
                    batch_start = offset;
                }
            }
        }
    }
    if offset < data.len() {
        batch.push(offset..data.len());
    }
    timings.parser += clock.elapsed().as_secs_f64();
    if !batch.is_empty() {
        let wait = Instant::now();
        let _ = tx.send(Ok(ReadBatch::Encoded(batch)));
        timings.producer_wait += wait.elapsed().as_secs_f64();
    }
    ensure!(timings.records > 0, "Read input contains no records");
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
    let mut buckets: Vec<Vec<(usize, Span)>> = (0..gaps.len()).map(|_| Vec::new()).collect();
    let mapped = if args.stream_reads {
        None
    } else {
        input::mmap_uncompressed(&args.reads)?
    };
    let data = mapped.as_ref().map(InputData::data);
    let timings = std::thread::scope(|scope| -> Result<ReadTimings> {
        let (tx, rx) = bounded(2);
        let producer = scope.spawn(move || {
            let result = if let Some(data) = data {
                frame_reads(args, data, &tx, pool.current_num_threads())
            } else {
                produce_reads(args, &tx)
            };
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
                let incoming = batch?;
                let clock = Instant::now();
                let (batch, spans): (Vec<Read>, Vec<Vec<(usize, Span)>>) = match incoming {
                    ReadBatch::Records(batch) => {
                        let spans = pool.install(|| {
                            batch
                                .par_iter()
                                .map(|read| {
                                    READ_SCANNER.with_borrow_mut(|scanner| {
                                        index.extract(&read.sequence, scanner, args.max_span)
                                    })
                                })
                                .collect()
                        });
                        (batch, spans)
                    }
                    ReadBatch::Encoded(ranges) => {
                        let data = data.expect("Encoded batches require a mapped input");
                        let chunks = pool.install(|| {
                            ranges
                                .par_iter()
                                .map(|range| -> Result<Vec<_>> {
                                    let mut parser = FastxParser::<{ input::CONFIG }>::from_slice(
                                        &data[range.clone()],
                                    )?;
                                    let mut records = Vec::new();
                                    while parser.next().is_some() {
                                        let (id, id_hash) = read_id(&parser, args)?;
                                        let found = READ_SCANNER.with_borrow_mut(|scanner| {
                                            index.extract(
                                                parser.get_dna_string(),
                                                scanner,
                                                args.max_span,
                                            )
                                        });
                                        records.push((
                                            Read {
                                                id,
                                                id_hash,
                                                sequence: Vec::new(),
                                            },
                                            found,
                                        ));
                                    }
                                    Ok(records)
                                })
                                .collect::<Result<Vec<_>>>()
                        })?;
                        chunks.into_iter().flatten().unzip()
                    }
                };
                main_timings.extraction += clock.elapsed().as_secs_f64();
                let clock = Instant::now();
                let mut potential_bytes = 0usize;
                let mut hits = 0;
                for (read, found) in batch.iter().zip(&spans) {
                    hits += found.len();
                    for (_, span) in found {
                        if let Span::Sequence(sequence) = span {
                            potential_bytes = potential_bytes
                                .saturating_add(sequence.len())
                                .saturating_add(read.id.len());
                        }
                    }
                }
                // If even the upper bound fits, global capacity cannot affect
                // any decision. Merge disjoint gaps concurrently in read order.
                // Near the global limit use the original ordered serial path.
                if hits >= 4096
                    && pool.current_num_threads() > 1
                    && potential_bytes <= args.max_evidence_bytes.saturating_sub(total)
                {
                    for (read_index, found) in spans.into_iter().enumerate() {
                        for (gap, span) in found {
                            buckets[gap].push((read_index, span));
                        }
                    }
                    pool.install(|| {
                        evidence
                            .par_iter_mut()
                            .zip(buckets.par_iter_mut())
                            .for_each(|(e, bucket)| {
                                let mut local_total = e.bytes;
                                for (read_index, span) in bucket.drain(..) {
                                    e.add(span, &batch[read_index], args, &mut local_total);
                                }
                            })
                    });
                    total = evidence.iter().map(|e| e.bytes).sum();
                } else {
                    for (read, found) in batch.iter().zip(spans) {
                        for (gap, span) in found {
                            evidence[gap].add(span, read, args, &mut total);
                        }
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
    let max_len = evidence.sequences.iter().map(|s| s.len()).max().unwrap();
    let sum_len: usize = evidence.sequences.iter().map(|s| s.len()).sum();
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
    let consensus =
        consensus::build_with_backend(&mut evidence.sequences, args.preset, args.msa_backend)?;
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
    let report = normalized_destination(&args.report_path())?;
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
            "# yagaf={} preset={:?} k={} minimizer={} min_reads={} min_quality_percent={} max_span={} max_reads_per_gap={} max_gap_bytes={} max_evidence_bytes={} max_msa_cells={} msa_backend={}",
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
            args.max_msa_cells,
            args.msa_backend.name()
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
    let msa_threads = args.msa_threads.unwrap_or(threads).min(threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    mark_phase(&mut phases, &mut phase_clock, "setup");
    let scaffolds = input::scaffolds(&args.scaffolds)?;
    mark_phase(&mut phases, &mut phase_clock, "scaffold_read");
    let mut gaps = discover(&scaffolds, args.k);
    mark_phase(&mut phases, &mut phase_clock, "gap_discovery");
    let m = args.minimizer.unwrap_or(args.k.min(15));
    let mut index = pool.install(|| Index::new(&gaps, args.k, m));
    mark_phase(&mut phases, &mut phase_clock, "index_build");
    pool.install(|| index.reject_repeated_anchors(&scaffolds, &mut gaps));
    mark_phase(&mut phases, &mut phase_clock, "anchor_uniqueness");
    index.retain_eligible(&gaps);
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
    eprintln!(
        "Scanned {} reads in {:.3}s",
        reads.records,
        scan_start.elapsed().as_secs_f64()
    );
    mark_phase(&mut phases, &mut phase_clock, "read_scan");
    let msa_start = Instant::now();
    let smaller_pool = if msa_threads != threads {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(msa_threads)
                .build()?,
        )
    } else {
        None
    };
    let msa_pool = smaller_pool.as_ref().unwrap_or(&pool);
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
                e.id_index.clear();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colliding_read_id_hashes_keep_distinct_support_and_detect_conflicts() {
        let args = Args::parse_from([
            "yagaf",
            "--scaffolds",
            "scaffolds.fa",
            "--reads",
            "reads.fa",
            "--output",
            "output.fa",
            "--report",
            "gaps.tsv",
            "--preset",
            "hifi",
            "--k",
            "31",
            "--min-reads",
            "2",
            "--min-quality",
            "100",
        ]);
        let read = |id: &[u8]| Read {
            id: Arc::from(id),
            id_hash: 7,
            sequence: Vec::new(),
        };
        let first = read(b"first");
        let second = read(b"second");
        let mut evidence = Evidence::default();
        let mut total = 0;
        evidence.add(
            Span::Sequence(Sequence::from_slice(b"A")),
            &first,
            &args,
            &mut total,
        );
        evidence.add(
            Span::Sequence(Sequence::from_slice(b"C")),
            &second,
            &args,
            &mut total,
        );
        evidence.add(
            Span::Sequence(Sequence::from_slice(b"C")),
            &second,
            &args,
            &mut total,
        );
        assert_eq!(evidence.sequences.len(), 2);
        assert_eq!(evidence.duplicate_ids, 1);
        assert!(evidence.rejection.is_none());
        evidence.add(
            Span::Sequence(Sequence::from_slice(b"T")),
            &second,
            &args,
            &mut total,
        );
        assert_eq!(evidence.rejection, Some("conflicting_duplicate_read_id"));
        assert_eq!(total, 0);
        assert!(evidence.sequences.is_empty());
    }
}
