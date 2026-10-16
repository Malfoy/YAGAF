//! Reproducible masking/recovery experiment on a supplied reference.
//! Approximate HiFi-like Q30 reads, not an empirical PacBio CCS model.
//! Every FASTA read goes through helicase.
use anyhow::{Context, Result, ensure};
use clap::Parser;
use gapfil::{discover, input, is_base, revcomp};
use rustc_hash::FxHashSet;
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{BufWriter, Write},
    path::PathBuf,
    process::Command,
    time::Instant,
};

#[derive(Parser, Debug)]
struct Options {
    #[arg(long)]
    reference: PathBuf,
    #[arg(long)]
    directory: PathBuf,
    #[arg(long, default_value_t = 50)]
    coverage: usize,
    #[arg(long, default_value_t = 10_000)]
    mask_count: usize,
    #[arg(long, default_value_t = 15000.0)]
    mean_length: f64,
    #[arg(long, default_value_t = 3000.0)]
    sd_length: f64,
    #[arg(long, default_value_t = 0.001)]
    error_rate: f64,
    #[arg(long, default_value_t = 42)]
    seed: u64,
    #[arg(long, default_value_t = 4)]
    threads: usize,
    #[arg(long, default_value_t = 100.0)]
    min_quality: f64,
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
        ((self.next() >> 11) as f64 + 0.5) / ((1u64 << 53) as f64)
    }
    fn index(&mut self, n: usize) -> usize {
        (self.next() % (n as u64)) as usize
    }
    fn normal(&mut self) -> f64 {
        (-2.0 * self.uniform().ln()).sqrt() * (std::f64::consts::TAU * self.uniform()).cos()
    }
}

fn record(writer: &mut impl Write, header: &[u8], seq: &[u8]) -> Result<()> {
    writer.write_all(b">")?;
    writer.write_all(header)?;
    writer.write_all(b"\n")?;
    for line in seq.chunks(80) {
        writer.write_all(line)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let options = Options::parse();
    ensure!(
        options.seed != 0 && options.coverage > 0 && options.threads > 0,
        "seed, coverage and threads must be positive"
    );
    ensure!(
        options.error_rate.is_finite() && (0.0..=0.1).contains(&options.error_rate),
        "error-rate must be in [0,0.1]"
    );
    ensure!(
        options.mean_length.is_finite()
            && options.sd_length.is_finite()
            && options.mean_length >= 1000.0
            && options.sd_length >= 0.0,
        "Invalid length distribution"
    );
    let records = input::scaffolds(&options.reference)?;
    ensure!(
        records.len() == 1,
        "Supply a single circular bacterial chromosome"
    );
    let reference = &records[0];
    let sequence = reference.sequence.to_ascii_uppercase();
    ensure!(
        sequence.iter().all(|&b| is_base(b)),
        "Reference contains ambiguous bases"
    );
    ensure!(
        options.mask_count > 0 && options.mask_count < sequence.len(),
        "Invalid mask count"
    );
    fs::create_dir(&options.directory)?;
    let genome_len = sequence.len();
    let mut mask_rng = Rng(options.seed ^ 0xd1b54a32d192ed03);
    let mut positions = FxHashSet::default();
    while positions.len() < options.mask_count {
        positions.insert(mask_rng.index(genome_len));
    }
    let mut positions: Vec<_> = positions.into_iter().collect();
    positions.sort_unstable();
    let mut masked = sequence.clone();
    let mut masks = BufWriter::new(File::create(options.directory.join("masks.tsv"))?);
    writeln!(masks, "position_0based\ttruth_base")?;
    for &pos in &positions {
        masked[pos] = b'N';
        writeln!(masks, "{pos}\t{}", sequence[pos] as char)?;
    }
    masks.flush()?;
    record(
        &mut BufWriter::new(File::create(options.directory.join("masked.fa"))?),
        &reference.header,
        &masked,
    )?;
    let gap_records = input::fasta_bytes(&[b">masked\n".as_slice(), &masked, b"\n"].concat())?;
    let gaps = discover(&gap_records, 31);
    // Reads are simulated from the original, unmasked sequence. Masks and reads
    // use independent PRNG streams so mask placement cannot favor read coverage.
    let mut rng = Rng(options.seed);
    let mut reads = BufWriter::new(File::create(options.directory.join("reads.fa"))?);
    let mut origins = BufWriter::new(File::create(options.directory.join("read_origins.tsv"))?);
    writeln!(
        origins,
        "read_id\tstart_0based\ttemplate_length\treverse\tsubstitutions\tinsertions\tdeletions"
    )?;
    let target = genome_len
        .checked_mul(options.coverage)
        .context("Coverage overflows usize")?;
    let mut consumed = 0;
    let mut read_count = 0;
    let mut read_bases = 0;
    let mut total_errors = [0usize; 3];
    let generation = Instant::now();
    while consumed < target {
        let length = (options.mean_length + options.sd_length * rng.normal())
            .round()
            .clamp(1000.0, options.mean_length * 2.0) as usize;
        let length = length.min(target - consumed);
        let start = rng.index(genome_len);
        let mut read = Vec::with_capacity(length + length / 100);
        let mut errors = [0usize; 3];
        for offset in 0..length {
            let base = sequence[(start + offset) % genome_len];
            if rng.uniform() < options.error_rate {
                match rng.index(3) {
                    0 => {
                        // uniformly select a different base
                        let current = b"ACGT".iter().position(|&b| b == base).unwrap();
                        read.push(b"ACGT"[(current + 1 + rng.index(3)) % 4]);
                        errors[0] += 1;
                    }
                    1 => {
                        read.push(b"ACGT"[rng.index(4)]);
                        read.push(base);
                        errors[1] += 1;
                    }
                    _ => {
                        errors[2] += 1;
                    }
                }
            } else {
                read.push(base);
            }
        }
        let reverse = rng.index(2) == 1;
        if reverse {
            read = revcomp(&read);
        }
        let id = format!("sim{read_count}");
        record(&mut reads, id.as_bytes(), &read)?;
        writeln!(
            origins,
            "{id}\t{start}\t{length}\t{reverse}\t{}\t{}\t{}",
            errors[0], errors[1], errors[2]
        )?;
        for i in 0..3 {
            total_errors[i] += errors[i];
        }
        read_bases += read.len();
        consumed += length;
        read_count += 1;
    }
    reads.flush()?;
    origins.flush()?;
    let mut parameters = File::create(options.directory.join("parameters.txt"))?;
    writeln!(
        parameters,
        "{options:#?}\nreference_header={}\ngenome_length={genome_len}\nmodel=independent Q30-like errors, equal S/I/D; truncated normal lengths; uniform circular starts; random strands\ncalibrated_hifi_model=false\ntemplate_bases={consumed}\nread_bases={read_bases}\nread_count={read_count}\nsubstitutions={}\ninsertions={}\ndeletions={}\ngeneration_seconds={:.6}",
        String::from_utf8_lossy(&reference.header),
        total_errors[0],
        total_errors[1],
        total_errors[2],
        generation.elapsed().as_secs_f64()
    )?;
    let executable = std::env::current_exe()?
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("gapfil");
    let rss = options.directory.join("peak_rss_kib.txt");
    let start = Instant::now();
    let process = Command::new("/usr/bin/time")
        .args(["-f", "%M", "-o"])
        .arg(&rss)
        .arg(executable)
        .arg("--scaffolds")
        .arg(options.directory.join("masked.fa"))
        .arg("--reads")
        .arg(options.directory.join("reads.fa"))
        .arg("--output")
        .arg(options.directory.join("corrected.fa"))
        .arg("--report")
        .arg(options.directory.join("gaps.tsv"))
        .args([
            "--preset",
            "hifi",
            "--k",
            "31",
            "--min-reads",
            "5",
            "--min-quality",
            &options.min_quality.to_string(),
            "--threads",
            &options.threads.to_string(),
        ])
        .output()?;
    let seconds = start.elapsed().as_secs_f64();
    fs::write(options.directory.join("gapfil.stderr.log"), &process.stderr)?;
    ensure!(
        process.status.success(),
        "Gapfil failed: {}",
        String::from_utf8_lossy(&process.stderr)
    );
    let corrected = input::scaffolds(&options.directory.join("corrected.fa"))?;
    ensure!(
        corrected.len() == 1 && corrected[0].header == reference.header,
        "Output header/record count changed"
    );
    let corrected = &corrected[0].sequence;
    let report = fs::read_to_string(options.directory.join("gaps.tsv"))?;
    let mut rows = report.lines().filter(|l| !l.starts_with('#')).skip(1);
    let mut correct = 0;
    let mut wrong = 0;
    let mut unchanged = 0;
    let mut exact_gaps = 0;
    let mut wrong_gaps = 0;
    let mut input_pos = 0;
    let mut output_pos = 0;
    let mut reasons = BTreeMap::new();
    for gap in &gaps {
        let row = rows.next().context("Missing report row")?;
        let fields: Vec<_> = row.split('\t').collect();
        ensure!(
            fields[3].parse::<usize>()? == gap.start && fields[4].parse::<usize>()? == gap.end,
            "Report coordinates disagree"
        );
        let prefix = gap.start - input_pos;
        ensure!(
            corrected.get(output_pos..output_pos + prefix) == Some(&sequence[input_pos..gap.start]),
            "Unmasked sequence was changed"
        );
        output_pos += prefix;
        let filled = fields[15] == "filled";
        let length = if filled {
            fields[11].parse::<usize>()?
        } else {
            gap.end - gap.start
        };
        let observed = corrected
            .get(output_pos..output_pos + length)
            .context("Truncated output")?;
        *reasons.entry(fields[16].to_string()).or_insert(0usize) += 1;
        if filled {
            if observed == &sequence[gap.start..gap.end] {
                exact_gaps += 1;
            } else {
                wrong_gaps += 1;
            }
            if length == gap.end - gap.start {
                for (&base, &truth) in observed.iter().zip(&sequence[gap.start..gap.end]) {
                    if base == truth {
                        correct += 1;
                    } else {
                        wrong += 1;
                    }
                }
            } else {
                wrong += gap.end - gap.start;
            }
        } else {
            ensure!(
                observed.iter().all(|&b| b == b'N'),
                "Rejected gap was modified"
            );
            unchanged += gap.end - gap.start;
        }
        input_pos = gap.end;
        output_pos += length;
    }
    ensure!(rows.next().is_none(), "Extra report rows");
    ensure!(
        corrected.get(output_pos..) == Some(&sequence[input_pos..]),
        "Unmasked suffix was changed"
    );
    ensure!(
        correct + wrong + unchanged == options.mask_count,
        "Masked-base accounting mismatch"
    );
    let peak = fs::read_to_string(rss)?;
    let mut summary = File::create(options.directory.join("summary.tsv"))?;
    let header = "genome_length\tcoverage\tread_count\tread_bases\tmasked_bases\tgap_runs\texact_filled_gaps\twrong_filled_gaps\tcorrect_recovered_bases\twrong_recovered_bases\tunchanged_masked_bases\trecovery_percent\tfill_precision_percent\twall_seconds\tmegabases_per_second\tpeak_rss_kib";
    let precision = if correct + wrong == 0 {
        0.0
    } else {
        100.0 * correct as f64 / (correct + wrong) as f64
    };
    let row = format!(
        "{genome_len}\t{}\t{read_count}\t{read_bases}\t{}\t{}\t{exact_gaps}\t{wrong_gaps}\t{correct}\t{wrong}\t{unchanged}\t{:.6}\t{precision:.6}\t{seconds:.6}\t{:.3}\t{}",
        options.coverage,
        options.mask_count,
        gaps.len(),
        100.0 * correct as f64 / options.mask_count as f64,
        read_bases as f64 / seconds / 1e6,
        peak.trim()
    );
    writeln!(summary, "{header}\n{row}")?;
    println!("{header}\n{row}");
    let mut rejected = File::create(options.directory.join("reasons.tsv"))?;
    writeln!(rejected, "reason\tgap_runs")?;
    for (reason, count) in reasons {
        writeln!(rejected, "{reason}\t{count}")?;
    }
    Ok(())
}
