//! Re-evaluate existing single-reference masking outputs without regenerating reads.
use anyhow::{Context, Result, ensure};
use clap::Parser;
use std::{fs, path::PathBuf};
use yagaf::{discover, input};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    reference: PathBuf,
    #[arg(long)]
    masked: PathBuf,
    #[arg(long)]
    corrected: PathBuf,
    #[arg(long)]
    report: PathBuf,
}
fn main() -> Result<()> {
    let options = Options::parse();
    let reference = input::scaffolds(&options.reference)?;
    let masked = input::scaffolds(&options.masked)?;
    let corrected = input::scaffolds(&options.corrected)?;
    ensure!(
        reference.len() == 1 && masked.len() == 1 && corrected.len() == 1,
        "Expected one reference record"
    );
    ensure!(
        reference[0].header == masked[0].header && reference[0].header == corrected[0].header,
        "Headers changed"
    );
    let truth = &reference[0].sequence;
    let masked_sequence = &masked[0].sequence;
    let output = &corrected[0].sequence;
    ensure!(
        truth.len() == masked_sequence.len(),
        "Masked input length changed"
    );
    ensure!(
        truth
            .iter()
            .zip(masked_sequence)
            .all(|(a, b)| b.eq_ignore_ascii_case(&b'N') || a == b),
        "Input has changes outside masks"
    );
    let gaps = discover(&masked, 31);
    let report = fs::read_to_string(options.report)?;
    let mut rows = report.lines().filter(|l| !l.starts_with('#')).skip(1);
    let (mut input_pos, mut output_pos) = (0, 0);
    let (mut exact_gaps, mut wrong_gaps, mut correct, mut wrong, mut unchanged) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    for gap in &gaps {
        let fields: Vec<_> = rows
            .next()
            .context("Missing report row")?
            .split('\t')
            .collect();
        ensure!(fields.len() >= 17, "Unexpected report columns");
        ensure!(
            fields[3].parse::<usize>()? == gap.start && fields[4].parse::<usize>()? == gap.end,
            "Report coordinates differ"
        );
        let prefix = gap.start - input_pos;
        ensure!(
            output.get(output_pos..output_pos + prefix) == Some(&truth[input_pos..gap.start]),
            "Unmasked bases changed"
        );
        output_pos += prefix;
        let filled = fields[15] == "filled";
        let length = if filled {
            fields[11].parse()?
        } else {
            gap.end - gap.start
        };
        let observed = output
            .get(output_pos..output_pos + length)
            .context("Truncated output")?;
        if filled {
            if observed == &truth[gap.start..gap.end] {
                exact_gaps += 1;
            } else {
                wrong_gaps += 1;
                eprintln!(
                    "Incorrect fill at {}..{}: observed {} bases, expected {}",
                    gap.start,
                    gap.end,
                    length,
                    gap.end - gap.start
                );
            }
            if length == gap.end - gap.start {
                for (a, b) in observed.iter().zip(&truth[gap.start..gap.end]) {
                    if a == b {
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
                observed == &masked_sequence[gap.start..gap.end],
                "Rejected gap changed"
            );
            unchanged += gap.end - gap.start;
        }
        input_pos = gap.end;
        output_pos += length;
    }
    ensure!(rows.next().is_none(), "Extra report rows");
    ensure!(
        output.get(output_pos..) == Some(&truth[input_pos..]),
        "Unmasked suffix changed"
    );
    let total: usize = gaps.iter().map(|g| g.end - g.start).sum();
    ensure!(
        correct + wrong + unchanged == total,
        "Mask accounting differs"
    );
    println!(
        "masked_bases\tgap_runs\texact_filled_gaps\twrong_filled_gaps\tcorrect_recovered_bases\twrong_recovered_bases\tunchanged_masked_bases\trecovery_percent"
    );
    println!(
        "{total}\t{}\t{exact_gaps}\t{wrong_gaps}\t{correct}\t{wrong}\t{unchanged}\t{:.6}",
        gaps.len(),
        100.0 * correct as f64 / total as f64
    );
    Ok(())
}
