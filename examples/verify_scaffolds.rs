use std::{env, fs};
use yagaf::{discover, input, is_base};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().collect();
    if args.len() != 4 {
        return Err("Usage: verify_output INPUT_FASTA OUTPUT_FASTA GAP_REPORT".into());
    }
    let original = input::scaffolds(std::path::Path::new(&args[1]))?;
    let output = input::scaffolds(std::path::Path::new(&args[2]))?;
    assert_eq!(original.len(), output.len(), "Scaffold count changed");
    let gaps = discover(&original, 31);
    let report = fs::read_to_string(&args[3])?;
    let parameters = report.lines().next().ok_or("Missing parameters")?;
    let min_reads: usize = parameters
        .split_ascii_whitespace()
        .find_map(|s| s.strip_prefix("min_reads="))
        .ok_or("Missing read threshold")?
        .parse()?;
    let min_quality: f64 = parameters
        .split_ascii_whitespace()
        .find_map(|s| s.strip_prefix("min_quality_percent="))
        .ok_or("Missing quality threshold")?
        .parse()?;
    let rows: Vec<Vec<&str>> = report
        .lines()
        .filter(|l| !l.starts_with('#'))
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), gaps.len(), "Report gap count differs");
    let (mut cursor, mut filled, mut removed_ns, mut inserted_bases, mut changed_scaffolds) =
        (0, 0, 0usize, 0usize, 0);
    for (scaffold_id, (before, after)) in original.iter().zip(&output).enumerate() {
        assert_eq!(before.header, after.header, "Header changed");
        let (mut in_pos, mut out_pos, mut changed) = (0, 0, false);
        while cursor < gaps.len() && gaps[cursor].scaffold == scaffold_id {
            let gap = &gaps[cursor];
            let row = &rows[cursor];
            assert!(row.len() >= 17, "Unexpected report width");
            assert_eq!(row[0].parse::<usize>()?, cursor);
            assert_eq!(row[1].parse::<usize>()?, scaffold_id);
            assert_eq!(row[3].parse::<usize>()?, gap.start);
            assert_eq!(row[4].parse::<usize>()?, gap.end);
            assert_eq!(row[5].parse::<usize>()?, gap.end - gap.start);
            let prefix = gap.start - in_pos;
            assert_eq!(
                after.sequence.get(out_pos..out_pos + prefix),
                Some(&before.sequence[in_pos..gap.start]),
                "Bases outside a gap changed"
            );
            out_pos += prefix;
            match row[15] {
                "filled" => {
                    let support = row[10].parse::<usize>()?;
                    let quality = row[14].parse::<f64>()?;
                    assert!(support >= min_reads && (min_quality..=100.0).contains(&quality));
                    assert_eq!(row[16], "accepted");
                    let length = row[11].parse::<usize>()?;
                    assert!(length > 0);
                    let seq = after
                        .sequence
                        .get(out_pos..out_pos + length)
                        .ok_or("Truncated fill")?;
                    assert!(seq.iter().all(|&b| is_base(b)), "Invalid fill bases");
                    filled += 1;
                    removed_ns += gap.end - gap.start;
                    inserted_bases += length;
                    out_pos += length;
                    changed = true;
                }
                "unchanged" => {
                    let length = gap.end - gap.start;
                    assert_eq!(
                        after.sequence.get(out_pos..out_pos + length),
                        Some(&before.sequence[gap.start..gap.end]),
                        "Rejected gap changed"
                    );
                    out_pos += length;
                }
                _ => return Err("Unknown gap status".into()),
            }
            in_pos = gap.end;
            cursor += 1;
        }
        assert_eq!(
            after.sequence.get(out_pos..),
            Some(&before.sequence[in_pos..]),
            "Scaffold suffix changed"
        );
        changed_scaffolds += usize::from(changed);
    }
    assert_eq!(cursor, gaps.len());
    let original_bases: usize = original.iter().map(|s| s.sequence.len()).sum();
    let output_bases: usize = output.iter().map(|s| s.sequence.len()).sum();
    let original_ns: usize = gaps.iter().map(|g| g.end - g.start).sum();
    let output_ns: usize = output
        .iter()
        .map(|s| {
            s.sequence
                .iter()
                .filter(|b| b.eq_ignore_ascii_case(&b'N'))
                .count()
        })
        .sum();
    assert_eq!(output_bases, original_bases - removed_ns + inserted_bases);
    assert_eq!(output_ns, original_ns - removed_ns);
    println!("metric\tvalue");
    for (key, value) in [
        ("scaffolds", original.len()),
        ("input_bases", original_bases),
        ("output_bases", output_bases),
        ("gap_runs", gaps.len()),
        ("filled_gaps", filled),
        ("changed_scaffolds", changed_scaffolds),
        ("original_ns", original_ns),
        ("remaining_ns", output_ns),
        ("replaced_ns", removed_ns),
        ("inserted_consensus_bases", inserted_bases),
    ] {
        println!("{key}\t{value}");
    }
    println!("unchanged_regions_verified\ttrue");
    println!("headers_and_record_order_verified\ttrue");
    Ok(())
}
