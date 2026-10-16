use clap::Parser;
use std::{fs, io::Write};
use tempfile::TempDir;
use yagaf::{
    Scaffold, discover,
    index::{Index, Scanner, Span},
    input,
    pipeline::{Args, run},
    revcomp,
};

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

fn joined(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

struct Fixture {
    dir: TempDir,
    left: Vec<u8>,
    right: Vec<u8>,
}

impl Fixture {
    fn new() -> Self {
        let mut seed = 0x123456789abcdef;
        let f = Self {
            dir: tempfile::tempdir().unwrap(),
            left: dna(&mut seed, 23),
            right: dna(&mut seed, 23),
        };
        f.write_scaffolds(&[joined(&[&f.left, b"nnn", &f.right])]);
        f
    }
    fn write_scaffolds(&self, records: &[Vec<u8>]) {
        let mut file = fs::File::create(self.dir.path().join("scaffolds.fa")).unwrap();
        for (i, seq) in records.iter().enumerate() {
            writeln!(
                file,
                ">s{i} full header {i}\n{}",
                String::from_utf8_lossy(seq)
            )
            .unwrap();
        }
    }
    fn span(&self, middle: &[u8]) -> Vec<u8> {
        joined(&[&self.left, middle, &self.right])
    }
    fn reads(&self, middles: &[&[u8]]) {
        let reads: Vec<_> = middles
            .iter()
            .enumerate()
            .map(|(i, s)| (format!("r{i}"), self.span(s)))
            .collect();
        self.write_reads(&reads);
    }
    fn write_reads(&self, reads: &[(String, Vec<u8>)]) {
        let mut file = fs::File::create(self.dir.path().join("reads.fa")).unwrap();
        for (id, seq) in reads {
            writeln!(file, ">{id} description\n{}", String::from_utf8_lossy(seq)).unwrap();
        }
    }
    fn args(&self, name: &str, extra: &[&str]) -> Args {
        let mut values: Vec<String> = vec![
            "yagaf".into(),
            "--scaffolds".into(),
            self.dir.path().join("scaffolds.fa").display().to_string(),
            "--reads".into(),
            self.dir.path().join("reads.fa").display().to_string(),
            "--output".into(),
            self.dir
                .path()
                .join(format!("{name}.fa"))
                .display()
                .to_string(),
            "--report".into(),
            self.dir
                .path()
                .join(format!("{name}.tsv"))
                .display()
                .to_string(),
            "--preset".into(),
            "hifi".into(),
            "--k".into(),
            "23".into(),
        ];
        // Allow tests to replace any configurable value without duplicate flags.
        for (flag, default) in [
            ("--min-reads", "3"),
            ("--min-quality", "100"),
            ("--threads", "2"),
        ] {
            if !extra.contains(&flag) {
                values.extend([flag.into(), default.into()]);
            }
        }
        values.extend(extra.iter().map(|s| s.to_string()));
        Args::parse_from(values)
    }
    fn output(&self, name: &str) -> Vec<Scaffold> {
        input::scaffolds(&self.dir.path().join(format!("{name}.fa"))).unwrap()
    }
    fn report(&self, name: &str) -> String {
        fs::read_to_string(self.dir.path().join(format!("{name}.tsv"))).unwrap()
    }
}

#[test]
fn exact_and_reverse_complement_reads_fill_and_preserve_headers_and_case() {
    let f = Fixture::new();
    f.write_scaffolds(&[joined(&[&f.left.to_ascii_lowercase(), b"nnn", &f.right])]);
    let seq = f.span(b"ACTGAC");
    f.write_reads(&[
        ("forward".into(), seq.clone()),
        ("reverse".into(), revcomp(&seq)),
        ("lower".into(), seq.to_ascii_lowercase()),
    ]);
    run(&f.args("out", &[])).unwrap();
    let out = f.output("out");
    assert_eq!(out[0].header, b"s0 full header 0");
    assert_eq!(
        out[0].sequence,
        joined(&[&f.left.to_ascii_lowercase(), b"ACTGAC", &f.right])
    );
    assert!(f.report("out").contains("\t100.00000000\tfilled\taccepted"));
}

#[test]
fn quality_threshold_is_inclusive_and_does_not_count_anchor_context() {
    let f = Fixture::new();
    f.reads(&[b"AACCGGTT", b"AACCGGTT", b"AACCGGTA"]);
    run(&f.args("pass", &["--min-quality", "87.5"])).unwrap();
    run(&f.args("fail", &["--min-quality", "87.500001"])).unwrap();
    assert_eq!(f.output("pass")[0].sequence, f.span(b"AACCGGTT"));
    assert_eq!(f.output("fail")[0].sequence, f.span(b"nnn"));
    assert!(f.report("pass").contains("\t7\t8\t87.50000000\tfilled"));
    assert!(f.report("fail").contains("low_msa_quality"));
}

#[test]
fn indel_msa_uses_all_columns_and_ties_are_rejected() {
    let f = Fixture::new();
    f.reads(&[b"ACGTAC", b"ACTGTAC", b"ACTGTAC"]);
    run(&f.args("indel", &["--min-quality", "85"])).unwrap();
    assert_eq!(f.output("indel")[0].sequence, f.span(b"ACTGTAC"));
    assert!(f.report("indel").contains("\t6\t7\t85.71428571\tfilled"));
    f.reads(&[b"ACGTAC", b"ACGTTC"]);
    run(&f.args("tie", &["--min-reads", "2", "--min-quality", "1"])).unwrap();
    assert_eq!(f.output("tie")[0].sequence, f.span(b"nnn"));
    assert!(f.report("tie").contains("consensus_tie"));
}

#[test]
fn repeated_read_ids_do_not_inflate_support() {
    let f = Fixture::new();
    f.write_reads(&[
        ("r".into(), f.span(b"ACTG")),
        ("r".into(), f.span(b"ACTG")),
        ("r".into(), f.span(b"ACTG")),
    ]);
    run(&f.args("duplicate", &[])).unwrap();
    assert!(f.report("duplicate").contains("insufficient_support"));
    assert_eq!(f.output("duplicate")[0].sequence, f.span(b"nnn"));
    f.write_reads(&[
        ("r".into(), f.span(b"ACTG")),
        ("other".into(), f.span(b"ACTG")),
        ("r".into(), f.span(b"ACCG")),
    ]);
    run(&f.args("conflict", &["--min-reads", "2"])).unwrap();
    assert!(
        f.report("conflict")
            .contains("conflicting_duplicate_read_id")
    );
}

#[test]
fn repeated_scaffold_anchors_are_rejected() {
    let f = Fixture::new();
    f.write_scaffolds(&[f.span(b"NNN"), f.left.clone()]);
    f.reads(&[b"ACTG", b"ACTG", b"ACTG"]);
    run(&f.args("repeat", &[])).unwrap();
    assert!(f.report("repeat").contains("nonunique_scaffold_anchor"));
    assert_eq!(f.output("repeat")[0].sequence, f.span(b"NNN"));
}

#[test]
fn repeated_read_anchors_and_wrong_order_do_not_support_a_gap() {
    let f = Fixture::new();
    f.write_reads(&[
        ("repeat".into(), joined(&[&f.left, &f.span(b"ACTG")])),
        ("wrong".into(), joined(&[&f.right, b"ACTG", &f.left])),
        ("partial".into(), f.left.clone()),
    ]);
    run(&f.args("out", &["--min-reads", "2"])).unwrap();
    assert!(f.report("out").contains("insufficient_support"));
    assert_eq!(f.output("out")[0].sequence, f.span(b"nnn"));
}

#[test]
fn ambiguous_spans_are_not_silently_removed_from_evidence() {
    let f = Fixture::new();
    f.reads(&[b"ACTG", b"ACTG", b"ACTG", b"ACRG"]);
    run(&f.args("out", &[])).unwrap();
    assert!(f.report("out").contains("ambiguous_read_span"));
    assert_eq!(f.output("out")[0].sequence, f.span(b"nnn"));
}

#[test]
fn terminal_and_adjacent_gaps_keep_original_sequence() {
    let f = Fixture::new();
    let seq = joined(&[b"nn", &f.left, b"NN", b"A", b"n", &f.right, b"N"]);
    f.write_scaffolds(std::slice::from_ref(&seq));
    f.reads(&[b"ACTG"]);
    run(&f.args("out", &[])).unwrap();
    assert_eq!(f.output("out")[0].sequence, seq);
    assert_eq!(
        f.report("out")
            .matches("missing_or_ambiguous_flank")
            .count(),
        4
    );
}

#[test]
fn resource_limits_reject_whole_gaps() {
    let f = Fixture::new();
    f.reads(&[b"ACTG", b"ACTG", b"ACTG", b"ACTG"]);
    run(&f.args("cap", &["--max-reads-per-gap", "3"])).unwrap();
    assert!(f.report("cap").contains("read_support_limit"));
    assert_eq!(f.output("cap")[0].sequence, f.span(b"nnn"));
    f.reads(&[b"ACTG", b"ACTG", b"ACTG", b"ACTGACTGACTG"]);
    run(&f.args("span", &["--max-span", "8"])).unwrap();
    assert!(f.report("span").contains("span_limit"));
    f.reads(&[b"ACTG", b"ACCG", b"ACTG"]);
    run(&f.args("msa", &["--max-msa-cells", "1"])).unwrap();
    assert!(f.report("msa").contains("msa_size_limit"));
    run(&f.args("memory", &["--max-evidence-bytes", "5"])).unwrap();
    assert!(f.report("memory").contains("global_evidence_limit"));
}

#[test]
fn forty_gaps_and_length_changing_replacements_are_deterministic() {
    let f = Fixture::new();
    let mut seed = 9321805748291;
    let mut scaffold = Vec::new();
    let mut expected = Vec::new();
    let mut reads = Vec::new();
    for id in 0..40 {
        let left = dna(&mut seed, 23);
        let right = dna(&mut seed, 23);
        let middle = dna(&mut seed, 4 + id % 11);
        scaffold.extend(joined(&[&left, b"NNN", &right]));
        expected.extend(joined(&[&left, &middle, &right]));
        let seq = joined(&[&left, &middle, &right]);
        for replicate in 0..3 {
            reads.push((
                format!("r{id}_{replicate}"),
                if replicate == 1 {
                    revcomp(&seq)
                } else {
                    seq.clone()
                },
            ));
        }
    }
    f.write_scaffolds(&[scaffold]);
    f.write_reads(&reads);
    run(&f.args("one", &["--threads", "1", "--batch-bases", "47"])).unwrap();
    run(&f.args("four", &["--threads", "4", "--batch-bases", "113"])).unwrap();
    assert_eq!(f.output("one")[0].sequence, expected);
    assert_eq!(f.output("four")[0].sequence, expected);
    assert_eq!(f.report("one"), f.report("four"));
    assert_eq!(f.report("one").matches("\tfilled\taccepted").count(), 40);
}

#[test]
fn parallel_merge_preserves_capacity_decisions_and_read_order() {
    let f = Fixture::new();
    let mut seed = 21934572581;
    let mut masked = Vec::new();
    let mut truth = Vec::new();
    for _ in 0..40 {
        let left = dna(&mut seed, 23);
        let right = dna(&mut seed, 23);
        masked.extend(joined(&[&left, b"NNN", &right]));
        truth.extend(joined(&[&left, b"ACTG", &right]));
    }
    f.write_scaffolds(&[masked]);
    f.write_reads(
        &(0..128)
            .map(|i| (format!("r{i}"), truth.clone()))
            .collect::<Vec<_>>(),
    );
    for (label, limits) in [
        ("room", vec![]),
        ("tight", vec!["--max-evidence-bytes", "1000"]),
    ] {
        let one = format!("one-{label}");
        let many = format!("many-{label}");
        let mut a = vec!["--threads", "1"];
        a.extend(&limits);
        let mut b = vec!["--threads", "32"];
        b.extend(&limits);
        run(&f.args(&one, &a)).unwrap();
        run(&f.args(&many, &b)).unwrap();
        assert_eq!(f.output(&one)[0].sequence, f.output(&many)[0].sequence);
        assert_eq!(f.report(&one), f.report(&many));
        if label == "room" {
            assert_eq!(f.output(&many)[0].sequence, truth);
        }
    }
}

#[test]
fn compressed_fasta_and_fastq_are_parsed_by_helicase() {
    let f = Fixture::new();
    use deko::{
        Format,
        write::{AnyEncoder, Compression},
    };
    let seq = f.span(b"ACTG");
    let scaffold = joined(&[b">s0 full header 0\n", &f.span(b"nnn"), b"\n"]);
    for format in [
        Format::Verbatim,
        Format::Gz,
        Format::Zstd,
        Format::Bz,
        Format::Xz,
    ] {
        let encode = |bytes: &[u8]| {
            let mut encoder = AnyEncoder::new(Vec::new(), format, Compression::Fast).unwrap();
            encoder.write_all(bytes).unwrap();
            encoder.finish().unwrap()
        };
        // Deliberately use opaque filenames: decoding must use content magic.
        let scaffolds = f.dir.path().join(format!("scaffolds-{format:?}.data"));
        fs::write(&scaffolds, encode(&scaffold)).unwrap();
        for fastq in [false, true] {
            let mut reads = Vec::new();
            for i in 0..3 {
                if fastq {
                    writeln!(
                        &mut reads,
                        "@q{i}\n{}\n+\n{}",
                        String::from_utf8_lossy(&seq),
                        "I".repeat(seq.len())
                    )
                    .unwrap();
                } else {
                    writeln!(&mut reads, ">r{i}\n{}", String::from_utf8_lossy(&seq)).unwrap();
                }
            }
            let path = f.dir.path().join(format!("reads-{format:?}-{fastq}.data"));
            fs::write(&path, encode(&reads)).unwrap();
            let name = format!("{format:?}-{fastq}");
            let mut args = f.args(&name, &[]);
            args.scaffolds = scaffolds.clone();
            args.reads = path;
            run(&args).unwrap();
            assert_eq!(f.output(&name)[0].sequence, seq);
        }
    }
}

#[test]
fn invalid_arguments_and_existing_outputs_fail_without_writing() {
    let f = Fixture::new();
    f.reads(&[b"ACTG", b"ACTG", b"ACTG"]);
    for extra in [
        vec!["--min-quality", "NaN"],
        vec!["--min-quality", "101"],
        vec!["--min-reads", "1"],
        vec!["--threads", "0"],
        vec!["--minimizer", "24"],
    ] {
        assert!(run(&f.args("bad", &extra)).is_err());
        assert!(!f.dir.path().join("bad.fa").exists());
    }
    run(&f.args("exists", &[])).unwrap();
    assert!(run(&f.args("exists", &[])).is_err());
    let mut args = f.args("same", &[]);
    args.report = args.output.clone();
    assert!(run(&args).is_err());
    let mut args = f.args("input", &[]);
    args.output = args.scaffolds.clone();
    assert!(run(&args).is_err());
}

#[test]
fn malformed_fastq_and_oversized_reads_fail_before_outputs() {
    let f = Fixture::new();
    fs::write(f.dir.path().join("bad.fq"), b"@r\nACGT\n+\n!!\n").unwrap();
    let mut args = f.args("badq", &[]);
    args.reads = f.dir.path().join("bad.fq");
    assert!(run(&args).is_err());
    assert!(!args.output.exists() && !args.report.exists());
    f.reads(&[b"ACTG", b"ACTG", b"ACTG"]);
    let args = f.args("oversized", &["--max-read-bases", "40", "--max-span", "8"]);
    assert!(run(&args).is_err());
    assert!(!args.output.exists() && !args.report.exists());
}

#[test]
fn a_genome_without_eligible_gaps_does_not_scan_reads() {
    let f = Fixture::new();
    f.write_scaffolds(&[f.span(b"ACTG")]);
    // No read file is created: there is no gap evidence to collect.
    run(&f.args("nogaps", &[])).unwrap();
    assert_eq!(f.output("nogaps")[0].sequence, f.span(b"ACTG"));
}

#[test]
fn simd_filter_matches_exact_search_across_k_sizes_and_ambiguous_runs() {
    let mut seed = 0x123456789abcdef;
    for k in [3, 5, 17, 31, 33, 63] {
        let left = dna(&mut seed, k);
        let right = dna(&mut seed, k);
        let scaffold = Scaffold {
            header: b"s".to_vec(),
            sequence: joined(&[&left, b"NNN", &right]),
        };
        let gaps = discover(&[scaffold], k);
        for m in [1, k.min(15), k.min(29)] {
            let index = Index::new(&gaps, k, m);
            let middle = dna(&mut seed, 19);
            let read = joined(&[b"NNR", &left, &middle, &right, b"NRY"]);
            // Some small random k-mers can repeat in the span; compare exact
            // matching counts and orientation instead of assuming uniqueness.
            for read in [read.clone(), revcomp(&read)] {
                let forward = exact_interval(&read, &left, &right);
                let reverse = exact_interval(&read, &revcomp(&right), &revcomp(&left));
                let found = index.extract(&read, &mut Scanner::default(), 1000);
                if let Some((start, end)) = forward.filter(|_| reverse.is_none()) {
                    assert!(
                        matches!(&found[..],[(0,Span::Sequence(seq))] if *seq==read[start..end])
                    );
                } else if let Some((start, end)) = reverse.filter(|_| forward.is_none()) {
                    assert!(
                        matches!(&found[..],[(0,Span::Sequence(seq))] if *seq==revcomp(&read[start..end]))
                    );
                }
            }
            // Inserting an ambiguity into an anchor must never reconstruct the
            // original anchor by concatenating the valid runs.
            if k >= 17 {
                let broken = joined(&[&left[..k / 2], b"N", &left[k / 2..], &middle, &right]);
                assert!(
                    index
                        .extract(&broken, &mut Scanner::default(), 1000)
                        .is_empty()
                );
            }
        }
    }
}

fn exact_interval(read: &[u8], left: &[u8], right: &[u8]) -> Option<(usize, usize)> {
    let l: Vec<_> = read
        .windows(left.len())
        .enumerate()
        .filter(|(_, w)| *w == left)
        .map(|(p, _)| p)
        .collect();
    let r: Vec<_> = read
        .windows(right.len())
        .enumerate()
        .filter(|(_, w)| *w == right)
        .map(|(p, _)| p)
        .collect();
    (l.len() == 1 && r.len() == 1 && l[0] + left.len() < r[0]).then(|| (l[0] + left.len(), r[0]))
}
