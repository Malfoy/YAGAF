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
fn repeated_scaffold_anchors_remain_eligible() {
    let f = Fixture::new();
    let repeated = joined(&[&f.left, b"ACGT", &revcomp(&f.right)]);
    f.write_scaffolds(&[f.span(b"NNN"), repeated.clone()]);
    f.reads(&[b"ACTG", b"ACTG", b"ACTG"]);
    run(&f.args("repeat", &[])).unwrap();
    assert!(f.report("repeat").contains("\tfilled\taccepted"));
    assert_eq!(f.output("repeat")[0].sequence, f.span(b"ACTG"));
    assert_eq!(f.output("repeat")[1].sequence, repeated);
}

#[test]
fn shared_assembly_flanks_keep_gap_evidence_separate() {
    let f = Fixture::new();
    let mut seed = 721432194213;
    let other_right = dna(&mut seed, 23);
    f.write_scaffolds(&[f.span(b"NNN"), joined(&[&f.left, b"NNNN", &other_right])]);
    let mut reads = Vec::new();
    for i in 0..3 {
        reads.push((format!("a{i}"), f.span(b"ACTG")));
        reads.push((
            format!("b{i}"),
            revcomp(&joined(&[&f.left, b"GGTTAC", &other_right])),
        ));
    }
    f.write_reads(&reads);
    run(&f.args("shared", &[])).unwrap();
    let output = f.output("shared");
    assert_eq!(output[0].sequence, f.span(b"ACTG"));
    assert_eq!(
        output[1].sequence,
        joined(&[&f.left, b"GGTTAC", &other_right])
    );
    assert_eq!(f.report("shared").matches("\tfilled\taccepted").count(), 2);
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
    args.report = Some(args.output.clone());
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
    assert!(!args.output.exists() && !args.report_path().exists());
    f.reads(&[b"ACTG", b"ACTG", b"ACTG"]);
    let args = f.args("oversized", &["--max-read-bases", "40", "--max-span", "8"]);
    assert!(run(&args).is_err());
    assert!(!args.output.exists() && !args.report_path().exists());
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
                        matches!(&found[..],[(0,Span::Sequence(seq))] if seq.as_slice()==&read[start..end])
                    );
                } else if let Some((start, end)) = reverse.filter(|_| forward.is_none()) {
                    assert!(
                        matches!(&found[..],[(0,Span::Sequence(seq))] if seq.as_slice()==revcomp(&read[start..end]))
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

#[test]
fn shared_defaults_and_automatic_report_work_with_three_required_paths() {
    let f = Fixture::new();
    let mut seed = 819241;
    let left = dna(&mut seed, 31);
    let right = dna(&mut seed, 31);
    f.write_scaffolds(&[joined(&[&left, b"NNN", &right])]);
    let read = joined(&[&left, b"ACTGAC", &right]);
    f.write_reads(&[
        ("a".into(), read.clone()),
        ("b".into(), read.clone()),
        ("c".into(), read.clone()),
    ]);
    let args = Args::parse_from([
        "yagaf".into(),
        "--scaffolds".into(),
        f.dir.path().join("scaffolds.fa").into_os_string(),
        "--reads".into(),
        f.dir.path().join("reads.fa").into_os_string(),
        "--output".into(),
        f.dir.path().join("defaults.fa").into_os_string(),
    ]);
    run(&args).unwrap();
    assert_eq!(f.output("defaults")[0].sequence, read);
    assert_eq!(
        args.report_path(),
        f.dir.path().join("defaults.fa.gaps.tsv")
    );
    let report = fs::read_to_string(args.report_path()).unwrap();
    assert!(report.contains("min_reads=3 min_quality_percent=50"));
    assert!(report.contains("msa_backend=auto"));
    assert!(run(&args).is_err());
    for scoring in ["hifi", "ont"] {
        let a = Args::parse_from([
            "yagaf", "-s", "s.fa", "-r", "r.fa", "-o", "o.fa", "--preset", scoring,
        ]);
        assert_eq!((a.k, a.min_reads, a.min_quality), (31, 3, 50.0));
        assert_eq!(a.max_reads_per_gap, 4096);
        assert_eq!(a.recruitment, yagaf::pipeline::RecruitmentMode::Adaptive);
        assert_eq!(
            (a.flank_window, a.max_flank_edits, a.max_msa_cells),
            (100, Some(1), 268_435_456)
        );
        assert!(matches!(a.msa_backend, yagaf::consensus::MsaBackend::Auto));
    }
}

#[test]
fn mapped_and_streamed_parsers_match_across_record_boundaries() {
    let f = Fixture::new();
    let middle = b"ACGTACGTTGCAACGTTACGATCGATGCTAGCTACGATCGT";
    let mut seed = 721937127397;
    for fastq in [false, true] {
        let path = f.dir.path().join("reads.input");
        let mut file = fs::File::create(&path).unwrap();
        for i in 0..120 {
            let padding = dna(&mut seed, if i == 0 { 80_000 } else { 1024 + i });
            let seq = joined(&[&padding, &f.span(middle), b"NN"]);
            let seq = if i % 2 == 0 { seq } else { revcomp(&seq) };
            writeln!(
                file,
                "{}r{} description > @@++",
                if fastq { '@' } else { '>' },
                i % 100
            )
            .unwrap();
            if fastq {
                writeln!(file, "{}\n+", String::from_utf8_lossy(&seq)).unwrap();
                let qualities: String =
                    (0..seq.len()).map(|p| b"@>+IACGT"[p % 8] as char).collect();
                writeln!(file, "{qualities}").unwrap();
            } else {
                for line in seq.chunks(79) {
                    writeln!(file, "{}", String::from_utf8_lossy(line)).unwrap();
                }
            }
        }
        drop(file);
        let label = if fastq { "fastq" } else { "fasta" };
        let mut mapped = f.args(
            &format!("mapped-{label}"),
            &["--threads", "32", "--batch-bases", "65536"],
        );
        mapped.reads = path.clone();
        let mut streamed = f.args(
            &format!("streamed-{label}"),
            &["--threads", "1", "--stream-reads", "--batch-bases", "65536"],
        );
        streamed.reads = path;
        run(&mapped).unwrap();
        run(&streamed).unwrap();
        assert_eq!(
            f.output(&format!("mapped-{label}"))[0].sequence,
            f.span(middle)
        );
        assert_eq!(
            f.report(&format!("mapped-{label}")),
            f.report(&format!("streamed-{label}"))
        );
        assert!(f.report(&format!("mapped-{label}")).contains("\t20\t100\t"));
    }
}

#[test]
fn flank_windows_recover_through_read_errors_without_copying_known_bases() {
    let mut f = Fixture::new();
    let mut seed = 824732674329;
    f.left = dna(&mut seed, 100);
    f.right = dna(&mut seed, 100);
    f.write_scaffolds(&[f.span(b"nnn")]);
    let mut left = f.left.clone();
    left.insert(83, b'A');
    left[94] = if left[94] == b'A' { b'C' } else { b'A' };
    let mut right = f.right.clone();
    right.remove(12);
    right[5] = if right[5] == b'A' { b'C' } else { b'A' };
    let read = joined(&[&left, b"GATTACA", &right]);
    let scaffolds = input::scaffolds(&f.dir.path().join("scaffolds.fa")).unwrap();
    let gaps = discover(&scaffolds, 23);
    assert!(
        Index::new(&gaps, 23, 15)
            .extract(&read, &mut Scanner::default(), 1000)
            .is_empty()
    );
    let index = Index::with_windows(&scaffolds, &gaps, 23, 15, 100, 11);
    for sequence in [&read, &revcomp(&read)] {
        let spans = index.extract(sequence, &mut Scanner::default(), 1000);
        assert_eq!(spans.len(), 1, "{spans:?}");
        assert!(
            matches!(&spans[0].1, Span::Sequence(s) if s.as_slice() == b"GATTACA"),
            "{spans:?}"
        );
    }
    f.write_reads(&[
        ("a".into(), read.clone()),
        ("b".into(), revcomp(&read)),
        ("c".into(), read),
    ]);
    run(&f.args("windows-one", &["--threads", "1", "--max-flank-edits", "2"])).unwrap();
    run(&f.args(
        "windows-four",
        &["--threads", "4", "--max-flank-edits", "2"],
    ))
    .unwrap();
    assert_eq!(f.output("windows-one")[0].sequence, f.span(b"GATTACA"));
    assert_eq!(f.output("windows-four")[0].sequence, f.span(b"GATTACA"));
    assert_eq!(f.report("windows-one"), f.report("windows-four"));
    let report = f.report("windows-one");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    // Every seed in a read still contributes just one spanning record and read ID.
    assert!(report.contains("\t3\t3\t"), "{row:?}");
}

#[test]
fn flank_context_resolves_repeated_seeds_but_rejects_distinct_intervals() {
    let mut seed = 241876538796;
    let left = dna(&mut seed, 100);
    let right = dna(&mut seed, 100);
    let scaffolds = vec![Scaffold {
        header: b"s".to_vec(),
        sequence: joined(&[&left, b"NNN", &right]),
    }];
    let gaps = discover(&scaffolds, 31);
    let index = Index::with_windows(&scaffolds, &gaps, 31, 15, 100, 15);
    let decoy = joined(&[&dna(&mut seed, 69), &left[69..], &dna(&mut seed, 40)]);
    let read = joined(&[&decoy, &left, b"GATTACA", &right]);
    assert!(matches!(
        Index::new(&gaps, 31, 15).extract(&read, &mut Scanner::default(), 1000)[0].1,
        Span::AmbiguousPair
    ));
    let spans = index.extract(&read, &mut Scanner::default(), 1000);
    assert!(
        matches!(&spans[0].1, Span::Sequence(s) if s.as_slice() == b"GATTACA"),
        "{spans:?}"
    );
    let repeated = joined(&[&left, b"CC", &right, &left, b"GATTACA", &right]);
    assert!(matches!(
        index.extract(&repeated, &mut Scanner::default(), 1000)[0].1,
        Span::AmbiguousPair
    ));
}

#[test]
fn window_arguments_reject_invalid_values() {
    let f = Fixture::new();
    for args in [
        &["--flank-window", "22"][..],
        &["--anchor-step", "0"],
        &["--anchor-step", "24"],
        &["--rescue-flank-window", "22"],
        &["--rescue-flank-window", "101"],
    ] {
        assert!(f.args("invalid-window", args).validate().is_err());
    }
}

#[test]
fn best_context_requires_a_unique_score_advantage_and_preserves_true_repeats() {
    use yagaf::index::Recruitment;
    let mut seed = 242117483579;
    let left = dna(&mut seed, 100);
    let right = dna(&mut seed, 100);
    let scaffolds = vec![Scaffold {
        header: b"s".to_vec(),
        sequence: joined(&[&left, b"NNN", &right]),
    }];
    let gaps = discover(&scaffolds, 31);
    let strict = Index::with_windows(&scaffolds, &gaps, 31, 15, 100, 15);
    let best = Index::with_windows(&scaffolds, &gaps, 31, 15, 100, 15)
        .with_recruitment(Recruitment::BestContext, 2);
    for mismatches in [0, 1, 3] {
        let mut decoy = left.clone();
        for &pos in &[5, 35, 65][..mismatches] {
            decoy[pos] = if decoy[pos] == b'A' { b'C' } else { b'A' };
        }
        let read = joined(&[&decoy, &dna(&mut seed, 40), &left, b"GATTACA", &right]);
        assert!(matches!(
            strict.extract(&read, &mut Scanner::default(), 1000)[0].1,
            Span::AmbiguousPair
        ));
        for oriented in [&read, &revcomp(&read)] {
            let spans = best.extract(oriented, &mut Scanner::default(), 1000);
            if mismatches >= 2 {
                assert!(
                    matches!(&spans[0].1,Span::Sequence(s) if s.as_slice()==b"GATTACA"),
                    "{spans:?}"
                );
            } else {
                assert!(matches!(spans[0].1, Span::AmbiguousPair), "{spans:?}");
            }
        }
    }
}

#[test]
fn length_quality_bound_rejects_impossible_scores_without_discarding_reads() {
    let f = Fixture::new();
    f.reads(&[b"ACTG", b"ACTGA", b"ACTGA"]);
    run(&f.args("bound-inclusive", &["--min-quality", "80"])).unwrap();
    assert_eq!(f.output("bound-inclusive")[0].sequence, f.span(b"ACTGA"));
    assert!(
        f.report("bound-inclusive")
            .contains("\t4\t5\t80.00000000\tfilled")
    );
    run(&f.args(
        "bound-rejected",
        &["--min-quality", "80.000001", "--max-msa-cells", "1"],
    ))
    .unwrap();
    assert_eq!(f.output("bound-rejected")[0].sequence, f.span(b"nnn"));
    assert!(f.report("bound-rejected").contains("msa_quality_bound"));
    let report = f.report("bound-rejected");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    assert_eq!(&row[22..27], &["4", "5", "14", "70", "80.00000000"]);
    assert!(
        f.report("bound-rejected")
            .contains("\t3\tNA\tNA\tNA\tNA\tunchanged")
    );
}

#[test]
fn msa_engines_and_orders_preserve_evidence_and_thread_determinism() {
    let f = Fixture::new();
    f.reads(&[b"ACTGAC", b"ACTTGAC", b"ACTTGAC", b"ACTTGAC"]);
    for engine in ["poasta", "spoars", "spoars-scalar"] {
        for order in ["lexical", "longest-first", "median-first"] {
            let mut previous = None;
            for threads in ["1", "4"] {
                let name = format!("{engine}-{order}-{threads}");
                run(&f.args(
                    &name,
                    &[
                        "--msa-backend",
                        engine,
                        "--msa-order",
                        order,
                        "--threads",
                        threads,
                        "--min-quality",
                        "85",
                    ],
                ))
                .unwrap();
                assert_eq!(f.output(&name)[0].sequence, f.span(b"ACTTGAC"));
                assert!(f.report(&name).contains("\t6\t7\t85.71428571\tfilled"));
                let report = f.report(&name);
                if let Some(previous) = &previous {
                    assert_eq!(&report, previous);
                }
                previous = Some(report);
            }
        }
    }
}

#[test]
fn single_edit_flank_extension_recovers_indels_and_obeys_the_cap() {
    let mut seed = 318764387321;
    let left = dna(&mut seed, 100);
    let right = dna(&mut seed, 100);
    let scaffold = Scaffold {
        header: b"s".to_vec(),
        sequence: joined(&[&left, b"NNN", &right]),
    };
    let scaffolds = vec![scaffold];
    let gaps = discover(&scaffolds, 31);
    let mut read_left = left.clone();
    read_left.insert(88, b'A');
    let mut read_right = right.clone();
    read_right.remove(12);
    let read = joined(&[&read_left, b"GATTACA", &read_right]);
    let index =
        Index::with_windows(&scaffolds, &gaps, 31, 15, 100, 15).with_flank_edit_limit(Some(1));
    for sequence in [&read, &revcomp(&read)] {
        let spans = index.extract(sequence, &mut Scanner::default(), 1000);
        assert_eq!(spans.len(), 1);
        assert!(
            matches!(&spans[0].1,Span::Sequence(s) if s.as_slice()==b"GATTACA"),
            "{spans:?}"
        );
    }
    let exact =
        Index::with_windows(&scaffolds, &gaps, 31, 15, 100, 15).with_flank_edit_limit(Some(0));
    assert!(
        exact
            .extract(&read, &mut Scanner::default(), 1000)
            .is_empty()
    );
}

#[test]
fn report_counts_seed_candidates_rejected_by_flank_context_or_order() {
    let mut f = Fixture::new();
    let mut seed = 817264819372;
    f.left = dna(&mut seed, 100);
    f.right = dna(&mut seed, 100);
    f.write_scaffolds(&[f.span(b"NNN")]);
    let mut left = f.left.clone();
    left[20] = if left[20] == b'A' { b'C' } else { b'A' };
    f.write_reads(&[
        ("context".into(), joined(&[&left, b"ACTG", &f.right])),
        ("order".into(), joined(&[&f.right, b"ACTG", &f.left])),
        ("one-side".into(), f.left.clone()),
    ]);
    run(&f.args(
        "candidate-audit",
        &["--max-flank-edits", "0", "--recruitment", "strict"],
    ))
    .unwrap();
    let report = f.report("candidate-audit");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    assert_eq!(&row[6..11], &["2", "0", "0", "0", "0"]);
    assert_eq!(f.output("candidate-audit")[0].sequence, f.span(b"NNN"));
}

fn long_flank_fixture() -> Fixture {
    let mut f = Fixture::new();
    let mut seed = 781264817392;
    f.left = dna(&mut seed, 100);
    f.right = dna(&mut seed, 100);
    f.write_scaffolds(&[f.span(b"nnn")]);
    f
}
fn rescue_read(f: &Fixture, middle: &[u8]) -> Vec<u8> {
    let mut left = f.left.clone();
    let mut right = f.right.clone();
    for pos in [7, 27] {
        left[pos] = if left[pos] == b'A' { b'C' } else { b'A' };
    }
    for pos in [72, 92] {
        right[pos] = if right[pos] == b'A' { b'C' } else { b'A' };
    }
    joined(&[&left, middle, &right])
}

#[test]
fn adaptive_rescue_retains_stronger_reads_and_is_deterministic() {
    let f = long_flank_fixture();
    for (label, weak, expected) in [
        ("agree", b"AAAA".as_slice(), b"AAAA".as_slice()),
        ("disagree", b"CCCC".as_slice(), b"nnn".as_slice()),
    ] {
        f.write_reads(&[
            ("primary".into(), f.span(b"AAAA")),
            ("weak1".into(), rescue_read(&f, weak)),
            ("weak2".into(), revcomp(&rescue_read(&f, weak))),
        ]);
        let mut previous = None;
        for threads in ["1", "4"] {
            let name = format!("adaptive-{label}-{threads}");
            run(&f.args(&name, &["--threads", threads, "--min-quality", "50"])).unwrap();
            assert_eq!(f.output(&name)[0].sequence, f.span(expected));
            let report = f.report(&name);
            let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
            assert_eq!(row[10], "3");
            assert_eq!(&row[18..22], &["rescue", "1", "2", "1"]);
            if label == "disagree" {
                assert_eq!(row[16], "low_msa_quality");
                assert_eq!(row[14], "0.00000000");
            }
            if let Some(previous) = &previous {
                assert_eq!(&report, previous);
            }
            previous = Some(report);
        }
    }
}

#[test]
fn adaptive_never_uses_rescue_to_bypass_rejected_primary_consensus() {
    let f = long_flank_fixture();
    let mut reads = vec![
        ("p0".into(), f.span(b"AAAA")),
        ("p1".into(), f.span(b"AAAA")),
        ("p2".into(), f.span(b"TTTT")),
    ];
    reads.extend((0..10).map(|i| (format!("weak{i}"), rescue_read(&f, b"AAAA"))));
    f.write_reads(&reads);
    run(&f.args("no-quality-fallback", &["--min-quality", "50"])).unwrap();
    let report = f.report("no-quality-fallback");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    assert_eq!(row[10], "3");
    assert_eq!(row[16], "low_msa_quality");
    assert_eq!(&row[18..22], &["primary", "3", "10", "0"]);
    assert_eq!(f.output("no-quality-fallback")[0].sequence, f.span(b"nnn"));
}

#[test]
fn primary_ambiguity_blocks_an_apparently_unique_weaker_placement() {
    let f = long_flank_fixture();
    let mut left = f.left.clone();
    left[75] = if left[75] == b'A' { b'C' } else { b'A' };
    let mut right = f.right.clone();
    right[25] = if right[25] == b'A' { b'C' } else { b'A' };
    let read = joined(&[&left, b"CC", &right, &f.left, b"GATTACA", &f.right]);
    let scaffolds = input::scaffolds(&f.dir.path().join("scaffolds.fa")).unwrap();
    let gaps = discover(&scaffolds, 31);
    let weak = Index::with_windows(&scaffolds, &gaps, 31, 15, 50, 15)
        .with_recruitment(yagaf::index::Recruitment::BestContext, 2);
    assert!(
        matches!(&weak.extract(&read,&mut Scanner::default(),1000)[0].1,Span::Sequence(s) if s.as_slice()==b"GATTACA")
    );
    f.write_reads(
        &(0..3)
            .map(|i| (format!("r{i}"), read.clone()))
            .collect::<Vec<_>>(),
    );
    let mut args = f.args("blocked-ambiguity", &[]);
    args.k = 31;
    run(&args).unwrap();
    assert_eq!(f.output("blocked-ambiguity")[0].sequence, f.span(b"nnn"));
    let report = f.report("blocked-ambiguity");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    assert_eq!(row[10], "0");
    assert_eq!(row[8], "3");
}

#[test]
fn adaptive_cross_tier_ids_are_counted_once_and_conflicts_reject() {
    let f = long_flank_fixture();
    for (name, middle, reason) in [
        ("duplicate", b"AAAA".as_slice(), "insufficient_support"),
        (
            "conflict",
            b"CCCC".as_slice(),
            "conflicting_duplicate_read_id",
        ),
    ] {
        f.write_reads(&[
            ("same".into(), f.span(b"AAAA")),
            ("same".into(), rescue_read(&f, middle)),
            ("other".into(), rescue_read(&f, b"AAAA")),
        ]);
        run(&f.args(name, &[])).unwrap();
        let report = f.report(name);
        let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
        assert_eq!(row[16], reason);
        assert_eq!(row[9], "1");
        assert_eq!(f.output(name)[0].sequence, f.span(b"nnn"));
        if name == "duplicate" {
            assert_eq!(row[10], "2");
        }
    }
}

#[test]
fn adaptive_checks_cross_tier_id_conflicts_when_primary_is_sufficient() {
    let f = long_flank_fixture();
    f.write_reads(&[
        ("same".into(), f.span(b"AAAA")),
        ("p1".into(), f.span(b"AAAA")),
        ("p2".into(), f.span(b"AAAA")),
        ("same".into(), rescue_read(&f, b"CCCC")),
    ]);
    run(&f.args("primary-id-conflict", &[])).unwrap();
    let report = f.report("primary-id-conflict");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    assert_eq!(row[16], "conflicting_duplicate_read_id");
    assert_eq!(row[9], "1");
    assert_eq!(&row[18..22], &["primary", "3", "1", "0"]);
    assert_eq!(f.output("primary-id-conflict")[0].sequence, f.span(b"nnn"));
}

#[test]
fn adaptive_combined_evidence_limits_reject_whole_gaps() {
    let f = long_flank_fixture();
    f.write_reads(&[
        ("p0".into(), f.span(b"AAAA")),
        ("p1".into(), f.span(b"AAAA")),
        ("w0".into(), rescue_read(&f, b"AAAA")),
        ("w1".into(), rescue_read(&f, b"AAAA")),
    ]);
    for (name, options, reason) in [
        (
            "union-reads",
            vec!["--max-reads-per-gap", "3"],
            "read_support_limit",
        ),
        (
            "union-bytes",
            vec!["--max-gap-bytes", "18"],
            "gap_evidence_limit",
        ),
    ] {
        run(&f.args(name, &options)).unwrap();
        assert!(f.report(name).contains(reason));
        assert_eq!(f.output(name)[0].sequence, f.span(b"nnn"));
    }
    let mut reads = (0..4)
        .map(|i| (format!("p{i}"), f.span(b"AAAA")))
        .collect::<Vec<_>>();
    reads.extend((0..3).map(|i| (format!("w{i}"), rescue_read(&f, b"AAAA"))));
    f.write_reads(&reads);
    run(&f.args("no-resource-fallback", &["--max-reads-per-gap", "3"])).unwrap();
    let report = f.report("no-resource-fallback");
    let row: Vec<_> = report.lines().nth(2).unwrap().split('\t').collect();
    assert_eq!(row[16], "read_support_limit");
    assert_eq!(row[18], "primary");
    assert_eq!(f.output("no-resource-fallback")[0].sequence, f.span(b"nnn"));
}
