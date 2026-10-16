# Gapfil

Conservative, parallel scaffold gap filling in Rust, using long reads. This is an
initial implementation with synthetic correctness tests; HiFi/ONT acceptance
settings and real-data accuracy have not yet been calibrated.

Gapfil finds each `N`/`n` run, takes one immediately adjacent k-mer per side, and
adapts K2Rmini's multi-query minimizer/k-mer indexing to find spanning reads in one
streamed pass. Exact bases, strand, order, and unique pairing are verified after
hash matches. Both anchors must also be unique across the scaffold input, considering
both strands. Unsupported or ambiguous gaps remain unchanged.

Spans are oriented to the scaffold and aligned with the vendored POASTA Rust
library. The consensus is the unique majority symbol in each MSA column, including
deletions; any tied column rejects the fill. All FASTA/FASTQ reading uses **helicase**,
including parsing POASTA's in-memory MSA export.

## Build

```sh
cargo build --release --locked
```

The binary is `target/release/gapfil`. `.cargo/config.toml` enables
`target-cpu=native` for helicase and SIMD indexing. Build on the machine where you
will run the binary; a native binary may not run on another CPU. AVX2 on x86-64 or
NEON on aarch64 is required by this dependency combination. For older AMD CPUs with
slow PDEP, build with `--features no-pdep`. This repository does not currently
provide a baseline portable binary.

Pinned upstream sources, adaptations, and licenses are described in
[THIRD_PARTY.md](THIRD_PARTY.md).

## Run

Choose `--preset hifi` or `--preset ont`, and explicitly provide `--k`,
`--min-reads`, and `--min-quality`. There are no implicit acceptance thresholds.
For example, this illustrates a strict configuration; these values are **not
validated recommendations**:

```sh
target/release/gapfil \
  --scaffolds scaffolds.fa --reads reads.fq.gz \
  --output corrected.fa --report gaps.tsv \
  --preset hifi --k 31 --min-reads 5 --min-quality 100 \
  --threads 8
```

FASTA scaffolds and FASTA/FASTQ reads may be uncompressed, gzip (`.gz`), zstd
(`.zst`), bzip2 (`.bz`/`.bz2`), or xz (`.xz`). Compression is detected from file
contents, independently of the filename extension; no format flags are needed.
Paths must refer to files; stdin/stdout streaming is not yet exposed. Both output
paths must be new files in existing directories. Input files are never overwritten.
Headers, record order, and all bases outside replaced gaps are preserved, including
case. FASTA sequence wrapping is normalized to one line per record.

Add `--timings phases.tsv` to save phase timings. The TSV separates sequential
wall-time phases from overlapping parser, extraction, and evidence-merge timings;
overlapping values should not be added to the total. Phase times are also logged
to stderr. The timing output must be a new file separate from the other outputs.

The two presets initially change global affine alignment penalties:

| Preset | Mismatch | Gap open | Gap extend |
| --- | ---: | ---: | ---: |
| HiFi | 4 | 6 | 2 |
| ONT | 4 | 4 | 1 |

These scoring choices are provisional. ONT still requires exact matches to both
flank k-mers; errors in those anchors reduce sensitivity. The Ns' original count
is never treated as an accurate expected biological length.

## Acceptance and quality

Both minimum support and minimum quality must pass:

`msa_quality_percent = 100 × unanimous_columns / evaluated_columns`

A column is unanimous only if every retained spanning sequence has the same
unambiguous A/C/G/T. Insertions, deletions, and missing coverage prevent unanimity.
All-gap columns are excluded. Only the gap sequence is evaluated, excluding both
anchors. A score of 100 therefore means complete agreement in the extracted spans.
Identical spans use the equivalent trivial MSA fast path.

Read IDs are the first whitespace-delimited header field. Repeated IDs with the
same span count once; conflicting spans for the same ID reject the gap. Multiple
anchor occurrences or multiple valid strand pairings in one read do not support a
gap. A uniquely paired span containing ambiguous bases or no bases rejects the
gap, even when other reads agree; it is not silently discarded to raise quality.

The TSV report contains original zero-based, half-open coordinates, scaffold record
index and ID, candidate/spanning record counts, ambiguous pairing and duplicate-ID
counts, retained independent support, consensus length, unanimous/total MSA columns,
quality percentage, status, and reason. Missing metrics are `NA`. The leading `#`
line records parameters; ID bytes requiring escaping use `%XX`. `spanning_records`
counts uniquely paired input records before duplicate-ID consolidation. Evidence
is released on resource rejection, so `retained_reads` is then zero. Initial flank
rejections are excluded from the read index and have zero candidate counts.

## Resource controls

Read scanning uses bounded batches, a two-batch producer queue, and a shared
immutable index. Extraction runs in parallel; results merge in input order.
Independent MSAs run with at most `min(threads,2)` workers by default; use
`--msa-threads` to change that within the read-thread budget. Alignment inputs are
sorted so results do not depend on worker scheduling or batch boundaries.

| Flag | Default | Behavior |
| --- | ---: | --- |
| `--max-span` | 50,000 bases | Reject gap if any unique span exceeds limit |
| `--max-reads-per-gap` | 256 | Reject whole gap on overflow; never subsample |
| `--max-gap-bytes` | 4 MiB | Retained sequence + read-ID bytes per gap |
| `--max-evidence-bytes` | 512 MiB | Retained evidence payload across all gaps |
| `--max-msa-cells` | 20,000,000 | Bound `sum(span lengths) × max(span length)` before nontrivial MSA |
| `--batch-bases` | 4 MiB | Target batch bytes, counting sequences and IDs |
| `--max-read-bases` | 10,000,000 | Fail run on oversized read |

These are computational guardrails, not biological length filters. If a gap crosses
a capacity limit, all its retained evidence is released and the original gap is
kept. Global capacity decisions follow input order. They can therefore differ when
the read file is reordered, while remaining reproducible across worker counts.
Limits do not constitute a hard RSS cap: scaffolds, index/postings, record buffers,
temporary extraction results, allocator overhead, and POASTA states also use memory.
Helicase must parse a complete record before the read-size check. Individual
scaffolds are currently limited to `2^31-1` bases by the indexing implementation.

Processing errors occur before either output is published. Each output is published
complete without overwriting an existing file. Publication of two paths is not a
single atomic transaction; if FASTA publication fails after report publication, the
error identifies the report that was already written.

## Verification and synthetic benchmark

```sh
cargo fmt --all --check
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
cargo test --manifest-path vendor/poasta/Cargo.toml --lib
cargo build --release --locked --examples
target/release/examples/synthetic_bench --directory /tmp/gapfil-benchmark
```

The benchmark creates deterministic synthetic scaffolds and reads, runs separate
one- and four-worker processes, verifies every output sequence against synthetic
truth through helicase, and reports wall time, input-base throughput, and peak RSS
on Linux using `/usr/bin/time`. The destination directory must not exist. The
default workload has 100 Mb of reads and 128 supported gaps, with one mismatch in
one of five spans per gap to exercise POASTA. This is a reproducible smoke benchmark,
not a substitute for HiFi/ONT assembly benchmarks or a real-data accuracy estimate.

The requested E. coli recovery experiment, with 50x simulated reads and random
reference bases replaced by N, is documented in
[benchmarks/ecoli/README.md](benchmarks/ecoli/README.md), including results, seeds,
simulation assumptions, reference accessions, and reproduction commands.
