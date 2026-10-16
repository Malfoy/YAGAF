# YAGAF

Conservative, parallel scaffold gap filling in Rust using long reads.

YAGAF finds each `N`/`n` run, takes one immediately adjacent k-mer per side, and
adapts K2Rmini's multi-query minimizer/k-mer indexing to find spanning reads in one
streamed pass. Exact bases, strand, order, and unique pairing are verified after
hash matches. Flanking k-mers may repeat elsewhere in the assembly. Ambiguous
pairings within a read are excluded. Unsupported gaps remain unchanged.

Spans are oriented to the scaffold and aligned with the vendored POASTA Rust
library. The consensus is the unique majority symbol in each MSA column, including
deletions; any tied column rejects the fill. All FASTA/FASTQ reading uses **helicase**,
with alignment rows passed directly in memory.

## Build

```sh
cargo build --release --locked
```

The binary is `target/release/yagaf`. Release builds use optimization level 3, fat
LTO, one codegen unit, aborting panics, and stripped symbols. `.cargo/config.toml` enables
`target-cpu=native` for helicase and SIMD indexing. Build on the machine where you
will run the binary; a native binary may not run on another CPU. AVX2 on x86-64 or
NEON on aarch64 is required by this dependency combination. For older AMD CPUs with
slow PDEP, build with `--features no-pdep`. This repository does not currently
provide a baseline portable binary.

Pinned upstream sources, adaptations, and licenses are described in
[THIRD_PARTY.md](../THIRD_PARTY.md).

## Run

Only scaffold, read, and corrected-output paths are required. Shared defaults are
`k=31`, a minimum of 3 independent reads, and at least 50% unanimous MSA columns.
All available CPU workers are used. These settings apply independently of read
technology. Override any setting when needed:

```sh
target/release/yagaf \
  --scaffolds scaffolds.fa --reads reads.fq.gz \
  --output corrected.fa
```

FASTA scaffolds and FASTA/FASTQ reads may be uncompressed, gzip (`.gz`), zstd
(`.zst`), bzip2 (`.bz`/`.bz2`), or xz (`.xz`). Compression is detected from file
contents, independently of the filename extension; no format flags are needed.
The report defaults to `OUTPUT.gaps.tsv`; use `--report` to set another path.
`-h` shows basic options, while `--help` includes all advanced overrides.
Paths must refer to files; stdin/stdout streaming is not yet exposed. Both output
paths must be new files in existing directories. Input files are never overwritten.
Headers, record order, and all bases outside replaced gaps are preserved, including
case. FASTA sequence wrapping is normalized to one line per record.

POASTA is the default MSA backend. Use `--msa-backend spoars` for native Rust SIMD
POA or `--msa-backend spoars-scalar` for its scalar comparison engine. All backends
use the same unanimity metric and conservative acceptance rules.

Add `--timings phases.tsv` to save phase timings. The TSV separates sequential
wall-time phases from overlapping parser, extraction, and evidence-merge timings;
overlapping values should not be added to the total. In mapped mode,
`read_parse_active` measures producer framing, and `read_extraction` includes worker
parsing and span extraction. Phase times are also logged to stderr. The timing output must be a new file separate from the other outputs.

The optional legacy `--preset hifi|ont` flag changes only global affine alignment
penalties; `hifi` scoring is the default. It does not change k, backend, acceptance
thresholds, or resource settings:

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
immutable index. Uncompressed regular read files use helicase record framing
followed by parallel helicase parsing and span extraction from a read-only memory
map. Complete read sequences stay local to each parser; only extracted spans and
read IDs enter the evidence pipeline. Spans of up to 16 bases use inline storage.
Compressed inputs use streamed helicase parsing. `--stream-reads` forces that path
for uncompressed inputs too, reducing mapped-file RSS. Keep mapped inputs unchanged
while the program runs. Evidence merges across independent gaps in parallel,
preserving each gap's read order. Near the global memory limit, merging follows
input order serially so capacity decisions remain deterministic.
Independent MSAs use the configured thread count by default; use
`--msa-threads` to change that within the read-thread budget. Alignment inputs are
sorted so results do not depend on worker scheduling or batch boundaries.

| Flag | Default | Behavior |
| --- | ---: | --- |
| `--max-span` | 50,000 bases | Reject gap if any unique span exceeds limit |
| `--max-reads-per-gap` | 256 | Reject whole gap on overflow; never subsample |
| `--max-gap-bytes` | 4 MiB | Retained sequence + read-ID bytes per gap |
| `--max-evidence-bytes` | 512 MiB | Retained evidence payload across all gaps |
| `--max-msa-cells` | 20,000,000 | Bound `sum(span lengths) × max(span length)` before nontrivial MSA |
| `--batch-bases` | 16 MiB | Target batch bytes (encoded bytes for mapped input; sequence + ID bytes when streaming) |
| `--max-read-bases` | 10,000,000 | Fail run on oversized read |

These are computational guardrails, not biological length filters. If a gap crosses
a capacity limit, all its retained evidence is released and the original gap is
kept. Global capacity decisions follow input order. They can therefore differ when
the read file is reordered, while remaining reproducible across worker counts.
Limits do not constitute a hard RSS cap: scaffolds, index/postings, record buffers,
temporary extraction results, allocator overhead, and POASTA states also use memory.
Mapped input pages also contribute to RSS, so mapped-mode peak RSS can be larger
than with `--stream-reads`, particularly for warm cached input.
Helicase must parse a complete record before the read-size check. Individual
scaffolds are currently limited to `2^31-1` bases by the indexing implementation.

Processing errors occur before either output is published. Each output is published
complete without overwriting an existing file. Publication of two paths is not a
single atomic transaction; if FASTA publication fails after report publication, the
error identifies the report that was already written.

## Development checks

```sh
cargo fmt --all --check
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
```
