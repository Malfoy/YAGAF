# YAGAF

Conservative, parallel scaffold gap filling in Rust using long reads.

YAGAF finds each `N`/`n` run and indexes overlapping exact seeds from up to 100
known bases on each side. It adapts K2Rmini's multi-query minimizer/k-mer indexing
to find spanning reads in one streamed pass. Exact seed bases, strand, order,
flank context, and unique spanning intervals are verified after hash matches. Flanking k-mers may repeat elsewhere in the assembly. Ambiguous
pairings within a read are excluded. Unsupported gaps remain unchanged.

Spans are oriented to the scaffold and aligned through POASTA or SPOARS Rust APIs. The consensus is the unique majority symbol in each MSA column, including
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
`k=31`, adaptive 100-base flanks with 50-base rescue, a minimum of 3 independent reads, and at least
50% unanimous MSA columns.
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

The default `--msa-backend auto` chooses POASTA for clean spans and native Rust
SPOARS SIMD for more divergent spans. Routing compares up to 64 evenly spaced
exact 15-mers per span against the most frequent span (ties prefer median length,
then lexical order). If more than 10% are absent, it selects SPOARS; shorter spans
use shorter words. This is an engine-routing heuristic, not an identity estimate
or acceptance gate. Every recruited read remains in the MSA with its original
weight. Identical spans bypass alignment. Use `--msa-backend poasta`, `spoars`, or
`spoars-scalar` to force an engine. All backends
use the same unanimity metric and conservative acceptance rules. `--msa-order`
compares lexical, longest-first, and median-first progressive alignment orders.
Every order retains every recruited sequence; lexical is the default.

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

These scoring choices are provisional. At least one exact seed per side is required;
errors affecting all indexed seeds still reduce sensitivity. The Ns' original count
is never treated as an accurate expected biological length.

## Flank recruitment

`--flank-window 100` sets the maximum known context per side. Windows stop at
scaffold ends and ambiguous bases, including another gap. Two complete unambiguous
adjacent k-mers must be available. Seed starts are spaced by `--anchor-step`
(default `k/2`, or 15 for k=31); both window ends are always included, giving six
seeds per side for a full 100-base window. Set `--anchor-step 1` to index every
overlapping k-mer, or `--flank-window 31` with k=31 for adjacent seeds only.

A seed away from the gap is extended through the intervening known scaffold
sequence with anchored edit alignment. Projection requires a unique best read
endpoint. The complete available flank is then checked outward from that boundary.
Each extension permits at most `floor(context_length/5)+2` edits.
`--max-flank-edits` caps primary/single-tier extensions further: the default is
**1 edit**; 0 requires exact context in that tier. Adaptive rescue uses the normal
budget unless `--rescue-max-flank-edits` is set. This mapping
tolerance is independent of the MSA quality gate. Flanks up to 128 bases use cached
bit-parallel edit masks; longer requested windows use banded dynamic programming.
Known flank bases are excluded from both consensus and quality scoring. Repeated
seeds are allowed when they establish one validated interval; two distinct valid
intervals reject that read's evidence. Multiple seed hits never increase read support.

The default `--recruitment adaptive` indexes two tiers together and reads the input
file once:

- Primary: up to 100 known bases per flank, at most one edit per extension,
  and one validated spanning interval per read.
- Rescue: up to 50 known bases per flank (`min(flank-window,50)`, at least k),
  the normal `floor(context_length/5)+2` edit budget, and unique best-context
  placement with the configured score margin. `--rescue-flank-window` and
  `--rescue-max-flank-edits` override its defaults.

A gap with enough independent primary reads uses that tier. Otherwise it may use
rescue-only reads plus every retained primary read. A read mapped in the primary
tier never contributes a second vote in rescue; known primary ambiguity also
blocks that read's rescue placement. Duplicate IDs are consolidated across tiers,
and conflicting spans reject the gap. A rejected primary tier never triggers
rescue. Tier selection happens **before any MSA** and uses only mapping criteria
and support count: failed quality, a tied consensus, or an MSA resource limit
cannot trigger a weaker retry. No selection uses agreement with a proposed fill.
The report shows both input support counts, the selected tier, and primary reads
added to the rescue MSA. All contributing sequences enter the reported quality.

`--recruitment strict` disables rescue. The single-tier `best-context` alternative
ranks intervals by the sum of the two known-flank edit costs, charging an additional
edit per unavailable context base. The best interval must beat every alternative by
at least `--recruitment-margin` edits (default 2); ties and competing strands remain
ambiguous. Gap sequence, gap length, consensus, and MSA quality never participate
in recruitment scoring.

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
distinct validated intervals or multiple valid strand pairings in one read do not
support a gap. A uniquely paired span containing ambiguous bases or no bases rejects the
gap, even when other reads agree; it is not silently discarded to raise quality.

The TSV report contains original zero-based, half-open coordinates, scaffold record
index and ID, candidate/spanning record counts, ambiguous pairing and duplicate-ID
counts, retained independent support, consensus length, unanimous/total MSA columns,
quality percentage, status, reason, and the actual MSA backend. `identical` denotes
the trivial exact-span path; `NA` denotes no alignment. The backend column is
appended after the original 17 columns. Schema 4 also appends `recruitment_tier`,
`primary_reads`, `rescue_reads`, and `supplemental_primary_reads`, followed by
shortest/longest span, retained span bases, estimated MSA work, and the quality
upper bound. The bound is distinct from measured MSA quality; lengths make it
auditable even when no alignment was run. Missing metrics are
`NA`. The leading `#`
line records parameters; ID bytes requiring escaping use `%XX`. `spanning_records`
counts uniquely paired input records before duplicate-ID consolidation.
`candidate_records` counts input records with exact seeds from both sides in a
common orientation, including candidates rejected by flank context or ordering.
Candidate/spanning/ambiguity counts describe the **selected tier** before primary
supplementation. Their difference counts candidates not contributing a span in
that tier, including primary records excluded from rescue. `primary_reads` and
`rescue_reads` report independent input support in the two exclusive tiers (`NA`
when rejection released it). `supplemental_primary_reads` records added
primary IDs; consequently final `retained_reads` can exceed the selected tier's
`spanning_records`. Repeated IDs never increase final support. Report comments
identify this scheme as `selected-tier-paired-seeds`. Evidence
is released on resource rejection, so `retained_reads` is then zero. Initial flank
rejections are excluded from the read index and have zero candidate counts.

## Resource controls

Read scanning uses bounded batches, a two-batch producer queue, and a shared
immutable index. A compact side-presence search avoids storing seed-hit lists;
only reads with paired seed candidates need the second in-memory index search for
placement. Proven ambiguity stops further alignment for that read/gap. The input
file is still read once. Uncompressed regular read files use helicase record framing
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
| `--max-reads-per-gap` | 4,096 | Reject whole gap on overflow; never subsample |
| `--max-gap-bytes` | 4 MiB | Sequence + ID bytes per gap/tier; also enforced on the final combined evidence |
| `--max-evidence-bytes` | 512 MiB | Retained evidence payload across all gaps |
| `--max-msa-cells` | 268,435,456 | Bound `sum(span lengths) × max(span length)` before nontrivial MSA |
| `--batch-bases` | 16 MiB | Target batch bytes (encoded bytes for mapped input; sequence + ID bytes when streaming) |
| `--max-read-bases` | 10,000,000 | Fail run on oversized read |

Before running an MSA, YAGAF uses the bound
`quality <= 100 * shortest_span / longest_span`. Every unanimous column needs a base
from every read, and the MSA must contain every base of the longest span. If this
bound is below the requested threshold, `msa_quality_bound` rejects the gap without
alignment or a claimed measured quality. All reads remain counted; none are trimmed
or removed to change the score.

These are computational guardrails, not biological length filters. If a gap crosses
a capacity limit, all its retained evidence is released and the original gap is
kept. Global capacity decisions follow input order. They can therefore differ when
the read file is reordered, while remaining reproducible across worker counts.
Adaptive recruitment temporarily retains two bounded evidence tiers per gap; the
512 MiB total payload limit is shared across both, and unused tiers are released
before alignment. The final combined evidence also obeys the per-gap limits.
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
