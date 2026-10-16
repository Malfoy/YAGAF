# Parallel parsing and shared CLI defaults

Measured 2026-10-16 on the same i9-13950HX, native release builds, 32 read/MSA
workers, warm cached input. Baseline is commit `a4d0d9d`. Nine repetitions per
dataset alternate baseline/current order; no compilations or other benchmarks
run concurrently. Input generation and simulation caveats are documented in
[the recovery experiment](../../ecoli/README.md).

## Implementation

Uncompressed regular reads are mapped read-only. A lightweight helicase parser
identifies complete record chunks; workers parse those chunks through helicase
and extract spans before advancing their local parser. This eliminates the
producer copy/ownership transfer of complete reads. Batches preserve input order
and use the existing two-batch queue. Compressed input retains streamed helicase
parsing; `--stream-reads` forces that path for uncompressed files too.

Extracted spans up to 16 bases use inline storage. Consensus operates directly
on either inline or heap-backed evidence, retaining every contributing read.
All anchor, duplicate-ID, capacity, and consensus-quality rules are preserved.

An attempted read-buffer recycling queue was removed after seven-repeat tests
showed it was slower. Inline storage reduced the streamed path's evidence
allocation cost and memory use. The default mapped path trades larger RSS
(including file-backed input pages) for throughput. No unsafe code was added.

## Performance with unchanged acceptance thresholds

These runs explicitly retain the old settings: k=31, 5 reads, 100% unanimity,
POASTA and HiFi alignment penalties, 16 MiB batches. The new shared CLI defaults
are evaluated separately below. Every corrected FASTA byte and gap-report data
row matches between versions within each dataset across all 72 comparison
processes.

| Strain | Masked bases | Baseline process (s) | Current process (s) | Speedup | Baseline / current peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| mg1655 | 1000 | 0.1084 | 0.0719 | 1.51× | 70.9 / 237.9 |
| mg1655 | 10000 | 0.1493 | 0.1093 | 1.37× | 135.1 / 286.6 |
| w3110 | 1000 | 0.1124 | 0.0756 | 1.49× | 73.4 / 238.2 |
| w3110 | 10000 | 0.1551 | 0.1121 | 1.38× | 134.5 / 286.6 |

MG1655, 10,000 masked bases: independent phase medians in milliseconds.
Internal total excludes process startup/shutdown and timing-file publication;
independent phase medians need not sum to the median total.

| Phase | Baseline (ms) | Current (ms) |
| --- | ---: | ---: |
| setup | 1.250 | 1.430 |
| scaffold_read | 9.875 | 8.701 |
| gap_discovery | 3.037 | 3.040 |
| index_build | 7.910 | 7.637 |
| anchor_uniqueness | 1.823 | 1.864 |
| index_filter | 0.691 | 0.808 |
| read_scan | 99.018 | 67.897 |
| consensus | 11.823 | 3.803 |
| output | 4.648 | 4.792 |
| total | 138.497 | 100.496 |
| read_parse_active | 89.907 | 32.735 |
| read_extraction | 33.142 | 41.194 |
| read_merge | 25.617 | 19.163 |

`read_parse_active`, `read_extraction`, and `read_merge` overlap: do not add
them to the sequential phases. In mapped mode the first measures producer
record framing; the second includes parallel full parsing and span extraction.

A separate alternating seven-repeat MG1655/10,000-mask comparison using
`--stream-reads` reduced median internal total from 137.65 to 123.02 ms and
median peak RSS from 138.4 to 124.5 MiB. Raw data is in `streaming.tsv`.
System conditions differ between sessions; compare versions within each run set.

## Shared defaults

The CLI now requires only `--scaffolds`, `--reads`, and `--output`. Shared defaults
are k=31, 3 independent reads, 50% unanimous columns, POASTA, and all available
CPU workers. Reports default to `OUTPUT.gaps.tsv`. `-h` lists basic options;
`--help` lists all overrides. Optional legacy alignment-scoring presets do not
change these shared values.

Reusing the exact same four prepared read/mask datasets with the new defaults
gives the results below. `examples/verify_recovery.rs` checks reference equality
inside accepted fills, retention of rejected gaps, and preservation of every
unmasked base and header. All FASTA parsing uses helicase.

| Strain | Masked bases | Correct bases | Wrong bases | Unchanged Ns | Recovery |
| --- | ---: | ---: | ---: | ---: | ---: |
| mg1655 | 1000 | 928 | 0 | 72 | 92.80% |
| mg1655 | 10000 | 8175 | 0 | 1825 | 81.75% |
| w3110 | 1000 | 911 | 0 | 89 | 91.10% |
| w3110 | 10000 | 8146 | 0 | 1854 | 81.46% |

The strict and shared-default runs differ because acceptance thresholds changed,
not because of an optimization changing alignment or evidence. These controlled
tests do not establish real-data accuracy or calibrate the defaults.

## Reproduction

Build baseline commit `a4d0d9d` separately, then build this tree using
`cargo build --release --locked --bins --examples` and prepare the four datasets
as documented in `../../ecoli`.

```sh
python3 benchmarks/performance/run.py \
  --baseline /path/to/a4d0d9d/yagaf --binary target/release/yagaf \
  --data benchmarks/ecoli/runs --output /tmp/yagaf-parsing --repeats 9

target/release/yagaf \
  -s benchmarks/ecoli/runs/mg1655-50x-mask10000-seed42/masked.fa \
  -r benchmarks/ecoli/runs/mg1655-50x-mask10000-seed42/reads.fa \
  -o /tmp/mg1655-defaults.fa

target/release/examples/verify_recovery \
  --reference benchmarks/ecoli/references/NC_000913.3.fna.gz \
  --masked benchmarks/ecoli/runs/mg1655-50x-mask10000-seed42/masked.fa \
  --corrected /tmp/mg1655-defaults.fa --report /tmp/mg1655-defaults.fa.gaps.tsv
```
