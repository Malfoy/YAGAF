# Initial performance measurements

This report compares the first optimization (`a4d0d9d`) against `5c92e89`.
See [parallel parsing](parallel-parsing/README.md) for the next optimization and
shared CLI defaults.

Measured 2026-10-16 on Intel Core i9-13950HX (24 cores, 32 logical CPUs), Linux,
native release builds. Runs execute serially with 32 read and 32 MSA workers.
Inputs are warm in the filesystem cache. Other system activity and CPU frequency
vary; these subsecond timings are workload-specific, not throughput guarantees.

## End-to-end comparison

Baseline: commit `5c92e89`, thin LTO. Optimized: native CPU, opt-level 3, fat LTO,
one codegen unit, aborting panics, stripped symbols. Both use POASTA, HiFi scoring,
`k=31`, minimum 5 reads, 100% unanimous-column acceptance. Baseline batches are
4 MiB; optimized batches are 16 MiB. Each entry is the median of seven runs,
alternating baseline/current order. Process time includes startup and shutdown.
RSS is sampled by `/usr/bin/time` and is the median of per-process maxima.

The datasets and independent Q30-like error model are documented in
[the recovery experiment](../ecoli/README.md). Each input contains 50x coverage,
about 232 Mb of reads. This comparison reuses identical inputs. Every corrected
FASTA byte and every gap-report data row matches across all runs and versions;
parameter comments differ. Previously verified recovery and zero erroneous fills
are preserved. All sequence parsing in generation/evaluation uses helicase.

| Strain | Masked bases | Baseline (s) | Optimized (s) | Speedup | Baseline / optimized RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| mg1655 | 1000 | 0.1922 | 0.1174 | 1.64× | 33.3 / 80.7 |
| mg1655 | 10000 | 0.4390 | 0.1678 | 2.62× | 87.1 / 156.6 |
| w3110 | 1000 | 0.1945 | 0.1194 | 1.63× | 33.2 / 79.3 |
| w3110 | 10000 | 0.4740 | 0.1796 | 2.64× | 87.0 / 154.3 |

## Phase times

MG1655, 10,000 masked bases; milliseconds, independently computed seven-run medians.
Phase medians need not sum to the median total. Internal total excludes process
startup/shutdown and timing-file publication. Parser, extraction, and merge overlap
and must not be added to the sequential phases.

| Phase | Baseline (ms) | Optimized (ms) |
| --- | ---: | ---: |
| setup | 1.194 | 1.441 |
| scaffold_read | 9.828 | 9.049 |
| gap_discovery | 3.719 | 3.626 |
| index_build | 14.993 | 8.092 |
| anchor_uniqueness | 23.157 | 2.224 |
| index_filter | 14.389 | 1.009 |
| read_scan | 231.499 | 110.290 |
| consensus | 127.156 | 14.682 |
| output | 5.257 | 5.628 |
| total | 429.834 | 153.759 |
| read_parse_active | 105.329 | 98.949 |
| read_extraction | 139.641 | 44.454 |
| read_merge | 86.696 | 37.577 |

The remaining read-scan bottleneck is largely serial parsing and read delivery.
Optimizations include minimizer-offset postings with exact verification, parallel
index preparation and scaffold uniqueness scans, pruning instead of rebuilding
the index, reused worker buffers, shared read IDs with collision-checked hash
lookup, parallel evidence merging, and direct in-memory POASTA MSA rows.
Parallel merging falls back to input-order merging near the global evidence cap.

A randomized seven-repeat 1/2/4/8/16/32 MiB batch sweep selected 16 MiB locally.
The larger bounded queue trades additional memory for throughput; use
`--batch-bases 4194304` for the previous lower-memory batching. Raw sweep timings
are separate from the alternating comparison and should not be compared across
sessions as if system conditions were identical.

## MSA engines

Native Rust POASTA, SPOARS SIMD, and SPOARS scalar receive the same sorted,
reference-derived families: 32 independent gaps, 50 reads at lengths 128/500,
16 reads at length 1000, seed 42, 32 workers, five repetitions. Generation and
cloning are excluded from timed alignment/consensus. All families contain at least
one disagreement, so the identical-span fast path cannot dominate. Each input
base has an independent error probability split equally among substitution,
insertion, and deletion. HiFi scoring uses 0.1% errors; ONT scoring uses 5%.
These are controlled stress tests, not empirical platform error models.

| Preset | Length | POASTA (ms) | SPOARS SIMD (ms) | SPOARS scalar (ms) | Exact consensus families POASTA / SIMD |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hifi | 128 | 5.24 | 3.41 | 8.51 | 32 / 32 |
| Hifi | 500 | 20.49 | 41.32 | 136.54 | 32 / 32 |
| Hifi | 1000 | 11.54 | 184.69 | 231.51 | 32 / 32 |
| Ont | 128 | 74.38 | 6.33 | 14.48 | 14 / 32 |
| Ont | 500 | 1737.25 | 202.65 | 330.50 | 3 / 29 |
| Ont | 1000 | 1249.13 | 258.08 | 292.75 | 15 / 29 |

SPOARS SIMD/scalar consensus counts and quality numerators/denominators match
in every repetition. POASTA and SPOARS can place indels differently and produce
different MSA quality and consensus results. POASTA remains the default; SPOARS
is optional via `--msa-backend spoars`. Exact counts above describe consensus
sequence equality to truth before applying a minimum-quality acceptance threshold.
See raw `msa.tsv` for rejected consensuses, column counts, and peak RSS.

## Reproduction

Prepare the four recovery datasets using `examples/recovery_bench.rs` and the
accessions/settings in `../ecoli`. Build baseline commit `5c92e89` separately and
the current tree with `cargo build --release --locked --bins --examples`.

```sh
python3 benchmarks/performance/run.py \
  --baseline /path/to/baseline/gapfil --binary target/release/yagaf \
  --data benchmarks/ecoli/runs --output /tmp/yagaf-comparison --repeats 7

python3 benchmarks/performance/msa.py \
  --binary target/release/examples/msa_bench \
  --reference benchmarks/ecoli/references/NC_000913.3.fna.gz \
  --output /tmp/yagaf-msa --repeats 5
```

Raw measurements: `runs.tsv`, `phases.tsv`, `msa.tsv`, `batch-sweep.tsv`.
No hardware performance counters were available (`perf_event_paranoid=4`).
