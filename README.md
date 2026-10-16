# YAGAF

Conservative, parallel long-read scaffold gap filling in Rust.

```sh
cargo build --release --locked
target/release/yagaf -s scaffolds.fa -r reads.fq.gz -o corrected.fa
```

Defaults: **k=31**, **3 independent reads**, **50% unanimous MSA columns**, and all
available CPU threads. The report is `corrected.fa.gaps.tsv`. Unsupported gaps stay
unchanged. All settings can be overridden.

All sequence input uses **helicase**, with automatic gzip, zstd, bzip2, and xz
detection. POASTA is the default MSA engine; `--msa-backend spoars` selects native
Rust SIMD POA. `--timings phases.tsv` saves phase runtimes.

Build on the target machine for native CPU optimization. Output files must be new.
Use `-h` for basic options or `--help` for all overrides. See [usage](docs/usage.md)
and [licenses](THIRD_PARTY.md).
