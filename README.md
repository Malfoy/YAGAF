# YAGAF

Conservative, parallel long-read scaffold gap filling in Rust.

```sh
cargo build --release --locked

target/release/yagaf \
  --scaffolds scaffolds.fa --reads reads.fq.gz \
  --output corrected.fa --report gaps.tsv \
  --preset hifi --k 31 --min-reads 5 --min-quality 100 --threads 32
```

Use `--preset hifi` or `ont`. Choose `k`, minimum read support, and minimum
percentage of unanimous MSA columns explicitly. Unsupported gaps stay unchanged.

All sequence input uses **helicase**, with automatic gzip, zstd, bzip2, and xz
detection. **POASTA** is the default MSA engine; `--msa-backend spoars` selects
native Rust SIMD POA. `--timings phases.tsv` saves phase runtimes.

Build on the target machine for native CPU optimization. Output files must be new.
See `yagaf --help`, [usage](docs/usage.md), and [licenses](THIRD_PARTY.md).
