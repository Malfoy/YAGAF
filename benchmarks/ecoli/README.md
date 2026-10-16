# E. coli reference recovery experiment

Run on 2026-10-16, using the initial Gapfil Rust implementation. Every FASTA read,
including simulation input and corrected-output evaluation, uses helicase.

## Inputs and simulation

Two E. coli K-12 laboratory-strain reference chromosomes were downloaded from NCBI:

- MG1655: [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3),
  4,641,652 bases, assembly GCF_000005845.2.
- W3110: [NC_007779.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_007779.1),
  4,646,332 bases, assembly GCF_000010245.2.

For each chromosome, `examples/recovery_bench.rs` generates exactly 50x source-base
coverage with random circular start positions and strands. Template lengths follow
a normal distribution with mean 15 kb and SD 3 kb, truncated to 1–30 kb (the final
read is shortened to hit the coverage budget exactly). Independent errors occur
with probability 0.001 per template base, split equally among substitutions,
single-base insertions, and single-base deletions. A substitution always changes
the nucleotide. This is an approximate **Q30 HiFi-like simulation**, not an
empirically calibrated PacBio HiFi/CCS error model: it omits sequence-context effects,
homopolymer biases, read-to-read accuracy variation, and systematic errors.

Masks replace either 1,000 or 10,000 unique, uniformly sampled original bases with
`N`. Sampling does not exclude repeat regions, chromosome ends, or nearby masks;
adjacent masks become one N run. Read simulation uses the original, unmasked
chromosome and an independent PRNG stream from mask placement. Seed 42 is used
throughout. The sparse and dense masking cases use identical reads per strain.

Gapfil settings: `--preset hifi --k 31 --min-reads 5 --min-quality 100 --threads 4`.
All computational limits use the documented defaults. Both flanking anchors must
be unique across the input scaffold, including reverse complements. A single
empty spanning interval rejects the whole gap. No read subsampling was used.

## Recovery

| Strain | Masked bases | Correctly recovered | Recovery | Wrong fills (bases) | Kept as N |
| --- | ---: | ---: | ---: | ---: | ---: |
| MG1655 | 1,000 | 895 | 89.50% | 0 | 105 |
| W3110 | 1,000 | 873 | 87.30% | 0 | 127 |
| MG1655 | 10,000 | 7,851 | 78.51% | 0 | 2,149 |
| W3110 | 10,000 | 7,796 | 77.96% | 0 | 2,204 |

Observed fill precision was 100% in all four runs. This is an observation on these
simulations, not a guarantee about real reads. Evaluation checked every filled
interval against the original reference, accounted for accepted replacements of
different lengths, and verified that every unmasked base, record header, and
rejected gap was preserved. All masked bases were accounted for exactly once.

The loss of sensitivity is consistent with the conservative rules. For 10,000
masked MG1655 bases, 1,230 gap runs lacked two complete unambiguous anchors because
of nearby masks, 232 had nonunique anchors, 465 failed unanimity, and 221 had an
empty spanning read interval. W3110 had 1,220, 283, 475, and 220 such runs,
respectively, plus one gap with insufficient support. Counts here refer to gap
runs; the recovery table refers to bases, so the totals can differ slightly.

With one-base gaps, a mismatch makes the sole MSA column nonunanimous: the score
is then zero, even when a strong majority supports the correct base. Lowering a
high threshold from 100% to 95% would not rescue that case. Nearby masks also
invalidate the one immediately adjacent anchor strategy; this explains the lower
recovery at the higher masking density without attributing those rejections to
incorrect consensus calls.

## Runtime and files

These are single-run process measurements on a shared Intel Core i9-13950HX
machine, Rust 1.97.0, release profile with native CPU features. Simulation and
evaluation are excluded from Gapfil runtime. The MG1655 dense run was much slower
than the other runs; these measurements should not be used to compare strains or
claim stable throughput without repetition and control of background load.

| Strain | Masks | Read bases | Gapfil wall time | Peak RSS |
| --- | ---: | ---: | ---: | ---: |
| MG1655 | 1,000 | 232,082,106 | 1.417 s | 30.8 MiB |
| W3110 | 1,000 | 232,316,147 | 1.034 s | 30.8 MiB |
| MG1655 | 10,000 | 232,082,106 | 6.383 s | 82.8 MiB |
| W3110 | 10,000 | 232,316,147 | 0.883 s | 84.2 MiB |

Small reproducible artifacts are retained in this directory:

- `results.tsv`: base/gap recovery, timing, throughput, and RSS.
- `parameters.tsv`: settings, realized read lengths and errors, reference SHA-256.
- `rejection_reasons.tsv`: complete per-case gap-run rejection counts.

Downloaded references and full run artifacts are under ignored `references/` and
`runs/`. Each successful run includes reads, their origins and error counts,
masked/corrected FASTA, mask positions and truth bases, per-gap reports, parameters,
stderr log, and the raw summary. The first attempted dense MG1655 run is preserved
under `runs/mg1655-50x-mask10000-seed42-initial-export-failure`.

That initial run exposed an upstream POASTA FASTA-export padding bug for leading
indels: MSA rows could have unequal lengths. Gapfil rejected the malformed alignment
before writing outputs. The vendored exporter now tracks zero-based output columns
correctly, and a regression test covers leading/trailing indels and single-base
substitutions. Successful results above were obtained after this fix.

## Reproduce

```sh
mkdir -p benchmarks/ecoli/references benchmarks/ecoli/runs
curl -L --fail \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz \
  -o benchmarks/ecoli/references/NC_000913.3.fna.gz
curl -L --fail \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/245/GCF_000010245.2_ASM1024v1/GCF_000010245.2_ASM1024v1_genomic.fna.gz \
  -o benchmarks/ecoli/references/NC_007779.1.fna.gz
cargo build --release --locked --bins --examples
target/release/examples/recovery_bench \
  --reference benchmarks/ecoli/references/NC_000913.3.fna.gz \
  --directory benchmarks/ecoli/runs/mg1655-repeat \
  --coverage 50 --mask-count 10000 --seed 42 --threads 4
```

Use a new run directory, switch the reference for W3110, and set `--mask-count 1000`
for the sparse control. `/usr/bin/time` is used to record peak process RSS on Linux.
