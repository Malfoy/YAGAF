# Recruitment and alignment comparisons

Build the release binary and examples with native CPU features:

```sh
cargo build --release --locked --bins --examples
```

`run.py` compares adjacent anchors, spaced and dense flank windows, seed lengths,
known-flank edit limits, and unique best-context placement. Its MSA matrix compares
POASTA, SPOARS SIMD/scalar, input orders, affine scoring, and computation budgets.
All runs use 32 workers, minimum support 3, and 50% unanimous columns by default.
Rejected reads are never selected by their proposed gap sequence or MSA agreement.

```sh
python3 benchmarks/recruitment/run.py \
  --scaffolds scaffolds.fa --reads reads.fa.gz --output-dir results/recruitment

python3 benchmarks/recruitment/run.py \
  --scaffolds scaffolds.fa --reads reads.fa.gz --output-dir results/msa \
  --matrix msa --only poasta-1000M spoars-1000M spoars-scalar-1000M \
  -- --max-flank-edits 1
```

Each configuration saves its exact command, binary SHA-256, per-phase timing,
process peak RSS, corrected FASTA, per-gap report, rejection counts, and validation.
Runs are sequential to avoid multiplying the thread budget. A timeout kills the
complete process group and records failure. Existing output directories are refused.
Single-run timings on a shared machine are exploratory measurements, not stable
speedup estimates. Parsing and extraction overlap; their times must not be summed.

Validation uses Rust/helicase tools. `verify_scaffolds` verifies all headers,
record order, unreplaced gaps, and every base outside replacements. With
`--reference`, `verify_recovery` additionally checks every filled interval against
known truth, including length-changing errors. A length-changing wrong fill counts
as one wrong gap, and all original masked bases of that interval are counted as
failed recovery; that base count is not an edit distance. No Python FASTA parser
is used.

Public reference experiments reuse the two E. coli chromosomes and 50x independent
error simulations documented in [the reference experiment](../ecoli/README.md).
`prepare_masks` supplies an additional MG1655 workload with 200 non-overlapping gaps
of 50, 100, 250, 500, and 1,000 bases (40 each, seed 42; 76,000 masked bases).
Reference repeats are retained; intervals are rejected during preparation only to
avoid overlapping or directly touching another mask. `recovery_bench --prepare-only`
can generate a second read workload at an explicitly chosen independent error rate.
These are controlled simulations, not calibrated empirical HiFi or ONT models.

Keep full inputs and run artifacts under ignored directories.

## Measured reference results

Measured on 2026-10-17 with 32 workers, an Intel Core i9-13950HX, Rust 1.97.0,
and the native release profile. The six final-default cases in
[truth-results.tsv](truth-results.tsv) had **zero incorrect fills**:

| Dataset | Exact filled gaps | Gap runs | Correct masked bases | Masked bases |
| --- | ---: | ---: | ---: | ---: |
| MG1655, sparse Q30-like | 949 | 999 | 950 | 1,000 |
| MG1655, dense Q30-like | 8,334 | 9,987 | 8,347 | 10,000 |
| W3110, sparse Q30-like | 926 | 1,000 | 926 | 1,000 |
| W3110, dense Q30-like | 8,351 | 9,986 | 8,361 | 10,000 |
| MG1655, 200 contiguous gaps, error 0.001 | 199 | 200 | 75,900 | 76,000 |
| MG1655, 200 contiguous gaps, error 0.01 | 195 | 200 | 75,100 | 76,000 |

The adjacent-anchor MG1655 dense control made one incorrect repeat fill. All tested
window strategies avoided it. Final defaults use 100-base primary context with a
one-edit cap, 50-base rescue, automatic MSA, minimum support 3, unanimity 50%,
maximum support 4,096, and an MSA work budget of 268,435,456. Evidence payload caps
remain 4 MiB per tier/final gap and 512 MiB globally across both tiers.

The MSA controls use the **same** 1%-error reads and permissive 100-base recruitment
(`--recruitment strict --max-flank-edits 22`), so they compare alignment rather than
read selection:

| Engine/order/scoring | Exact fills | Incorrect fills |
| --- | ---: | ---: |
| POASTA, lexical, HiFi | 180 | 4 |
| POASTA, longest-first, HiFi | 183 | 3 |
| POASTA, median-first, HiFi | 179 | 4 |
| SPOARS SIMD/scalar, lexical, HiFi | 190 | 0 |
| SPOARS SIMD, longest-first, HiFi | 186 | 1 |
| SPOARS SIMD, median-first, HiFi | 190 | 0 |
| POASTA, lexical, ONT | 188 | 0 |
| SPOARS SIMD, lexical, ONT | 191 | 0 |
| Automatic, lexical, HiFi | 190 | 0 |

Automatic selection retains POASTA's efficiency for clean evidence while using
SPOARS for higher span divergence. This finite test set does not establish general
accuracy guarantees. All individual results and process measurements are retained
in the TSV. Source input preparation and validation use helicase throughout.

To generate the long-gap/noisier controls, use new directories:

```sh
target/release/examples/prepare_masks \
  --reference benchmarks/ecoli/references/NC_000913.3.fna.gz \
  --directory results/long-masks --gaps 200 --seed 42
target/release/examples/recovery_bench \
  --reference benchmarks/ecoli/references/NC_000913.3.fna.gz \
  --directory results/error01 --coverage 50 --mask-count 1 \
  --error-rate 0.01 --seed 42 --prepare-only
python3 benchmarks/recruitment/run.py \
  --scaffolds results/long-masks/masked.fa --reads results/error01/reads.fa \
  --reference benchmarks/ecoli/references/NC_000913.3.fna.gz \
  --output-dir results/msa-error01 --matrix msa \
  --only poasta-200M spoars-200M auto-200M \
  -- --recruitment strict --flank-window 100 --max-flank-edits 22
```
