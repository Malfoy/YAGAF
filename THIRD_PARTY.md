# Third-party sources

Gapfil is GPL-3.0-only; the full license is in `LICENSE`.

## K2Rmini

- Source: <https://github.com/Malfoy/K2Rmini>
- Revision: `29d39e9a4c81d099bb89c4f39e705322affad508`
- License: GPL-3.0, retained at `vendor/K2Rmini/LICENSE`.
- A clone is retained under `vendor/K2Rmini` as the upstream reference. It is not
  built or invoked by Gapfil.
- `src/index.rs` adapts `src/bin/k2rminimulti.rs`: shared minimizer/k-mer indexes,
  SIMD minimizer prefilter, NtHash SIMD k-mer hashing, and read batching. Gapfil
  replaces fixed query masks with postings of gap/side/strand identifiers, verifies
  actual bases and positions, and extracts one interval per uniquely paired gap.
  Each ACTG run is processed independently; ambiguous bases are never removed to
  concatenate unrelated sequence runs. Reverse complements are indexed explicitly.

## POASTA

- Source: <https://github.com/broadinstitute/poasta>
- Revision: `cd53c50e10cd645138b766f82be53e3f9e96296f`
- License: BSD-3-Clause, retained at `vendor/poasta/LICENSE`.
- `vendor/poasta` is a library-only snapshot of the cloned source. The alignment
  algorithm is retained. CLI binaries and graph import/export functions unrelated
  to alignment were excluded from the compiled library. `Cargo.toml.upstream` and
  `README.upstream.md` retain upstream package information.
- The FASTA MSA exporter writes bytes directly instead of using noodles. Its output
  is parsed exclusively by helicase in Gapfil. No upstream FASTA reader is compiled.
- Fixed an upstream zero-based column padding bug: leading insertions/deletions
  could produce unequal MSA row lengths. Export now tracks the next output column
  explicitly and pads all rows through the last column.
- The small subset of `nonmax::NonMaxU32` used by POASTA is implemented with
  `std::num::NonZeroU32` in `src/nonmax.rs`, preserving numeric ordering and the
  four-byte optional representation. This removes an otherwise uncached dependency
  while keeping this snapshot buildable in the restricted environment. Its bounds,
  ordering, conversions, and representation are covered by the vendored unit test.
- Rust formatting was applied to the snapshot. Uncompiled upstream helper modules
  are retained for reference; Gapfil does not expose their original file readers.

## Helicase and SIMD crates

Every FASTA/FASTQ input, including in-memory MSA FASTA, is parsed by helicase 0.2.0:
<https://github.com/imartayan/helicase>. Gapfil retains Ns, case, headers, and sequence
coordinates. SIMD indexing uses simd-minimizers 2.3.1 and its matching packed-seq and
seq-hash versions. `Cargo.lock` pins the complete dependency graph; the vendored
POASTA snapshot does not require a network Git checkout during a build.
