# Third-party sources

YAGAF is GPL-3.0-only; the full license is in `LICENSE`.

## K2Rmini

- Source: <https://github.com/Malfoy/K2Rmini>
- Revision: `29d39e9a4c81d099bb89c4f39e705322affad508`
- License: GPL-3.0, retained at `vendor/K2Rmini/LICENSE`.
- A clone is retained under `vendor/K2Rmini` as the upstream reference. It is not
  built or invoked by YAGAF.
- `src/index.rs` adapts `src/bin/k2rminimulti.rs`: shared minimizer/k-mer indexes,
  SIMD minimizer selection, NtHash SIMD k-mer hashing, and read batching. YAGAF
  replaces fixed query masks with postings of gap/side/strand identifiers, verifies
  actual bases and positions, and extracts one interval per uniquely paired gap.
  Minimizer-offset postings locate candidate starts directly, with full k-mer
  hashing retained for short minimizers or large posting lists. Index construction
  runs in parallel; repeated assembly anchors remain indexed for every gap. Each
  ACTG run is processed independently; ambiguous bases are never removed to
  concatenate unrelated sequence runs. Reverse complements are indexed explicitly.
  Overlapping seeds from bounded flank windows share this index. Cached Myers
  bit-vector alignment projects seed hits through known flanks; unique validated
  read intervals contribute at most one span per gap and read. Primary and rescue
  queries share this index; primary ambiguity blocks weaker placement, and tier
  selection uses support before any alignment score is computed.

## POASTA

- Source: <https://github.com/broadinstitute/poasta>
- Revision: `cd53c50e10cd645138b766f82be53e3f9e96296f`
- License: BSD-3-Clause, retained at `vendor/poasta/LICENSE`.
- `vendor/poasta` is a library-only snapshot of the cloned source. The alignment
  algorithm is retained. CLI binaries and graph import/export functions unrelated
  to alignment were excluded from the compiled library. `Cargo.toml.upstream` and
  `README.upstream.md` retain upstream package information.
- Alignment rows pass directly from POASTA to YAGAF in memory, avoiding FASTA
  serialization and reparsing. The retained FASTA exporter wraps the same row
  construction. No upstream FASTA reader is compiled.
- Fixed an upstream zero-based column padding bug: leading insertions/deletions
  could produce unequal MSA row lengths. Export now tracks the next output column
  explicitly and pads all rows through the last column.
- The small subset of `nonmax::NonMaxU32` used by POASTA is implemented with
  `std::num::NonZeroU32` in `src/nonmax.rs`, preserving numeric ordering and the
  four-byte optional representation. This removes an otherwise uncached dependency
  while keeping this snapshot buildable in the restricted environment. Its bounds,
  ordering, conversions, and representation are covered by the vendored unit test.
- Rust formatting was applied to the snapshot. Uncompiled upstream helper modules
  are retained for reference; YAGAF does not expose their original file readers.

## SPOARS

- Source: <https://github.com/fg-labs/spoars>
- Published crate: `spoars` 0.1.4, downloaded from crates.io.
- Archive SHA-256: `501e9faa2a392e46084161ff4ec615ea13bdf20cef884891aa353d451dd7a80e`.
- License: MIT, retained at `vendor/spoars/LICENSE`; upstream provenance is in
  `vendor/spoars/THIRD-PARTY.md`.
- `vendor/spoars` contains its native Rust library, including SIMD and scalar POA
  engines. CLI, tests, oracle, and file parser dependencies are excluded. The
  library aligns extracted spans in memory; it does not read sequence files.
- Global affine scoring accounts for the different first-gap penalty convention.
  MSA rows use the same unanimity and consensus checks as POASTA.

## Helicase and SIMD crates

Every FASTA/FASTQ input is parsed by helicase 0.2.0:
<https://github.com/imartayan/helicase>. YAGAF retains Ns, case, headers, and sequence
coordinates. The uncompressed read path also uses helicase to frame record
chunks before parsing them in parallel; compressed files retain its streaming
decoders. SIMD indexing uses simd-minimizers 2.3.1 and its matching packed-seq and
seq-hash versions. `Cargo.lock` pins the complete dependency graph; the vendored
POASTA snapshot does not require a network Git checkout during a build.
