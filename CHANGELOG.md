# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- A `chainfile` command line tool, enabled with the `cli` feature, with
  `liftover bed`, `liftover vcf`, `view`, `stats`, and `check` subcommands.
- A `chainfile-generate` crate (in `tools/`) that generates synthetic genomes and
  chain files, and random BED and VCF records that exercise specific liftover
  situations (such as gaps, unaligned regions, multiple chains, and swapped
  alleles) at configurable percentages.
- A `chainfile-compare` crate (in `tools/`) that compares BED liftover with
  CrossMap.

### Changed

- The minimum supported Rust version is now 1.91.0 to support current
  transitive dependency releases and `noodles` 0.119.
- Dependencies are updated to their latest releases, including `nonempty`
  0.12, which is exposed by `alignment::Section::data()`, and `rust-lapper` 2.

### Fixed

- Malformed step-through iterators now terminate after reporting a data-section
  alignment error.
- Sequence construction now rejects alignment endpoints beyond the chromosome
  size.
- Liftover machine construction now returns an error for conflicting chromosome
  sizes.
- The alignment section iterator now stops after yielding an error instead of
  panicking when iteration continues.
- Liftover machine construction now ignores exact duplicate chain sections
  instead of returning their segments twice, and returns the new
  `InconsistentChainData` error for duplicate chain IDs with matching headers
  but different alignment data. This adds a variant to the public
  `liftover::machine::builder::Error` enum.
- Step-through iterators now stop after any error, and no longer log a
  warning about not being finished when dropped after an error.

### Removed

- The `chain_check`, `chain_count`, and `chain_view` examples, which are
  superseded by the `check`, `stats`, and `view` subcommands.
- The `compare-crossmap` binary and the `binaries` feature, which are
  superseded by `chainfile-compare crossmap`.

## 0.4.0 — 03-19-2026

### Added

- `Machine::liftover()` now returns results grouped by chain via `LiftoverResult`,
  making it possible to distinguish ambiguous mappings (multiple chains) from
  straddle splits (multiple segments within one chain)
  ([#9](https://github.com/stjude-rust-labs/chainfile/pull/9)).
- `LiftoverResult` exposes the full chain header via `chain()` and the mapping
  segments via `segments()` and `into_segments()`
  ([#9](https://github.com/stjude-rust-labs/chainfile/pull/9)).
- Duplicate chain IDs with inconsistent headers are now detected at build time
  in `liftover::machine::Builder`
  ([#9](https://github.com/stjude-rust-labs/chainfile/pull/9)).

### Changed

- All error types across the crate now use `thiserror`
  ([#9](https://github.com/stjude-rust-labs/chainfile/pull/9)).
- MSRV bumped to 1.85.0
  ([#9](https://github.com/stjude-rust-labs/chainfile/pull/9)).
- Bumped `omics` from `v0.2.0` to `v0.4.0`, adopting `Arc<str>`-backed
  `Contig` for `O(1)` clones, the new `Contig` validation API, and
  in-place coordinate move methods
  ([#10](https://github.com/stjude-rust-labs/chainfile/pull/10)).
- Optimized liftover hot paths: stepthrough now uses in-place
  `move_forward`/`move_backward` (cutting per-step clones from 8–12 to 4),
  lapper results are filtered by strand before cloning, `Contig` creation is
  hoisted out of the inner loop in the builder, and data record parsing uses
  `splitn` instead of `split().collect()`
  ([#10](https://github.com/stjude-rust-labs/chainfile/pull/10)).
- Refactored `AnnotatedPair` to store `chain_id: usize` instead of a full
  `Header`, with a `HashMap<usize, Header>` side table on `Machine`
  ([#10](https://github.com/stjude-rust-labs/chainfile/pull/10)).
- Liftover grouping replaced `HashMap` with a sort-and-linear-scan over a
  flat `Vec`
  ([#10](https://github.com/stjude-rust-labs/chainfile/pull/10)).
- `Machine::liftover()` now returns results in deterministic order: sorted by
  score descending (best chain first) then chain ID ascending, with segments
  sorted by reference start position
  ([#10](https://github.com/stjude-rust-labs/chainfile/pull/10)).
- Added Criterion benchmarks for single-position liftover, interval liftover,
  machine building, query-only throughput, and end-to-end throughput with
  fixed-seed RNG for reproducibility
  ([#10](https://github.com/stjude-rust-labs/chainfile/pull/10)).

## 0.3.0 - 01-03-2025

### Added

- Added a new `bin/compare-crossmap.rs` that compares the performance of `chainfile`
  with `CrossMap` (https://crossmap.sourceforge.net).

### Revise

- The inner data structures were replaced with facilities provided by the
  [`omics`](https://github.com/stjude-rust-labs/omics) family of crates
  ([#4](https://github.com/stjude-rust-labs/chainfile/pull/4)).
- In a subsequent pull request, the `omics` crates were updated to `v0.2.0`.
  This was a much more robust system for representing coordinates
  ([#5](https://github.com/stjude-rust-labs/chainfile/pull/5)).

## 0.2.1 — 09-09-2023

### Fixes

- Checks strand equality before performing liftover in `chainfile::liftover::Machine`.

## 0.2.0 — 04-24-2023

### Added

- Adds liftover capabilities and related examples.

## 0.1.0 — 03-30-2023

### Added

- Includes first version of the crate, including correct chain file parsing and
  iteration over alignment data sections within the chain file.
