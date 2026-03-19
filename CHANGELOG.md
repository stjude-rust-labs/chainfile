# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

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
