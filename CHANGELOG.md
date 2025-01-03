# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

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
