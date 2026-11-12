<p align="center">
  <h1 align="center">
    <code>chainfile</code>
  </h1>

  <p align="center">
    <a href="https://github.com/stjude-rust-labs/chainfile/actions/workflows/CI.yml" target="_blank">
      <img alt="CI: Status" src="https://github.com/stjude-rust-labs/chainfile/actions/workflows/CI.yml/badge.svg" />
    </a>
    <a href="https://github.com/stjude-rust-labs/chainfile/actions/workflows/comparison.yml" target="_blank">
      <img alt="Comparison: Status" src="https://github.com/stjude-rust-labs/chainfile/actions/workflows/comparison.yml/badge.svg" />
    </a>
    <a href="https://crates.io/crates/chainfile" target="_blank">
      <img alt="crates.io version" src="https://img.shields.io/crates/v/chainfile">
    </a>
    <img alt="crates.io downloads" src="https://img.shields.io/crates/d/chainfile">
    <a href="https://github.com/stjude-rust-labs/chainfile/blob/master/LICENSE-APACHE" target="_blank">
      <img alt="License: Apache 2.0" src="https://img.shields.io/badge/license-Apache 2.0-blue.svg" />
    </a>
    <a href="https://github.com/stjude-rust-labs/chainfile/blob/master/LICENSE-MIT" target="_blank">
      <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg" />
    </a>
  </p>

  <p align="center">
    A crate for working with genomics chain files.
    <br />
    <a href="https://docs.rs/chainfile"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/stjude-rust-labs/chainfile/issues/new?assignees=&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    ·
    <a href="https://github.com/stjude-rust-labs/chainfile/issues/new?assignees=&title=Descriptive%20Title&labels=bug">Report Bug</a>
    ·
    ⭐ Consider starring the repo! ⭐
    <br />
  </p>
</p>

## Guiding Principles

This crate is written in the style of
[noodles](https://github.com/zaeleus/noodles), as it was originally intended to
be included as a pull request. After discussion with the maintainer of noodles,
we decided this should be its own, complimentary crate.

## 📚 Getting Started

To include this crate in your project, simply use the following command.

```bash
cargo add chainfile
```

You can take a look at the
[examples](https://github.com/stjude-rust-labs/chainfile/tree/main/examples) to
get a sense of how to use the crate.

## 🧬 Command Line Tool

The crate also provides a `chainfile` command line tool behind the `cli`
feature. It reads plain or gzip-compressed chain files and uses `-` to mean
standard input or output.

```bash
cargo install chainfile --features cli

# Lift BED records over; unmapped records are written with their reason.
chainfile liftover bed --chain hg19ToHg38.over.chain.gz input.bed \
  --output lifted.bed --unmapped unmapped.bed

# Lift VCF records over, checking REF alleles against the indexed query genome.
chainfile liftover vcf --chain hg19ToHg38.over.chain.gz --query-fasta hg38.fa \
  input.vcf --output lifted.vcf --rejected rejected.vcf

# Print alignment blocks, chain headers, or summary statistics.
chainfile view hg19ToHg38.over.chain.gz
chainfile view --headers hg19ToHg38.over.chain.gz
chainfile stats hg19ToHg38.over.chain.gz

# Check how well a chain file's aligned bases agree between two genomes.
chainfile check --chain hg19ToHg38.over.chain.gz \
  --reference-fasta hg19.fa --query-fasta hg38.fa --min-identity 0.95
```

FASTA files need a `.fai` index, plus a `.gzi` index when BGZF-compressed.

Both liftover commands use all available cores by default and keep records in
input order. Use `--threads` to set the number of worker threads.

### BED liftover

- Each lifted segment is written as its own record.
- The strand column is flipped when the alignment is on the opposite strand.
- Records with more than six columns are an error unless you pass
  `--include-bed-plus` (copy the extra columns) or `--ignore-bed-plus` (drop
  them).

### VCF liftover

- Records are rejected when they do not map uniquely and contiguously, or when
  REF does not match the query genome. Rejected records get a
  `LiftoverRejectReason` INFO field.
- On opposite-strand alignments, alleles are reverse-complemented and
  re-anchored.
- When a biallelic substitution's ALT matches the query genome, REF and ALT are
  swapped and the record gets a `LiftoverSwappedAlleles` INFO flag.
  Allele-indexed fields (such as `GT`, `AD`, `PL`, `AF`, and `AC`) are
  rewritten, and fields that cannot be rewritten are dropped.

Run `chainfile help <command>` for the full set of options.

## 🖥️ Development

To bootstrap a development environment, please use the following commands.

```bash
# Clone the repository
git clone git@github.com:stjude-rust-labs/chainfile.git
cd chainfile

# Build the crate in release mode
cargo build --release

# List out the examples
cargo run --release --example
```

### Generating test data

The `chainfile-generate` crate in `tools/` generates random test data.

- `fixture` writes a synthetic reference genome, query genome, and chain file.
- `bed` and `vcf` write random records for a chain file. Flags such as
  `--spans-gap 5` set the percentage of records in each situation.

Each run prints its seed. Pass it back with `--seed` to regenerate the same
data.

```bash
cargo run --release -p chainfile-generate -- fixture -o fixture
cargo run --release -p chainfile-generate -- bed -c fixture/reference-to-query.chain \
    -n 1000 --spans-gap 5 --unaligned 5 --multiple-chains 5 --opposite-strand 5
```

### Comparing with other tools

The `chainfile-compare` crate in `tools/` lifts over random records with both
`chainfile` and another tool. It prints a table of outcomes and fails if any
record disagrees.

```bash
cargo build --release --workspace --all-features

# BED liftover with CrossMap (`pip install CrossMap`).
target/release/chainfile-compare --chainfile target/release/chainfile \
    crossmap bed -c fixture/reference-to-query.chain
```

The `Compare` workflow runs these on every push. Each run uses a new random
fixture and new random records, along with UCSC's `hg19ToHg38` and `hg38ToHs1`
chain files. The seeds are printed in the logs so failures can be reproduced.

## 🚧️ Tests

Before submitting any pull requests, please make sure the code passes the
following checks.

```bash
# Run the project's tests.
cargo test --workspace --all-features

# Ensure the project doesn't have any linting warnings.
cargo clippy --workspace --all-targets --all-features

# Ensure the project passes `cargo fmt`.
cargo +nightly fmt --check
```

## Minimum Supported Rust Version (MSRV)

This crate is designed to work with Rust version 1.91.0 or later. It may, by
happenstance, work with earlier versions of Rust.

## 🤝 Contributing

Contributions, issues and feature requests are welcome! Feel free to check
[issues page](https://github.com/stjude-rust-labs/chainfile/issues).

## 📝 License

This project is licensed as either [Apache 2.0][license-apache] or
[MIT][license-mit] at your discretion. Additionally, please see [the
disclaimer](https://github.com/stjude-rust-labs#disclaimer) that applies to all
crates and command line tools made available by St. Jude Rust Labs.

Copyright © 2023-Present [St. Jude Children's Research Hospital](https://github.com/stjude).

[license-apache]: https://github.com/stjude-rust-labs/chainfile/blob/master/LICENSE-APACHE
[license-mit]: https://github.com/stjude-rust-labs/chainfile/blob/master/LICENSE-MIT
