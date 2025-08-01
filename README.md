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

## 🚧️ Tests

Before submitting any pull requests, please make sure the code passes the
following checks.

```bash
# Run the project's tests.
cargo test --all-features

# Ensure the project doesn't have any linting warnings.
cargo clippy --all-features

# Ensure the project passes `cargo fmt`.
cargo fmt --check
```

## Minimum Supported Rust Version (MSRV)

This crate is designed to work with Rust version 1.82.0 or later. It may, by
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
