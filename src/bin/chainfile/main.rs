//! A command line tool for inspecting chain files and lifting over genomic
//! data between genome builds.
//!
//! ```shell
//! cargo run --release --features=cli --bin=chainfile -- --help
//! ```

mod commands;
mod fasta;
mod io;
mod sequence;

use anyhow::Result;
use clap::Parser;
use clap::Subcommand;
use clap_verbosity_flag::InfoLevel;
use clap_verbosity_flag::Verbosity;
use tracing_log::AsTrace as _;
use tracing_subscriber::EnvFilter;

/// Inspect chain files and lift over genomic data between genome builds.
#[derive(Parser, Debug)]
#[command(name = "chainfile", version, propagate_version = true)]
struct Cli {
    /// The subcommand to run.
    #[command(subcommand)]
    command: Command,

    /// The logging verbosity.
    #[command(flatten)]
    verbosity: Verbosity<InfoLevel>,
}

/// The available subcommands.
///
/// NOTE: the help for each subcommand comes from the doc comment on its `Args`
/// struct, so the variants are intentionally undocumented.
#[derive(Subcommand, Debug)]
enum Command {
    Liftover(commands::liftover::Args),

    View(commands::view::Args),

    Stats(commands::stats::Args),

    Check(commands::check::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // NOTE: logs go to stderr so that they never mix with data on stdout.
    match std::env::var("RUST_LOG") {
        Ok(_) => tracing_subscriber::fmt()
            .with_env_filter(EnvFilter::from_default_env())
            .with_writer(std::io::stderr)
            .init(),
        Err(_) => tracing_subscriber::fmt()
            .with_max_level(cli.verbosity.log_level_filter().as_trace())
            .with_writer(std::io::stderr)
            .init(),
    };

    match cli.command {
        Command::Liftover(args) => commands::liftover::liftover(args),
        Command::View(args) => commands::view::view(args),
        Command::Stats(args) => commands::stats::stats(args),
        Command::Check(args) => commands::check::check(args),
    }
}
