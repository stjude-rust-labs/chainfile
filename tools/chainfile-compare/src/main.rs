//! Compares `chainfile` liftover results with other liftover tools.

mod bed;
mod context;
mod crossmap;
mod report;

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use clap::Subcommand;

use crate::context::Context;

/// Compares `chainfile` liftover results with other liftover tools.
///
/// Random records are generated for the reference genome of a chain file (see
/// `chainfile-generate`), lifted over with both the `chainfile` command line
/// tool and another tool, and the results are compared record by record. A
/// table of outcomes by situation is printed, and the command fails if any
/// record disagrees.
#[derive(Parser)]
#[command(version)]
struct Cli {
    /// The seed for the random number generator.
    ///
    /// A random seed is chosen (and printed) when omitted.
    #[arg(long, global = true)]
    seed: Option<u64>,

    /// The `chainfile` executable.
    #[arg(long, global = true, default_value = "chainfile")]
    chainfile: PathBuf,

    /// A directory in which to keep the generated and lifted files.
    ///
    /// A temporary directory that is removed afterward is used when omitted.
    #[arg(long, global = true)]
    work_dir: Option<PathBuf>,

    /// The maximum number of disagreements to print.
    #[arg(long, global = true, default_value_t = 10)]
    examples: usize,

    /// The command.
    #[command(subcommand)]
    command: Command,
}

/// A command.
// NOTE: the variants are intentionally undocumented, as a variant's doc
// comment would override the help of its arguments struct.
#[derive(Subcommand)]
enum Command {
    #[allow(missing_docs)]
    Crossmap(crossmap::Args),
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let result = Context::new(cli.seed, cli.chainfile, cli.work_dir, cli.examples).and_then(
        |mut context| match cli.command {
            Command::Crossmap(args) => crossmap::run(&mut context, args),
        },
    );

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error:?}");
            ExitCode::FAILURE
        }
    }
}
