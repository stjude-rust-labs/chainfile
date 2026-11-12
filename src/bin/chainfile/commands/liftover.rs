//! The `liftover` subcommand.

pub mod bed;
pub mod vcf;

use anyhow::Result;
use clap::Subcommand;

/// Lift over records from the reference genome to the query genome.
#[derive(clap::Args, Debug)]
pub struct Args {
    /// The format of the records to lift over.
    #[command(subcommand)]
    format: Format,
}

/// The supported liftover formats.
///
/// NOTE: the help for each format comes from the doc comment on its `Args`
/// struct, so the variants are intentionally undocumented.
#[derive(Subcommand, Debug)]
enum Format {
    Bed(bed::Args),

    Vcf(vcf::Args),
}

/// Runs the `liftover` subcommand.
pub fn liftover(args: Args) -> Result<()> {
    match args.format {
        Format::Bed(args) => bed::liftover(args),
        Format::Vcf(args) => vcf::liftover(args),
    }
}
