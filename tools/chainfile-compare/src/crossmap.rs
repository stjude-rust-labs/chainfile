//! Comparisons with CrossMap.

use std::path::PathBuf;

use anyhow::Result;
use chainfile_generate::bed;
use clap::Subcommand;

use crate::bed::compare_exact;
use crate::bed::describe;
use crate::bed::read_names;
use crate::bed::read_segments;
use crate::context::Context;
use crate::context::create;
use crate::context::display;
use crate::context::execute;
use crate::report::Report;

/// Compares `chainfile` with CrossMap.
#[derive(clap::Args)]
pub struct Args {
    /// The format.
    #[command(subcommand)]
    format: Format,
}

/// A format.
// NOTE: the variants are intentionally undocumented, as a variant's doc
// comment would override the help of its arguments struct.
#[derive(Subcommand)]
enum Format {
    #[allow(missing_docs)]
    Bed(BedArgs),
}

/// Compares BED liftover with CrossMap.
///
/// Random BED6 records are lifted over with `chainfile liftover bed` and
/// `CrossMap bed`. Both tools report one record per aligned block, so the
/// lifted records for each input record must be identical.
#[derive(clap::Args)]
#[command(
    mut_arg("spans_gap", |arg| arg.default_value("5")),
    mut_arg("unaligned", |arg| arg.default_value("5")),
    mut_arg("multiple_chains", |arg| arg.default_value("5")),
    mut_arg("opposite_strand", |arg| arg.default_value("5")),
)]
pub struct BedArgs {
    /// The chain file (which may be gzip-compressed).
    #[arg(short, long)]
    chain: PathBuf,

    /// The CrossMap executable.
    #[arg(long, default_value = "CrossMap")]
    crossmap: PathBuf,

    /// The options for generating records.
    #[command(flatten)]
    options: bed::Options,
}

/// Runs a comparison with CrossMap.
pub fn run(context: &mut Context, args: Args) -> Result<()> {
    match args.format {
        Format::Bed(args) => compare_bed(context, args),
    }
}

/// Compares BED liftover with CrossMap.
fn compare_bed(context: &mut Context, args: BedArgs) -> Result<()> {
    let input = context.path("input.bed");
    create(&input, |writer| {
        bed::generate(&args.chain, &args.options, &mut context.rng, writer)
    })?;

    let lifted = context.path("chainfile.bed");
    let unmapped = context.path("chainfile.unmapped.bed");
    execute(
        context
            .chainfile()
            .args(["liftover", "bed", "-c"])
            .arg(&args.chain)
            .arg("-o")
            .arg(&lifted)
            .arg("-u")
            .arg(&unmapped)
            .arg(&input),
    )?;

    let crossmap = context.path("crossmap.bed");
    execute(
        std::process::Command::new(&args.crossmap)
            .arg("bed")
            .arg(&args.chain)
            .arg(&input)
            .arg(&crossmap),
    )?;

    let names = read_names(&input)?;
    let lifted = read_segments(&lifted)?;
    let crossmap = read_segments(&crossmap)?;

    let mut report = Report::new("crossmap", context.examples);

    for name in &names {
        let ours = lifted.get(name).map(Vec::as_slice).unwrap_or_default();
        let theirs = crossmap.get(name).map(Vec::as_slice).unwrap_or_default();

        report.record(name, compare_exact(ours, theirs), || {
            format!(
                "chainfile: {}; crossmap: {}",
                describe(ours),
                describe(theirs)
            )
        });
    }

    report.finish(&format!(
        "`chainfile` vs. CrossMap: {} BED records lifted over with `{}`",
        names.len(),
        display(&args.chain)
    ))
}
