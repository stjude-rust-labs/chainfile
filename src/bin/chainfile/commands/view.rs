//! The `view` subcommand.

use std::io::Write as _;
use std::path::PathBuf;

use anyhow::Context as _;
use anyhow::Result;

use crate::io::open_chain;
use crate::io::open_output;
use crate::sequence::forward_bounds;

/// Print the alignment blocks within a chain file.
///
/// Alignment blocks are printed as tab-separated values. Coordinates are
/// 0-based, half-open, and numbered from the start of the forward strand (as
/// in BED), regardless of the strand of the alignment.
#[derive(clap::Args, Debug)]
pub struct Args {
    /// The chain file (optionally gzipped); `-` reads from standard input.
    chain: PathBuf,

    /// Only show the chain with this ID; may be given more than once.
    #[arg(long = "id", value_name = "ID")]
    ids: Vec<usize>,

    /// Print the chain header lines instead of the alignment blocks.
    #[arg(long)]
    headers: bool,

    /// The output file; standard output is used when not provided.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// The columns of the alignment block table.
const COLUMNS: &[&str] = &[
    "chain_id",
    "reference_contig",
    "reference_strand",
    "reference_start",
    "reference_end",
    "query_contig",
    "query_strand",
    "query_start",
    "query_end",
    "size",
    "dt",
    "dq",
];

/// Runs the `view` subcommand.
pub fn view(args: Args) -> Result<()> {
    let mut reader = open_chain(&args.chain)?;
    let mut output = open_output(args.output.as_deref())?;

    if !args.headers {
        writeln!(output, "{}", COLUMNS.join("\t"))?;
    }

    for result in reader.sections() {
        let section = result.context("reading chain file")?;
        let header = section.header();

        if !args.ids.is_empty() && !args.ids.contains(&header.id()) {
            continue;
        }

        if args.headers {
            writeln!(output, "{header}")?;
            continue;
        }

        let steps = section
            .stepthrough_with_data()
            .with_context(|| format!("stepping through chain {}", header.id()))?;

        for result in steps {
            let (pair, data) =
                result.with_context(|| format!("stepping through chain {}", header.id()))?;
            let (reference_start, reference_end) = forward_bounds(pair.reference());
            let (query_start, query_end) = forward_bounds(pair.query());

            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                header.id(),
                pair.reference().contig(),
                pair.reference().strand(),
                reference_start,
                reference_end,
                pair.query().contig(),
                pair.query().strand(),
                query_start,
                query_end,
                data.size(),
                optional(data.dt()),
                optional(data.dq()),
            )?;
        }
    }

    output.flush()?;
    Ok(())
}

/// Formats an optional value, using `.` when the value is missing.
fn optional(value: Option<u64>) -> String {
    value.map(|v| v.to_string()).unwrap_or_else(|| ".".into())
}
