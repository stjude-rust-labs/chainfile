//! The `stats` subcommand.

use std::collections::HashSet;
use std::io::Write as _;
use std::path::PathBuf;

use anyhow::Context as _;
use anyhow::Result;
use omics::coordinate::Strand;

use crate::io::open_chain;
use crate::io::open_output;

/// Print summary statistics for a chain file.
///
/// Statistics are printed as tab-separated `name` and `value` pairs.
#[derive(clap::Args, Debug)]
pub struct Args {
    /// The chain file (optionally gzipped); `-` reads from standard input.
    chain: PathBuf,

    /// The output file; standard output is used when not provided.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// Runs the `stats` subcommand.
pub fn stats(args: Args) -> Result<()> {
    let mut reader = open_chain(&args.chain)?;

    let mut chains = 0u64;
    let mut positive_chains = 0u64;
    let mut negative_chains = 0u64;
    let mut blocks = 0u64;
    let mut aligned_bases = 0u64;
    let mut reference_gap_bases = 0u64;
    let mut query_gap_bases = 0u64;
    let mut reference_contigs = HashSet::new();
    let mut query_contigs = HashSet::new();

    for result in reader.sections() {
        let section = result.context("reading chain file")?;
        let header = section.header();

        chains += 1;
        match header.query_sequence().strand() {
            Strand::Positive => positive_chains += 1,
            Strand::Negative => negative_chains += 1,
        }

        reference_contigs.insert(header.reference_sequence().chromosome_name().to_string());
        query_contigs.insert(header.query_sequence().chromosome_name().to_string());

        for record in section.data().iter() {
            blocks += 1;
            aligned_bases += record.size();
            reference_gap_bases += record.dt().unwrap_or_default();
            query_gap_bases += record.dq().unwrap_or_default();
        }
    }

    let mut output = open_output(args.output.as_deref())?;
    writeln!(output, "name\tvalue")?;
    writeln!(output, "chains\t{chains}")?;
    writeln!(output, "positive_strand_chains\t{positive_chains}")?;
    writeln!(output, "negative_strand_chains\t{negative_chains}")?;
    writeln!(output, "blocks\t{blocks}")?;
    writeln!(output, "aligned_bases\t{aligned_bases}")?;
    writeln!(output, "reference_gap_bases\t{reference_gap_bases}")?;
    writeln!(output, "query_gap_bases\t{query_gap_bases}")?;
    writeln!(output, "reference_contigs\t{}", reference_contigs.len())?;
    writeln!(output, "query_contigs\t{}", query_contigs.len())?;
    output.flush()?;

    Ok(())
}
