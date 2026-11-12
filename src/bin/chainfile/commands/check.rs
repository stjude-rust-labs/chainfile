//! The `check` subcommand.

use std::io::Write as _;
use std::path::PathBuf;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use omics::coordinate::Strand;
use omics::coordinate::interval::interbase::Interval;

use crate::fasta::IndexedFasta;
use crate::io::open_chain;
use crate::io::open_output;
use crate::sequence::forward_bounds;
use crate::sequence::oriented;

/// Check how well the sequences aligned by a chain file agree.
///
/// Every aligned base in the chain file is compared between the reference and
/// query genomes. Only positions where both bases are unambiguous (`A`, `C`,
/// `G`, or `T`, ignoring case) are compared. The results are printed as
/// tab-separated values for all chains and for chains on each query strand.
#[derive(clap::Args, Debug)]
pub struct Args {
    /// The chain file (optionally gzipped); `-` reads from standard input.
    #[arg(short, long)]
    chain: PathBuf,

    /// The reference genome FASTA (the genome lifted over from), which must
    /// have a `.fai` index (and a `.gzi` index when BGZF-compressed).
    #[arg(long)]
    reference_fasta: PathBuf,

    /// The query genome FASTA (the genome lifted over to), which must have a
    /// `.fai` index (and a `.gzi` index when BGZF-compressed).
    #[arg(long)]
    query_fasta: PathBuf,

    /// Exit with an error if the overall identity is below this value
    /// (between 0 and 1).
    #[arg(long, value_parser = parse_identity)]
    min_identity: Option<f64>,

    /// The output file; standard output is used when not provided.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// Parses an identity threshold.
fn parse_identity(value: &str) -> Result<f64, String> {
    match value.parse::<f64>() {
        Ok(v) if (0.0..=1.0).contains(&v) => Ok(v),
        _ => Err(format!("`{value}` is not a number between 0 and 1")),
    }
}

/// Comparison counts.
#[derive(Clone, Copy, Debug, Default)]
struct Counts {
    /// The number of aligned bases.
    aligned: u64,

    /// The number of aligned bases where both bases are unambiguous.
    compared: u64,

    /// The number of compared bases that differ.
    mismatches: u64,
}

impl Counts {
    /// Gets the fraction of compared bases that match.
    fn identity(&self) -> Option<f64> {
        (self.compared > 0).then(|| 1.0 - (self.mismatches as f64 / self.compared as f64))
    }

    /// Adds counts to these counts.
    fn add(&mut self, other: Counts) {
        self.aligned += other.aligned;
        self.compared += other.compared;
        self.mismatches += other.mismatches;
    }
}

/// Runs the `check` subcommand.
pub fn check(args: Args) -> Result<()> {
    let mut reader = open_chain(&args.chain)?;
    let mut reference = IndexedFasta::open(&args.reference_fasta)?;
    let mut query = IndexedFasta::open(&args.query_fasta)?;

    let mut positive = Counts::default();
    let mut negative = Counts::default();

    for result in reader.sections() {
        let section = result.context("reading chain file")?;
        let id = section.header().id();

        let steps = section
            .stepthrough()
            .with_context(|| format!("stepping through chain {id}"))?;

        for result in steps {
            let pair = result.with_context(|| format!("stepping through chain {id}"))?;
            let reference_bases = fetch(&mut reference, pair.reference(), "reference")?;
            let query_bases = fetch(&mut query, pair.query(), "query")?;
            let counts = compare(&reference_bases, &query_bases);

            match pair.query().strand() {
                Strand::Positive => positive.add(counts),
                Strand::Negative => negative.add(counts),
            }
        }
    }

    let mut all = positive;
    all.add(negative);

    let mut output = open_output(args.output.as_deref())?;
    writeln!(
        output,
        "query_strand\taligned_bases\tcompared_bases\tmismatches\tidentity"
    )?;

    for (name, counts) in [("all", all), ("+", positive), ("-", negative)] {
        let identity = counts
            .identity()
            .map(|v| format!("{v:.6}"))
            .unwrap_or_else(|| String::from("NA"));

        writeln!(
            output,
            "{name}\t{}\t{}\t{}\t{identity}",
            counts.aligned, counts.compared, counts.mismatches
        )?;
    }

    output.flush()?;

    let Some(identity) = all.identity() else {
        bail!("no unambiguous aligned bases were compared");
    };

    if let Some(min_identity) = args.min_identity {
        if identity < min_identity {
            bail!("identity {identity:.6} is below the minimum of {min_identity}");
        }
    }

    Ok(())
}

/// Fetches the oriented bases for an interval.
fn fetch(fasta: &mut IndexedFasta, interval: &Interval, genome: &str) -> Result<Vec<u8>> {
    let contig = interval.contig().as_str();
    let (start, end) = forward_bounds(interval);

    let bases = fasta
        .fetch(contig, start, end)
        .with_context(|| format!("fetching {interval} from the {genome} FASTA"))?
        .with_context(|| {
            format!("contig `{contig}` is not in the {genome} FASTA; are the files correct?")
        })?;

    Ok(oriented(bases, interval.strand()))
}

/// Compares two aligned sequences of uppercased bases.
fn compare(reference: &[u8], query: &[u8]) -> Counts {
    let mut counts = Counts {
        aligned: reference.len() as u64,
        ..Default::default()
    };

    for (r, q) in reference.iter().zip(query) {
        if matches!(r, b'A' | b'C' | b'G' | b'T') && matches!(q, b'A' | b'C' | b'G' | b'T') {
            counts.compared += 1;

            if r != q {
                counts.mismatches += 1;
            }
        }
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compare_skips_ambiguous_bases() {
        let counts = compare(b"ACGTN", b"ACCTA");
        assert_eq!(counts.aligned, 5);
        assert_eq!(counts.compared, 4);
        assert_eq!(counts.mismatches, 1);
        assert_eq!(counts.identity(), Some(0.75));
    }

    #[test]
    fn identity_is_none_without_comparisons() {
        assert_eq!(compare(b"NN", b"AC").identity(), None);
    }

    #[test]
    fn identity_threshold_is_validated() {
        assert_eq!(parse_identity("0.5"), Ok(0.5));
        assert!(parse_identity("1.5").is_err());
        assert!(parse_identity("NaN").is_err());
        assert!(parse_identity("abc").is_err());
    }
}
