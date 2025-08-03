//! This example is a very simple approach to testing the validity of a chain
//! file. To run the example, you'll need a reference FASTA, a query FASTA, as
//! well as a chain file lifting over from the reference to the query. You can
//! call the program like so:
//!
//! ```
//! cargo run --release --example chain_check <CHAIN> <REFERENCE_FASTA> <QUERY_FASTA>
//! ```
//!
//! If all goes well, the program will print out the total number of positions
//! compared, the total number of positions with mismatches, and the overall
//! match rate between the two genomes (via the chain file). In some cases, you
//! will experience errors: this indicates a malformed chain file or incorrect
//! pairings of chain file/reference FASTA/query FASTA.
//!
//! Note that this program somewhat simplifies what a "match" between the two
//! genomes means. For example, all nucleotides within both FASTA files are
//! mapped to 'A', 'C', 'G', 'T', or 'Other'. Lowercase and uppercase letters in
//! the FASTA file have semantic meaning, but for the purposes of this
//! comparison, they are considered the same.

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::iter::zip;

use chain::liftover;
use chainfile as chain;
use flate2::read::GzDecoder;
use noodles::fasta;
use noodles::fasta::record::Sequence;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    ////////////////////////////////////////////////////////////////////////////////////
    // Read arguments
    ////////////////////////////////////////////////////////////////////////////////////

    let chain_file_path = env::args().nth(1).expect("missing chainfile path");
    let reference_fasta_path = env::args().nth(2).expect("missing reference fasta path");
    let query_fasta_path = env::args().nth(3).expect("missing query fasta path");

    ////////////////////////////////////////////////////////////////////////////////////
    // Build the machine
    ////////////////////////////////////////////////////////////////////////////////////

    let chain = File::open(chain_file_path)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    let machine = liftover::machine::builder::Builder.try_build_from(chain)?;

    ////////////////////////////////////////////////////////////////////////////////////
    // Gather intervals to search
    ////////////////////////////////////////////////////////////////////////////////////

    let mut reference = HashMap::new();
    for result in fasta::reader::Builder
        .build_from_path(reference_fasta_path)?
        .records()
    {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).to_string();
        reference.insert(name, record.sequence().clone());
    }

    let mut query = HashMap::new();
    for result in fasta::reader::Builder
        .build_from_path(query_fasta_path)?
        .records()
    {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).to_string();
        query.insert(name, record.sequence().clone());
    }

    let intervals = reference
        .iter()
        .map(|(name, reference_sequence)| {
            Interval::try_new(
                Coordinate::new(name.as_str(), Strand::Positive, 0_u64),
                Coordinate::new(
                    name.as_str(),
                    Strand::Positive,
                    reference_sequence.len() as u64,
                ),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();

    ////////////////////////////////////////////////////////////////////////////////////
    // Do the comparison
    ////////////////////////////////////////////////////////////////////////////////////

    let mut total_positions = 0usize;
    let mut total_mismatches = 0usize;

    let mut positive_positions = 0usize;
    let mut positive_mismatches = 0usize;

    let mut negative_positions = 0usize;
    let mut negative_mismatches = 0usize;

    for interval in intervals {
        let reference_sequence = reference
            .get(interval.contig().as_str())
            .unwrap_or_else(|| {
                panic!(
                    "reference FASTA did not contain necessary contig \"{}\"; are the FASTA files \
                     and chain file correct?",
                    interval.contig().as_str()
                )
            });

        let results = match machine.liftover(interval.clone()) {
            Some(results) => results,
            None => {
                println!(
                    "INFO: region {interval} is not included in the query genome."
                );
                continue;
            }
        };

        for region in results {
            let query_strand = region.query().strand();

            let query_sequence = query
                .get(region.query().contig().as_str())
                .unwrap_or_else(|| {
                    panic!(
                        "query FASTA did not contain necessary contig \"{}\": are the FASTA files \
                         and chain file correct?",
                        region.query().contig()
                    )
                });

            let (reference_interval, query_interval) = region.into_parts();
            let reference = get_sequence_for_interval(reference_sequence, reference_interval);
            let query = get_sequence_for_interval(query_sequence, query_interval);

            assert_eq!(reference.len(), query.len());

            let positions = query.len();
            let mismatches = count_mismatches(&reference, &query);

            total_positions += positions;
            total_mismatches += mismatches;

            if query_strand == Strand::Positive {
                positive_positions += positions;
                positive_mismatches += mismatches;
            } else {
                negative_positions += positions;
                negative_mismatches += mismatches;
            }
        }
    }

    println!();

    let total_matched = (1.0 - (total_mismatches as f64 / total_positions as f64)) * 100.0;
    println!("Total mismatches found: {total_mismatches}");
    println!("Total positions checked: {total_positions}");
    println!("Percentage match for all regions: {total_matched:.1}%");

    println!();

    let positive_matched = (1.0 - (positive_mismatches as f64 / positive_positions as f64)) * 100.0;
    println!(
        "> Positive strand mismatches found: {positive_mismatches}"
    );
    println!(
        "> Positive strand positions checked: {positive_positions}"
    );
    println!(
        "> Percentage match for positive stranded query regions: {positive_matched:.1}%"
    );

    println!();

    let negative_matched = (1.0 - (negative_mismatches as f64 / negative_positions as f64)) * 100.0;
    println!(
        "> Negative strand mismatches found: {negative_mismatches}"
    );
    println!(
        "> Negative strand positions checked: {negative_positions}"
    );
    println!(
        "> Percentage match for negative stranded regions: {negative_matched:.1}%"
    );

    Ok(())
}

fn get_sequence_for_interval(sequence: &Sequence, interval: Interval) -> Vec<Nucleotide> {
    let strand = interval.strand();
    let interval = parse_interval(interval);

    let nucleotides = sequence
        .slice(interval)
        .expect("sequence lookup did not succeed: are the FASTA files and chain file correct?")
        .as_ref()
        .iter()
        .map(Nucleotide::from)
        .collect::<Vec<_>>();

    match strand {
        Strand::Positive => nucleotides,
        Strand::Negative => nucleotides
            .into_iter()
            .rev()
            .map(|n| n.complement())
            .collect::<Vec<_>>(),
    }
}

fn parse_interval(interval: Interval) -> noodles::core::region::Interval {
    let interval = interval.into_equivalent_base();

    let (start, end) = match interval.strand() {
        Strand::Positive => (
            interval.start().position().get() as usize,
            interval.end().position().get() as usize,
        ),
        Strand::Negative => (
            interval.end().position().get() as usize,
            interval.start().position().get() as usize,
        ),
    };

    let start = noodles::core::Position::new(start).unwrap();
    let end = noodles::core::Position::new(end).unwrap();

    noodles::core::region::Interval::from(start..=end)
}

#[derive(Debug, Eq, PartialEq)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    Other,
}

impl Nucleotide {
    pub fn complement(self) -> Nucleotide {
        match self {
            Nucleotide::A => Nucleotide::T,
            Nucleotide::C => Nucleotide::G,
            Nucleotide::G => Nucleotide::C,
            Nucleotide::T => Nucleotide::A,
            Nucleotide::Other => Nucleotide::Other,
        }
    }
}

impl From<&u8> for Nucleotide {
    fn from(value: &u8) -> Self {
        match value {
            b'A' => Nucleotide::A,
            b'a' => Nucleotide::A,
            b'C' => Nucleotide::C,
            b'c' => Nucleotide::C,
            b'G' => Nucleotide::G,
            b'g' => Nucleotide::G,
            b'T' => Nucleotide::T,
            b't' => Nucleotide::T,
            _ => Nucleotide::Other,
        }
    }
}

pub fn count_mismatches(reference: &[Nucleotide], query: &[Nucleotide]) -> usize {
    let mut result = 0;

    for (r, q) in zip(reference.iter(), query.iter()) {
        if r != q {
            result += 1;
        }
    }

    result
}
