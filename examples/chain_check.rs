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
use omics::coordinate::Strand;
use omics::coordinate::interval::zero::Interval;
use omics::coordinate::position::Value;
use omics::coordinate::zero::Coordinate;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let chain_file_path = env::args().nth(1).expect("missing chainfile path");
    let reference_fasta_path = env::args().nth(2).expect("missing reference fasta path");
    let query_fasta_path = env::args().nth(3).expect("missing query fasta path");

    let chain = File::open(chain_file_path)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chain::Reader::new)?;

    let machine = liftover::machine::builder::Builder.try_build_from(chain)?;

    let mut reference = HashMap::new();
    for result in fasta::reader::Builder
        .build_from_path(reference_fasta_path)?
        .records()
    {
        let record = result?;
        reference.insert(record.name().to_string(), record.sequence().clone());
    }

    let mut query = HashMap::new();
    for result in fasta::reader::Builder
        .build_from_path(query_fasta_path)?
        .records()
    {
        let record = result?;
        query.insert(record.name().to_string(), record.sequence().clone());
    }

    let mut total_positions = 0usize;
    let mut mismatches = 0usize;

    for (name, reference_sequence) in reference {
        let interval = Interval::try_new(
            Coordinate::try_new(name.clone(), Strand::Positive, 0)?,
            Coordinate::try_new(name, Strand::Positive, reference_sequence.len())?,
        )?;

        let results = match machine.liftover(&interval) {
            Some(results) => results,
            None => {
                println!(
                    "INFO: region {} is not included in the query genome.",
                    interval
                );
                continue;
            }
        };

        for region in results {
            let query_sequence = query
                .get(region.query().contig().inner())
                .unwrap_or_else(|| {
                    panic!(
                        "query FASTA did not contain necessary contig \"{}\": are the FASTA files \
                         and chain file correct?",
                        region.query().contig()
                    )
                });

            let reference = get_sequence_for_interval(&reference_sequence, region.reference());
            let query = get_sequence_for_interval(query_sequence, region.query());

            assert_eq!(reference.len(), query.len());
            total_positions += reference.len();
            mismatches += count_mismatches(&reference, &query);
        }
    }

    let result = (1.0 - (mismatches as f64 / total_positions as f64)) * 100.0;
    println!("Total mismatches found: {}", mismatches);
    println!("Total positions checked: {}", total_positions);
    println!("Percentage match for lifted regions: {:.1}%", result);

    Ok(())
}

fn get_sequence_for_interval(
    sequence: &fasta::record::Sequence,
    interval: &Interval,
) -> Vec<Nucleotide> {
    let result = sequence
        .slice(parse_interval(interval))
        .expect("sequence lookup did not succeed: are the FASTA files and chain file correct?")
        .as_ref()
        .iter()
        .map(Nucleotide::from)
        .collect::<Vec<_>>();

    match interval.strand() {
        Strand::Positive => result,
        Strand::Negative => result
            .into_iter()
            .rev()
            .map(|n| n.complement())
            .collect::<Vec<_>>(),
    }
}

fn parse_interval(interval: &Interval) -> noodles::core::region::Interval {
    let (start, end) = match interval.strand() {
        Strand::Positive => {
            let start = match interval.start().position().inner() {
                Value::Usize(position) => noodles::core::Position::try_from(*position + 1).unwrap(),
                Value::LowerBound => unreachable!(),
            };

            let end = match interval.end().position().inner() {
                Value::Usize(position) => noodles::core::Position::try_from(*position).unwrap(),
                Value::LowerBound => unreachable!(),
            };

            (start, end)
        }
        Strand::Negative => {
            let start = match interval.end().position().inner() {
                Value::Usize(position) => noodles::core::Position::try_from(*position + 2).unwrap(),
                Value::LowerBound => noodles::core::Position::try_from(1).unwrap(),
            };

            let end = match interval.start().position().inner() {
                Value::Usize(position) => noodles::core::Position::try_from(*position + 1).unwrap(),
                Value::LowerBound => unreachable!(),
            };

            (start, end)
        }
    };

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
