//! A machine for lifting over coordinates within a reference genome.

use std::collections::HashMap;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interval::interbase::Interval;
use omics::coordinate::position::Number;
use rust_lapper as lapper;

use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;

pub mod builder;

pub use builder::Builder;

/// A dictionary of chromosome names and their size.
pub type ChromosomeDictionary = HashMap<String, Number>;

/// A machine for lifting over coordinates from a reference genome to a query
/// genome.
///
/// Generally, you will want to use a [`builder::Builder`] to construct one of
/// these.
#[derive(Debug)]
pub struct Machine {
    /// The inner lookup table of positions in the reference genome to positions
    /// in the query genome for each contig in the reference genome.
    inner: HashMap<Contig, lapper::Lapper<Number, ContiguousIntervalPair>>,

    /// The reference chromosome corpus.
    reference_chromosomes: ChromosomeDictionary,

    /// The query chromosome corpus.
    query_chromosomes: ChromosomeDictionary,
}

impl Machine {
    /// Gets all of the contigs of the reference genome along with the start and
    /// end positions for each contig.
    pub fn reference_chromosomes(&self) -> &ChromosomeDictionary {
        &self.reference_chromosomes
    }

    /// Gets all of the contigs of the query genome along with the start and
    /// end positions for each contig.
    pub fn query_chromosomes(&self) -> &ChromosomeDictionary {
        &self.query_chromosomes
    }

    /// Gets a reference to the inner hashmap.
    pub fn inner(&self) -> &HashMap<Contig, lapper::Lapper<Number, ContiguousIntervalPair>> {
        &self.inner
    }

    /// Performs a liftover from the specified `interval` to the query genome.
    pub fn liftover(&self, interval: Interval) -> Option<Vec<ContiguousIntervalPair>> {
        let entry = self.inner.get(interval.contig())?;

        let (start, stop) = match interval.strand() {
            Strand::Positive => {
                let start = interval.start().position().get();
                let end = interval.end().position().get();

                (start, end)
            }
            Strand::Negative => {
                let start = interval.end().position().get();
                let end = interval.start().position().get();

                (start, end)
            }
        };

        let strand = interval.strand();

        let results = entry
            .find(start, stop)
            .map(|e| e.val.clone())
            .filter(|i| i.reference().strand() == strand)
            .map(move |pair| pair.clamp(interval.clone()).unwrap())
            .collect::<Vec<_>>();

        match results.is_empty() {
            true => None,
            false => Some(results),
        }
    }
}

#[cfg(test)]
mod tests {
    use omics::coordinate::interbase::Coordinate;
    use omics::coordinate::position::interbase::Position;

    use super::*;
    use crate::Reader;
    use crate::liftover::machine;

    #[test]
    pub fn test_valid_positive_strand_liftover() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Positive, 1u64);
        let to = Coordinate::new("seq0", Strand::Positive, 7u64);

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().as_str(), "seq0");
        assert_eq!(result.reference().start().strand(), Strand::Positive);
        assert_eq!(result.reference().start().position(), &Position::new(1));
        assert_eq!(result.reference().end().strand(), Strand::Positive);
        assert_eq!(result.reference().end().position(), &Position::new(7));

        assert_eq!(result.query().contig().as_str(), "seq1");
        assert_eq!(result.query().start().strand(), Strand::Positive);
        assert_eq!(result.query().start().position(), &Position::new(1));
        assert_eq!(result.query().end().strand(), Strand::Positive);
        assert_eq!(result.query().end().position(), &Position::new(7));

        Ok(())
    }

    #[test]
    pub fn test_valid_negative_strand_liftover() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Negative, 7u64);
        let to = Coordinate::new("seq0", Strand::Negative, 1u64);

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().as_str(), "seq0");
        assert_eq!(result.reference().start().strand(), Strand::Negative);
        assert_eq!(result.reference().start().position(), &Position::new(7));
        assert_eq!(result.reference().end().strand(), Strand::Negative);
        assert_eq!(result.reference().end().position(), &Position::new(1));

        assert_eq!(result.query().contig().as_str(), "seq1");
        assert_eq!(result.query().start().strand(), Strand::Negative);
        assert_eq!(result.query().start().position(), &Position::new(7));
        assert_eq!(result.query().end().strand(), Strand::Negative);
        assert_eq!(result.query().end().position(), &Position::new(1));

        Ok(())
    }

    #[test]
    pub fn test_valid_positive_to_negative_strand_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Positive, 1u64);
        let to = Coordinate::new("seq0", Strand::Positive, 7u64);

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().as_str(), "seq0");
        assert_eq!(result.reference().start().strand(), Strand::Positive);
        assert_eq!(result.reference().start().position(), &Position::new(1));
        assert_eq!(result.reference().end().strand(), Strand::Positive);
        assert_eq!(result.reference().end().position(), &Position::new(7));

        assert_eq!(result.query().contig().as_str(), "seq1");
        assert_eq!(result.query().start().strand(), Strand::Negative);
        assert_eq!(result.query().start().position(), &Position::new(9));
        assert_eq!(result.query().end().strand(), Strand::Negative);
        assert_eq!(result.query().end().position(), &Position::new(3));

        Ok(())
    }

    #[test]
    pub fn test_valid_negative_to_positive_strand_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Negative, 7u64);
        let to = Coordinate::new("seq0", Strand::Negative, 1u64);

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().as_str(), "seq0");
        assert_eq!(result.reference().start().strand(), Strand::Negative);
        assert_eq!(result.reference().start().position(), &Position::new(7));
        assert_eq!(result.reference().end().strand(), Strand::Negative);
        assert_eq!(result.reference().end().position(), &Position::new(1));

        assert_eq!(result.query().contig().as_str(), "seq1");
        assert_eq!(result.query().start().strand(), Strand::Positive);
        assert_eq!(result.query().start().position(), &Position::new(3));
        assert_eq!(result.query().end().strand(), Strand::Positive);
        assert_eq!(result.query().end().position(), &Position::new(9));

        Ok(())
    }

    #[test]
    pub fn test_nonexistent_positive_strand_interval_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Positive, 1u64);
        let to = Coordinate::new("seq0", Strand::Positive, 7u64);

        let interval = Interval::try_new(from, to)?;
        let results = machine.liftover(interval);

        assert_eq!(results, None);

        Ok(())
    }

    #[test]
    pub fn test_nonexistent_negative_strand_interval_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Negative, 7u64);
        let to = Coordinate::new("seq0", Strand::Negative, 1u64);

        let interval = Interval::try_new(from, to)?;
        let results = machine.liftover(interval);

        assert_eq!(results, None);

        Ok(())
    }
}
