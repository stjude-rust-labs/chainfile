//! A machine for lifting over coordinates within a reference genome.

use std::collections::HashMap;

use rust_lapper as lapper;

use crate::core::coordinate::Position;
use crate::core::Contig;
use crate::core::Interval;
use crate::core::Strand;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;

pub mod builder;

pub use builder::Builder;

/// A machine for lifting over coordinates from a reference genome to a query
/// genome.
///
/// Generally, you will want to use a [`builder::Builder`] to construct one of
/// these.
pub struct Machine {
    inner: HashMap<Contig, lapper::Lapper<usize, ContiguousIntervalPair>>,
}

impl Machine {
    /// Performs a liftover from the specified `interval` to the query genome.
    pub fn liftover(&self, interval: &Interval) -> Option<Vec<ContiguousIntervalPair>> {
        let entry = self.inner.get(interval.contig())?;

        let (start, stop) = match interval.strand() {
            Strand::Positive => {
                let start = match interval.start().position() {
                    Position::ZeroBased(a) => *a,
                    Position::NegativeBound => {
                        unreachable!("negative bound not allowed on positive strand")
                    }
                };

                let stop = match interval.end().position() {
                    Position::ZeroBased(b) => *b,
                    Position::NegativeBound => {
                        unreachable!("negative bound not allowed on positive strand")
                    }
                };

                (start, stop)
            }
            Strand::Negative => {
                let start = match interval.end().position() {
                    Position::ZeroBased(a) => a + 1,
                    Position::NegativeBound => 0,
                };

                let stop = match interval.start().position() {
                    Position::ZeroBased(b) => b + 1,
                    Position::NegativeBound => unreachable!(
                        "negative bound will never be the start of a negative-stranded interval"
                    ),
                };

                (start, stop)
            }
        };

        let results = entry
            .find(start, stop)
            .map(|e| e.val.clone())
            .filter(|i| i.reference().strand() == interval.strand())
            .map(|pair| pair.clamp(interval).unwrap())
            .collect::<Vec<_>>();

        match results.is_empty() {
            true => None,
            false => Some(results),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::Coordinate;
    use crate::core::Interval;
    use crate::core::Position;
    use crate::core::Strand;
    use crate::liftover::machine;
    use crate::Reader;

    #[test]
    pub fn test_valid_positive_strand_liftover() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder::default().try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", 1, Strand::Positive)?;
        let to = Coordinate::try_new("seq0", 7, Strand::Positive)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Positive);
        assert_eq!(result.reference().start().position(), &Position::from(1));
        assert_eq!(result.reference().end().strand(), &Strand::Positive);
        assert_eq!(result.reference().end().position(), &Position::from(7));

        assert_eq!(result.query().contig(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Positive);
        assert_eq!(result.query().start().position(), &Position::from(1));
        assert_eq!(result.query().end().strand(), &Strand::Positive);
        assert_eq!(result.query().end().position(), &Position::from(7));

        Ok(())
    }

    #[test]
    pub fn test_valid_negative_strand_liftover() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder::default().try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", 7, Strand::Negative)?;
        let to = Coordinate::try_new("seq0", 1, Strand::Negative)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Negative);
        assert_eq!(result.reference().start().position(), &Position::from(7));
        assert_eq!(result.reference().end().strand(), &Strand::Negative);
        assert_eq!(result.reference().end().position(), &Position::from(1));

        assert_eq!(result.query().contig(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Negative);
        assert_eq!(result.query().start().position(), &Position::from(7));
        assert_eq!(result.query().end().strand(), &Strand::Negative);
        assert_eq!(result.query().end().position(), &Position::from(1));

        Ok(())
    }

    #[test]
    pub fn test_valid_positive_to_negative_strand_liftover(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder::default().try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", 1, Strand::Positive)?;
        let to = Coordinate::try_new("seq0", 7, Strand::Positive)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Positive);
        assert_eq!(result.reference().start().position(), &Position::from(1));
        assert_eq!(result.reference().end().strand(), &Strand::Positive);
        assert_eq!(result.reference().end().position(), &Position::from(7));

        assert_eq!(result.query().contig(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Negative);
        assert_eq!(result.query().start().position(), &Position::from(8));
        assert_eq!(result.query().end().strand(), &Strand::Negative);
        assert_eq!(result.query().end().position(), &Position::from(2));

        Ok(())
    }

    #[test]
    pub fn test_valid_negative_to_positive_strand_liftover(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder::default().try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", 7, Strand::Negative)?;
        let to = Coordinate::try_new("seq0", 1, Strand::Negative)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Negative);
        assert_eq!(result.reference().start().position(), &Position::from(7));
        assert_eq!(result.reference().end().strand(), &Strand::Negative);
        assert_eq!(result.reference().end().position(), &Position::from(1));

        assert_eq!(result.query().contig(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Positive);
        assert_eq!(result.query().start().position(), &Position::from(2));
        assert_eq!(result.query().end().strand(), &Strand::Positive);
        assert_eq!(result.query().end().position(), &Position::from(8));

        Ok(())
    }

    #[test]
    pub fn test_nonexistent_positive_strand_interval_liftover(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder::default().try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", 1, Strand::Positive)?;
        let to = Coordinate::try_new("seq0", 7, Strand::Positive)?;

        let interval = Interval::try_new(from, to)?;
        let results = machine.liftover(&interval);

        assert_eq!(results, None);

        Ok(())
    }

    #[test]
    pub fn test_nonexistent_negative_strand_interval_liftover(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder::default().try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", 7, Strand::Negative)?;
        let to = Coordinate::try_new("seq0", 1, Strand::Negative)?;

        let interval = Interval::try_new(from, to)?;
        let results = machine.liftover(&interval);

        assert_eq!(results, None);

        Ok(())
    }
}
