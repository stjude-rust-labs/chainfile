//! A machine for lifting over coordinates within a reference genome.

use std::collections::HashMap;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interval::zero::Interval;
use rust_lapper as lapper;

use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;

pub mod builder;

pub use builder::Builder;

/// A machine for lifting over coordinates from a reference genome to a query
/// genome.
///
/// Generally, you will want to use a [`builder::Builder`] to construct one of
/// these.
#[derive(Debug)]
pub struct Machine {
    /// The inner lookup table of positions in the reference genome to positions
    /// in the query genome for each contig in the reference genome.
    inner: HashMap<Contig, lapper::Lapper<usize, ContiguousIntervalPair>>,
}

impl Machine {
    /// Performs a liftover from the specified `interval` to the query genome.
    pub fn liftover(&self, interval: &Interval) -> Option<Vec<ContiguousIntervalPair>> {
        let entry = self.inner.get(interval.contig())?;

        let (start, stop) = match interval.strand() {
            Strand::Positive => {
                let start = interval
                    .start()
                    .position()
                    .inner()
                    .get()
                    .unwrap_or_else(|| unreachable!("lower bound not allowed on positive strand"));

                let end =
                    interval.end().position().inner().get().unwrap_or_else(|| {
                        unreachable!("lower bound not allowed on positive strand")
                    });

                (start, end)
            }
            Strand::Negative => {
                let start = interval
                    .end()
                    .position()
                    .get()
                    .map(|pos| pos + 1)
                    .unwrap_or(0);

                let end = interval
                    .start()
                    .position()
                    .get()
                    .map(|pos| pos + 1)
                    .unwrap_or_else(|| {
                        unreachable!(
                            "lower bound will never be the start of a negative-stranded interval"
                        )
                    });

                (start, end)
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
    use omics::coordinate::position::zero::Position;
    use omics::coordinate::zero::Coordinate;

    use super::*;
    use crate::Reader;
    use crate::liftover::machine;

    #[test]
    pub fn test_valid_positive_strand_liftover() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let to = Coordinate::try_new("seq0", Strand::Positive, 7)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().inner(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Positive);
        assert_eq!(result.reference().start().position(), &Position::from(1));
        assert_eq!(result.reference().end().strand(), &Strand::Positive);
        assert_eq!(result.reference().end().position(), &Position::from(7));

        assert_eq!(result.query().contig().inner(), "seq1");
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
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", Strand::Negative, 7)?;
        let to = Coordinate::try_new("seq0", Strand::Negative, 1)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().inner(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Negative);
        assert_eq!(result.reference().start().position(), &Position::from(7));
        assert_eq!(result.reference().end().strand(), &Strand::Negative);
        assert_eq!(result.reference().end().position(), &Position::from(1));

        assert_eq!(result.query().contig().inner(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Negative);
        assert_eq!(result.query().start().position(), &Position::from(7));
        assert_eq!(result.query().end().strand(), &Strand::Negative);
        assert_eq!(result.query().end().position(), &Position::from(1));

        Ok(())
    }

    #[test]
    pub fn test_valid_positive_to_negative_strand_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let to = Coordinate::try_new("seq0", Strand::Positive, 7)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().inner(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Positive);
        assert_eq!(result.reference().start().position(), &Position::from(1));
        assert_eq!(result.reference().end().strand(), &Strand::Positive);
        assert_eq!(result.reference().end().position(), &Position::from(7));

        assert_eq!(result.query().contig().inner(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Negative);
        assert_eq!(result.query().start().position(), &Position::from(8));
        assert_eq!(result.query().end().strand(), &Strand::Negative);
        assert_eq!(result.query().end().position(), &Position::from(2));

        Ok(())
    }

    #[test]
    pub fn test_valid_negative_to_positive_strand_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", Strand::Negative, 7)?;
        let to = Coordinate::try_new("seq0", Strand::Negative, 1)?;

        let interval = Interval::try_new(from, to)?;
        let mut results = machine.liftover(&interval).unwrap();

        assert_eq!(results.len(), 1);

        let result = results.pop().unwrap();

        assert_eq!(result.reference().contig().inner(), "seq0");
        assert_eq!(result.reference().start().strand(), &Strand::Negative);
        assert_eq!(result.reference().start().position(), &Position::from(7));
        assert_eq!(result.reference().end().strand(), &Strand::Negative);
        assert_eq!(result.reference().end().position(), &Position::from(1));

        assert_eq!(result.query().contig().inner(), "seq1");
        assert_eq!(result.query().start().strand(), &Strand::Positive);
        assert_eq!(result.query().start().position(), &Position::from(2));
        assert_eq!(result.query().end().strand(), &Strand::Positive);
        assert_eq!(result.query().end().position(), &Position::from(8));

        Ok(())
    }

    #[test]
    pub fn test_nonexistent_positive_strand_interval_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 - 0 10 seq1 10 - 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", Strand::Positive, 1)?;
        let to = Coordinate::try_new("seq0", Strand::Positive, 7)?;

        let interval = Interval::try_new(from, to)?;
        let results = machine.liftover(&interval);

        assert_eq!(results, None);

        Ok(())
    }

    #[test]
    pub fn test_nonexistent_negative_strand_interval_liftover()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::try_new("seq0", Strand::Negative, 7)?;
        let to = Coordinate::try_new("seq0", Strand::Negative, 1)?;

        let interval = Interval::try_new(from, to)?;
        let results = machine.liftover(&interval);

        assert_eq!(results, None);

        Ok(())
    }
}
