//! A machine for lifting over coordinates within a reference genome.

use std::collections::HashMap;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interval::interbase::Interval;
use omics::coordinate::position::Number;
use rust_lapper as lapper;

use crate::alignment::section::header::Record as Header;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;

pub mod builder;

pub use builder::Builder;

/// A dictionary of chromosome names and their size.
pub type ChromosomeDictionary = HashMap<String, Number>;

/// A mapping segment annotated with the chain from which it originated.
#[derive(Clone, Debug, Eq, PartialEq)]
struct AnnotatedPair {
    /// The chain from which this segment originated.
    chain: Header,

    /// The mapping segment.
    pair: ContiguousIntervalPair,
}

/// The liftover results from a single chain.
///
/// Each chain produces one or more contiguous mapping segments. When a
/// query interval spans a gap within a chain, the chain contributes
/// multiple segments (one per aligned block that overlaps the query).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LiftoverResult {
    /// The chain from which these liftover results originated.
    chain: Header,

    /// The mapping segments from this chain that overlap the query.
    segments: Vec<ContiguousIntervalPair>,
}

impl LiftoverResult {
    /// Gets the chain from which these liftover results originated.
    pub fn chain(&self) -> &Header {
        &self.chain
    }

    /// Gets the mapping segments.
    pub fn segments(&self) -> &[ContiguousIntervalPair] {
        &self.segments
    }

    /// Consumes `self` and returns the mapping segments.
    pub fn into_segments(self) -> Vec<ContiguousIntervalPair> {
        self.segments
    }
}

/// A machine for lifting over coordinates from a reference genome to a query
/// genome.
///
/// Generally, you will want to use a [`builder::Builder`] to construct one of
/// these.
#[derive(Debug)]
pub struct Machine {
    /// The inner lookup table of positions in the reference genome to positions
    /// in the query genome for each contig in the reference genome.
    inner: HashMap<Contig, lapper::Lapper<Number, AnnotatedPair>>,

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

    /// Performs a liftover from the specified `interval` to the query genome.
    ///
    /// Results are grouped by the chain from which each mapping segment
    /// originated. Each [`LiftoverResult`] contains one or more contiguous
    /// segments from a single chain, along with the chain's header.
    pub fn liftover(&self, interval: Interval) -> Option<Vec<LiftoverResult>> {
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

        let clamped = entry
            .find(start, stop)
            .map(|e| e.val.clone())
            .filter(|ap| ap.pair.reference().strand() == strand)
            .map(|ap| {
                // SAFETY: the lapper guarantees that `ap.pair` overlaps
                // `interval`, so `clamp()` will always succeed.
                let pair = ap.pair.clamp(interval.clone()).unwrap();
                (pair, ap.chain)
            })
            .collect::<Vec<_>>();

        if clamped.is_empty() {
            return None;
        }

        let mut groups = HashMap::<usize, (Header, Vec<ContiguousIntervalPair>)>::new();

        for (pair, chain) in clamped {
            groups
                .entry(chain.id())
                .or_insert_with(|| (chain, Vec::new()))
                .1
                .push(pair);
        }

        let results = groups
            .into_values()
            .map(|(chain, segments)| LiftoverResult { chain, segments })
            .collect::<Vec<_>>();

        Some(results)
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
        let results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 1);

        let result = &results[0].segments()[0];

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
        let results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 1);

        let result = &results[0].segments()[0];

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
        let results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 1);

        let result = &results[0].segments()[0];

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
        let results = machine.liftover(interval).unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 1);

        let result = &results[0].segments()[0];

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

    #[test]
    fn test_grouped_by_chain() -> Result<(), Box<dyn std::error::Error>> {
        // Two chains covering the same region on `seq0`, mapping to
        // `seq1` and `seq2` respectively.
        let data =
            b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 50 seq0 10 + 0 10 seq2 10 + 0 10 1\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Positive, 1u64);
        let to = Coordinate::new("seq0", Strand::Positive, 5u64);
        let interval = Interval::try_new(from, to)?;

        let mut results = machine.liftover(interval).unwrap();
        assert_eq!(results.len(), 2);

        // Sort by chain ID for deterministic assertions.
        results.sort_by_key(|r| r.chain().id());

        assert_eq!(results[0].chain().id(), 0);
        assert_eq!(results[0].chain().score(), 100);
        assert_eq!(results[0].segments().len(), 1);
        assert_eq!(results[0].segments()[0].query().contig().as_str(), "seq1");

        assert_eq!(results[1].chain().id(), 1);
        assert_eq!(results[1].chain().score(), 50);
        assert_eq!(results[1].segments().len(), 1);
        assert_eq!(results[1].segments()[0].query().contig().as_str(), "seq2");

        Ok(())
    }

    #[test]
    fn test_gapped_chain_returns_multiple_segments() -> Result<(), Box<dyn std::error::Error>> {
        // One chain with a gap: two blocks of 4, separated by a gap of
        // 2 in both reference and query.
        let data = b"chain 0 seq0 10 + 0 10 seq1 10 + 0 10 0\n4\t2\t2\n4";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new("seq0", Strand::Positive, 0u64);
        let to = Coordinate::new("seq0", Strand::Positive, 10u64);
        let interval = Interval::try_new(from, to)?;

        let results = machine.liftover(interval).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 2);
        assert_eq!(results[0].segments()[0].reference().count_entities(), 4);
        assert_eq!(results[0].segments()[1].reference().count_entities(), 4);

        Ok(())
    }
}
