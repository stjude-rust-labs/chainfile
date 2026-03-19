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

/// A mapping segment annotated with the chain ID from which it originated.
#[derive(Clone, Debug, Eq, PartialEq)]
struct AnnotatedPair {
    /// The chain ID from which this segment originated.
    chain_id: usize,

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

    /// A lookup table from chain ID to the chain header.
    chains: HashMap<usize, Header>,

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
    ///
    /// The returned `Vec<LiftoverResult>` is sorted by score descending
    /// (best chain first), with chain ID ascending as a tiebreaker. Within
    /// each result, `segments()` are sorted by reference interval start
    /// position ascending, with end position as a tiebreaker.
    pub fn liftover(&self, interval: Interval) -> Option<Vec<LiftoverResult>> {
        let entry = self.inner.get(interval.contig())?;

        let (start, stop) = match interval.strand() {
            Strand::Positive => (
                interval.start().position().get(),
                interval.end().position().get(),
            ),
            Strand::Negative => (
                interval.end().position().get(),
                interval.start().position().get(),
            ),
        };

        let strand = interval.strand();

        let mut hits = Vec::<(usize, ContiguousIntervalPair)>::new();

        for e in entry.find(start, stop) {
            if e.val.pair.reference().strand() != strand {
                continue;
            }

            // SAFETY: the lapper guarantees that `e.val.pair` overlaps
            // `interval`, so `clamp()` will always succeed.
            let pair = e.val.pair.clone().clamp(interval.clone()).unwrap();
            hits.push((e.val.chain_id, pair));
        }

        if hits.is_empty() {
            return None;
        }

        // Sort by chain ID so that consecutive entries with the same ID are
        // adjacent, then group them.
        hits.sort_by_key(|(id, _)| *id);

        let mut results = Vec::<LiftoverResult>::new();

        for (chain_id, pair) in hits {
            match results.last_mut() {
                Some(last) if last.chain().id() == chain_id => {
                    last.segments.push(pair);
                }
                _ => {
                    let chain = self.chains[&chain_id].clone();
                    results.push(LiftoverResult {
                        chain,
                        segments: vec![pair],
                    });
                }
            }
        }

        // Within each chain group, ensure segments are sorted by reference
        // start position.
        for result in &mut results {
            result.segments.sort_by(|a, b| {
                a.reference()
                    .start()
                    .cmp(b.reference().start())
                    .then_with(|| a.reference().end().cmp(b.reference().end()))
            });
        }

        // Sort results by score descending so the best chain comes first,
        // with chain ID ascending as a tiebreaker.
        results.sort_by(|a, b| {
            b.chain()
                .score()
                .cmp(&a.chain().score())
                .then_with(|| a.chain().id().cmp(&b.chain().id()))
        });

        Some(results)
    }
}

#[cfg(test)]
mod tests {
    use omics::coordinate::Contig;
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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 1u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 7u64);

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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 7u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 1u64);

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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 1u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 7u64);

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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 7u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 1u64);

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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 1u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 7u64);

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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 7u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 1u64);

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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 1u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 5u64);
        let interval = Interval::try_new(from, to)?;

        let results = machine.liftover(interval).unwrap();
        assert_eq!(results.len(), 2);

        // Results are sorted by score descending.
        assert_eq!(results[0].chain().score(), 100);
        assert_eq!(results[0].chain().id(), 0);
        assert_eq!(results[0].segments().len(), 1);
        assert_eq!(results[0].segments()[0].query().contig().as_str(), "seq1");

        assert_eq!(results[1].chain().score(), 50);
        assert_eq!(results[1].chain().id(), 1);
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

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 0u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 10u64);
        let interval = Interval::try_new(from, to)?;

        let results = machine.liftover(interval).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 2);
        assert_eq!(results[0].segments()[0].reference().count_entities(), 4);
        assert_eq!(results[0].segments()[1].reference().count_entities(), 4);

        Ok(())
    }

    #[test]
    fn test_liftover_ordering_is_deterministic() -> Result<(), Box<dyn std::error::Error>> {
        // Three chains with IDs 2, 0, 1 (intentionally out of order in the
        // input) covering the same region on `seq0`.
        let data = b"chain 30 seq0 10 + 0 10 seq3 10 + 0 10 2\n10\n\n\
                      chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\n\
                      chain 50 seq0 10 + 0 10 seq2 10 + 0 10 1\n10";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 0u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 10u64);
        let interval = Interval::try_new(from, to)?;

        // Run liftover multiple times and verify ordering is stable.
        for _ in 0..10 {
            let results = machine.liftover(interval.clone()).unwrap();
            assert_eq!(results.len(), 3);

            // Results must be sorted by score descending.
            assert_eq!(results[0].chain().score(), 100);
            assert_eq!(results[1].chain().score(), 50);
            assert_eq!(results[2].chain().score(), 30);
        }

        Ok(())
    }

    #[test]
    fn test_segments_ordered_by_reference_position() -> Result<(), Box<dyn std::error::Error>> {
        // One chain with three blocks separated by gaps, ensuring segments
        // come back in reference-position order. Layout: 5 + gap(2,2) + 5 +
        // gap(2,2) + 6 = 20 reference bases.
        let data = b"chain 0 seq0 20 + 0 20 seq1 20 + 0 20 0\n5\t2\t2\n5\t2\t2\n6";
        let reader = Reader::new(&data[..]);
        let machine = machine::Builder.try_build_from(reader)?;

        let from = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 0u64);
        let to = Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 20u64);
        let interval = Interval::try_new(from, to)?;

        let results = machine.liftover(interval).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 3);

        // Segments must be sorted by reference start position.
        let starts: Vec<_> = results[0]
            .segments()
            .iter()
            .map(|s| s.reference().start().position().get())
            .collect();
        assert_eq!(starts, vec![0, 7, 14]);

        Ok(())
    }
}
