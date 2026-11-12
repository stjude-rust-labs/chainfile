//! A builder for a [`Machine`].

use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::io::BufRead;
use std::ops::Range;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::position::Number;
use rust_lapper as lapper;
use thiserror::Error;

use super::AnnotatedPair;
use super::ChromosomeDictionary;
use super::Machine;
use crate::alignment;
use crate::alignment::section::header::Record as Header;
use crate::liftover;
use crate::reader;

/// The inner value of the liftover lookup data structure.
type Iv = lapper::Interval<Number, AnnotatedPair>;

/// An error related to building a [`Machine`].
#[derive(Debug, Error)]
pub enum Error {
    /// An error reading alignment sections.
    #[error("invalid data section: {0}")]
    InvalidSections(alignment::section::sections::Error),

    /// An error stepping through the liftover segments.
    #[error("stepthrough error: {0}")]
    StepthroughError(liftover::stepthrough::Error),

    /// Two chain sections share the same ID but have different headers.
    #[error(
        "duplicate chain id `{id}` with inconsistent headers: expected `{expected}`, found \
         `{found}`"
    )]
    DuplicateChainId {
        /// The chain ID.
        id: usize,

        /// The first header encountered for this chain ID.
        expected: Header,

        /// The conflicting header encountered for this chain ID.
        found: Header,
    },

    /// Two chain sections share the same ID and header but have different
    /// alignment data.
    #[error("duplicate chain id `{id}` with inconsistent alignment data")]
    InconsistentChainData {
        /// The chain ID.
        id: usize,
    },

    /// Two chain sections report different sizes for the same chromosome.
    #[error("conflicting sizes for chromosome `{contig}`: expected `{expected}`, found `{found}`")]
    ConflictingChromosomeSize {
        /// The chromosome name.
        contig: String,

        /// The first chromosome size encountered.
        expected: Number,

        /// The conflicting chromosome size.
        found: Number,
    },
}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

/// A builder for a [`Machine`].
#[allow(missing_debug_implementations)]
pub struct Builder;

impl Builder {
    /// Builds a [`Machine`] from the builder.
    ///
    /// # Examples
    ///
    /// ```
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let reader = chainfile::Reader::new(&data[..]);
    ///
    /// let result = chainfile::liftover::machine::Builder::default().try_build_from(reader)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[allow(clippy::result_large_err)]
    pub fn try_build_from<T>(&self, mut reader: reader::Reader<T>) -> Result<Machine>
    where
        T: BufRead,
    {
        let mut hm = HashMap::<Contig, Vec<Iv>>::default();
        let mut chains = HashMap::<usize, Header>::new();

        // The segments for each chain are pushed contiguously, so each chain
        // ID maps to the range of its segments within its contig's entries.
        let mut segment_ranges = HashMap::<usize, Range<usize>>::new();

        let mut reference_chromosomes = ChromosomeDictionaryBuilder::default();
        let mut query_chromosomes = ChromosomeDictionaryBuilder::default();

        for result in reader.sections() {
            let section = result.map_err(Error::InvalidSections)?;

            let header = section.header().clone();
            let contig = Contig::new_unchecked(header.reference_sequence().chromosome_name());

            match chains.entry(header.id()) {
                Entry::Occupied(e) => {
                    if *e.get() != header {
                        return Err(Error::DuplicateChainId {
                            id: header.id(),
                            expected: e.get().clone(),
                            found: header,
                        });
                    }

                    // Since the headers match, the segments match if and only
                    // if the alignment data match.
                    let existing = &hm[&contig][segment_ranges[&header.id()].clone()];
                    let mut duplicate = Vec::with_capacity(existing.len());
                    push_segments(&section, &mut duplicate)?;

                    if !existing
                        .iter()
                        .map(|iv| &iv.val)
                        .eq(duplicate.iter().map(|iv| &iv.val))
                    {
                        return Err(Error::InconsistentChainData { id: header.id() });
                    }

                    // The section is an exact duplicate of one that was
                    // already added, so adding its segments again would
                    // duplicate them in the liftover results.
                    continue;
                }
                Entry::Vacant(e) => {
                    e.insert(header.clone());
                }
            }

            query_chromosomes.update(
                header.query_sequence().chromosome_name().to_string(),
                header.query_sequence().chromosome_size(),
            )?;
            reference_chromosomes.update(
                header.reference_sequence().chromosome_name().to_string(),
                header.reference_sequence().chromosome_size(),
            )?;

            let entry = hm.entry(contig).or_default();
            let start = entry.len();
            push_segments(&section, entry)?;
            segment_ranges.insert(header.id(), start..entry.len());
        }

        let mut inner = HashMap::<Contig, lapper::Lapper<Number, AnnotatedPair>>::new();

        for (k, v) in hm.into_iter() {
            inner.insert(k, lapper::Lapper::new(v));
        }

        let reference_chromosomes = reference_chromosomes.build();
        let query_chromosomes = query_chromosomes.build();

        Ok(Machine {
            inner,
            chains,
            reference_chromosomes,
            query_chromosomes,
        })
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self
    }
}

/// Steps through a section and pushes its liftover segments onto `entries`.
#[expect(
    clippy::result_large_err,
    reason = "`Error` retains the conflicting chain headers"
)]
fn push_segments(section: &alignment::Section, entries: &mut Vec<Iv>) -> Result<()> {
    let chain_id = section.header().id();

    for pair_result in section.stepthrough().map_err(Error::StepthroughError)? {
        let pair = pair_result.map_err(Error::StepthroughError)?;

        let (start, stop) = match pair.reference().strand() {
            Strand::Positive => (
                pair.reference().start().position().get(),
                pair.reference().end().position().get(),
            ),
            Strand::Negative => (
                pair.reference().end().position().get(),
                pair.reference().start().position().get(),
            ),
        };

        entries.push(lapper::Interval {
            start,
            stop,
            val: AnnotatedPair { chain_id, pair },
        });
    }

    Ok(())
}

/// A builder for the size of the contigs in either the reference or query
/// sequences.
#[derive(Default)]
struct ChromosomeDictionaryBuilder(ChromosomeDictionary);

impl ChromosomeDictionaryBuilder {
    /// Updates the chromosome dictionary with a contig and its size.
    ///
    /// # Errors
    ///
    /// Returns [`Error::ConflictingChromosomeSize`] when the dictionary already
    /// contains the contig with a different size.
    #[expect(
        clippy::result_large_err,
        reason = "`Error` retains the conflicting chain headers"
    )]
    fn update(&mut self, contig: String, size: Number) -> Result<()> {
        match self.0.get(&contig) {
            Some(existing) => {
                if *existing != size {
                    return Err(Error::ConflictingChromosomeSize {
                        contig,
                        expected: *existing,
                        found: size,
                    });
                }
            }
            None => {
                self.0.insert(contig, size);
            }
        }

        Ok(())
    }

    /// Consumes `self` and returns the built [`ChromosomeDictionary`].
    fn build(self) -> ChromosomeDictionary {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use omics::coordinate::interval::interbase::Interval;

    use super::*;
    use crate::Reader;

    #[test]
    fn duplicate_chain_id_with_inconsistent_headers() {
        // Two chains with the same ID (`0`) but different scores.
        let data =
            b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 50 seq0 10 + 0 10 seq2 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        let err = Builder.try_build_from(reader).unwrap_err();

        assert!(
            matches!(err, Error::DuplicateChainId { id: 0, .. }),
            "expected DuplicateChainId error, got: {err}"
        );
    }

    #[test]
    fn duplicate_chain_id_with_consistent_headers_succeeds() {
        // Two identical chain sections with the same ID—this is fine.
        let data =
            b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10";
        let reader = Reader::new(&data[..]);
        assert!(Builder.try_build_from(reader).is_ok());
    }

    #[test]
    fn exact_duplicate_chains_are_not_double_counted() {
        let data = b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n4\t2\t2\n4\n\n\
                     chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n4\t2\t2\n4";
        let machine = Builder.try_build_from(Reader::new(&data[..])).unwrap();

        let results = machine
            .liftover("seq0:+:0-10".parse::<Interval>().unwrap())
            .unwrap();

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].segments().len(), 2);
    }

    #[test]
    fn identical_chains_with_different_ids_are_kept() {
        let data = b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\n\
                     chain 100 seq0 10 + 0 10 seq1 10 + 0 10 1\n10";
        let machine = Builder.try_build_from(Reader::new(&data[..])).unwrap();

        let results = machine
            .liftover("seq0:+:0-10".parse::<Interval>().unwrap())
            .unwrap();

        assert_eq!(results.len(), 2);
    }

    #[test]
    fn duplicate_chain_id_with_inconsistent_data() {
        let data = b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n3\t1\t1\n6\n\n\
                     chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n4\t1\t1\n5";
        let err = Builder.try_build_from(Reader::new(&data[..])).unwrap_err();

        assert!(matches!(err, Error::InconsistentChainData { id: 0 }));
        assert_eq!(
            err.to_string(),
            "duplicate chain id `0` with inconsistent alignment data"
        );
    }

    #[test]
    fn conflicting_reference_chromosome_sizes() {
        let data = b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 50 seq0 11 + 0 11 seq2 11 + 0 11 1\n11";
        let reader = Reader::new(&data[..]);
        let err = Builder.try_build_from(reader).unwrap_err();

        assert!(matches!(
            err,
            Error::ConflictingChromosomeSize {
                ref contig,
                expected: 10,
                found: 11,
            } if contig == "seq0"
        ));
        assert_eq!(
            err.to_string(),
            "conflicting sizes for chromosome `seq0`: expected `10`, found `11`"
        );
    }

    #[test]
    fn conflicting_query_chromosome_sizes() {
        let data = b"chain 100 seq0 10 + 0 10 seq1 10 + 0 10 0\n10\n\nchain 50 seq2 11 + 0 11 seq1 11 + 0 11 1\n11";
        let reader = Reader::new(&data[..]);
        let err = Builder.try_build_from(reader).unwrap_err();

        assert!(matches!(
            err,
            Error::ConflictingChromosomeSize {
                ref contig,
                expected: 10,
                found: 11,
            } if contig == "seq1"
        ));
    }
}
