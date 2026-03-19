//! A builder for a [`Machine`].

use std::collections::HashMap;
use std::io::BufRead;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::position::Number;
use rust_lapper as lapper;

use thiserror::Error;

use crate::alignment;
use crate::alignment::section::header::Record as Header;
use crate::liftover;
use super::AnnotatedPair;
use super::ChromosomeDictionary;
use super::Machine;
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
        "duplicate chain id `{id}` with inconsistent headers: expected `{expected}`, found `{found}`"
    )]
    DuplicateChainId {
        /// The chain ID.
        id: usize,

        /// The first header encountered for this chain ID.
        expected: Header,

        /// The conflicting header encountered for this chain ID.
        found: Header,
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

        let mut reference_chromosomes = ChromosomeDictionaryBuilder::default();
        let mut query_chromosomes = ChromosomeDictionaryBuilder::default();

        for result in reader.sections() {
            let section = result.map_err(Error::InvalidSections)?;

            let header = section.header().clone();

            match chains.get(&header.id()) {
                Some(existing) if *existing != header => {
                    return Err(Error::DuplicateChainId {
                        id: header.id(),
                        expected: existing.clone(),
                        found: header,
                    });
                }
                _ => {
                    chains.insert(header.id(), header.clone());
                }
            }

            query_chromosomes.update(
                header.query_sequence().chromosome_name().to_string(),
                header.query_sequence().chromosome_size(),
            );
            reference_chromosomes.update(
                header.reference_sequence().chromosome_name().to_string(),
                header.reference_sequence().chromosome_size(),
            );

            for pair_result in section.stepthrough().map_err(Error::StepthroughError)? {
                let pair = pair_result.map_err(Error::StepthroughError)?;
                let entry = hm.entry(pair.reference().contig().clone()).or_default();

                let (start, stop) = match pair.reference().strand() {
                    Strand::Positive => {
                        let start = pair.reference().start().position().get();
                        let end = pair.reference().end().position().get();

                        (start, end)
                    }
                    Strand::Negative => {
                        let start = pair.reference().end().position().get();
                        let end = pair.reference().start().position().get();

                        (start, end)
                    }
                };

                entry.push(lapper::Interval {
                    start,
                    stop,
                    val: AnnotatedPair {
                        chain: header.clone(),
                        pair,
                    },
                })
            }
        }

        let mut inner = HashMap::<Contig, lapper::Lapper<Number, AnnotatedPair>>::new();

        for (k, v) in hm.into_iter() {
            inner.insert(k, lapper::Lapper::new(v));
        }

        let reference_chromosomes = reference_chromosomes.build();
        let query_chromosomes = query_chromosomes.build();

        Ok(Machine {
            inner,
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

/// A builder for the size of the contigs in either the reference or query
/// sequences.
#[derive(Default)]
struct ChromosomeDictionaryBuilder(ChromosomeDictionary);

impl ChromosomeDictionaryBuilder {
    /// Updates the contig builder with the current start and end.
    fn update(&mut self, contig: String, size: Number) {
        match self.0.get(&contig) {
            Some(existing) => {
                assert!(
                    *existing == size,
                    "the current size conflicts with the previously set size"
                );
            }
            None => {
                self.0.insert(contig, size);
            }
        }
    }

    /// Consumes `self` and returns the built [`ChromosomeDictionary`].
    fn build(self) -> ChromosomeDictionary {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use crate::Reader;

    use super::*;

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
}
