//! A builder for a [`Machine`].

use std::collections::HashMap;
use std::io::BufRead;

use rust_lapper as lapper;

use crate::core::coordinate::Position;
use crate::core::Contig;
use crate::core::Strand;
use crate::liftover;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
use crate::liftover::Machine;
use crate::reader;

type Iv = lapper::Interval<usize, ContiguousIntervalPair>;

/// An error related to building a [`Machine`].
#[derive(Debug)]
pub enum Error {
    /// An error related to a malformed data section.
    InvalidDataSection(Box<reader::section::ParseError>),
    /// An error related to stepping through the liftover segments.
    StepthroughError(liftover::stepthrough::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidDataSection(err) => write!(f, "invalid data section: {}", err),
            Error::StepthroughError(err) => write!(f, "stepthrough error: {}", err),
        }
    }
}

impl std::error::Error for Error {}

/// A builder for a [`Machine`].
pub struct Builder;

impl Builder {
    /// Builds a [`Machine`] from the builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let reader = chain::Reader::new(&data[..]);
    ///
    /// let result = chain::liftover::machine::Builder::default().try_build_from(reader)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_build_from<T>(&self, mut reader: reader::Reader<T>) -> Result<Machine, Error>
    where
        T: BufRead,
    {
        let mut hm = HashMap::<Contig, Vec<Iv>>::default();

        for section_result in reader.sections() {
            let section = section_result.map_err(Error::InvalidDataSection)?;
            for pair_result in section.stepthrough().map_err(Error::StepthroughError)? {
                let pair = pair_result.map_err(Error::StepthroughError)?;
                let entry = hm.entry(pair.reference().contig().clone()).or_default();

                let (start, stop) = match pair.reference().strand() {
                    Strand::Positive => {
                        let start = match pair.reference().start().position() {
                            Position::ZeroBased(a) => *a,
                            Position::NegativeBound => {
                                unreachable!("negative bound not allowed on positive strand")
                            }
                        };

                        let stop = match pair.reference().end().position() {
                            Position::ZeroBased(b) => *b,
                            Position::NegativeBound => {
                                unreachable!("negative bound not allowed on positive strand")
                            }
                        };

                        (start, stop)
                    }
                    Strand::Negative => {
                        let start = match pair.reference().end().position() {
                            Position::ZeroBased(a) => a + 1,
                            Position::NegativeBound => 0,
                        };

                        let stop = match pair.reference().start().position() {
                            Position::ZeroBased(b) => b + 1,
                            Position::NegativeBound => unreachable!(
                                "negative bound will never be the start of a negative-stranded \
                                 interval"
                            ),
                        };

                        (start, stop)
                    }
                };

                entry.push(lapper::Interval {
                    start,
                    stop,
                    val: pair,
                })
            }
        }

        let mut inner = HashMap::<Contig, lapper::Lapper<usize, ContiguousIntervalPair>>::new();

        for (k, v) in hm.into_iter() {
            inner.insert(k, lapper::Lapper::new(v));
        }

        Ok(Machine { inner })
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self
    }
}
