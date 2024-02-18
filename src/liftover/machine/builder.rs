//! A builder for a [`Machine`].

use std::collections::HashMap;
use std::io::BufRead;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::position::Value;
use rust_lapper as lapper;

use crate::alignment;
use crate::liftover;
use crate::liftover::Machine;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
use crate::reader;

/// The inner value of the liftover lookup data structure.
type Iv = lapper::Interval<usize, ContiguousIntervalPair>;

/// An error related to building a [`Machine`].
#[derive(Debug)]
pub enum Error {
    /// An error reading alignment sections.
    InvalidSections(alignment::section::sections::Error),

    /// An error stepping through the liftover segments.
    StepthroughError(liftover::stepthrough::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidSections(err) => write!(f, "invalid data section: {}", err),
            Error::StepthroughError(err) => write!(f, "stepthrough error: {}", err),
        }
    }
}

impl std::error::Error for Error {}

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

        for section_result in reader.sections() {
            let section = section_result.map_err(Error::InvalidSections)?;
            for pair_result in section.stepthrough().map_err(Error::StepthroughError)? {
                let pair = pair_result.map_err(Error::StepthroughError)?;
                let entry = hm.entry(pair.reference().contig().clone()).or_default();

                let (start, end) = match pair.reference().strand() {
                    Strand::Positive => {
                        let start = match pair.reference().start().position().inner() {
                            Value::Usize(a) => *a,
                            Value::LowerBound => {
                                unreachable!("lower bound not allowed on positive strand")
                            }
                        };

                        let end = match pair.reference().end().position().inner() {
                            Value::Usize(b) => *b,
                            Value::LowerBound => {
                                unreachable!("lower bound not allowed on positive strand")
                            }
                        };

                        (start, end)
                    }
                    Strand::Negative => {
                        let start = match pair.reference().end().position().inner() {
                            Value::Usize(a) => a + 1,
                            Value::LowerBound => 0,
                        };

                        let end = match pair.reference().start().position().inner() {
                            Value::Usize(b) => b + 1,
                            Value::LowerBound => unreachable!(
                                "lower bound will never be the start of a negative-stranded \
                                 interval"
                            ),
                        };

                        (start, end)
                    }
                };

                entry.push(lapper::Interval {
                    start,
                    stop: end,
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
