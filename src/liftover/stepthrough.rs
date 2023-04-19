//! Utilities for stepping through a contiguous liftover segment.
//!
//! A contiguous liftover segment is defined as a pair of intervals (one from
//! the reference and one from the query) that have the same distance and match
//! to each other contiguously.

use std::slice::Iter;

use crate::core::coordinate;
use crate::core::interval;
use crate::core::Coordinate;
use crate::core::Interval;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
use crate::reader::AlignmentDataSection;
use crate::record::AlignmentDataRecord;

pub mod interval_pair;

/// An error related to stepping through liftover segments.
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// An error occurred when constructing a [`Coordinate`].
    InvalidCoordinate(coordinate::Error),
    /// An error occurred when constructing an [`Interval`].
    InvalidInterval(interval::Error),
    /// An error occurred when constructing a [`ContiguousIntervalPair`].
    InvalidIntervalPair(interval_pair::Error),
    /// The data section does not cumulative add up the specified end
    /// coordinates. This indicates a malformed data section.
    MisalignedDataSection,
    /// When stepping through the interval step through, moving the coordinates
    /// returned `None`, which indicates there was an out of bounds somewhere.
    IntervalStepthroughOutOfBounds(String, Coordinate, usize),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidCoordinate(err) => write!(f, "invalid coordinate: {}", err),
            Error::InvalidInterval(err) => write!(f, "invalid interval: {}", err),
            Error::InvalidIntervalPair(err) => {
                write!(f, "invalid interval pair: {}", err)
            }
            Error::MisalignedDataSection => write!(
                f,
                "misaligned data section, which indicates a malformed chain file"
            ),
            Error::IntervalStepthroughOutOfBounds(name, coordinate, magnitude) => {
                write!(
                    f,
                    "interval stepthrough out of bounds: moving {} from {} by {}",
                    name, coordinate, magnitude
                )
            }
        }
    }
}

impl std::error::Error for Error {}

/// The core struct using for stepping through contiguous liftover segments.
pub struct StepThrough<'a> {
    query_pointer: Coordinate,
    reference_pointer: Coordinate,
    expected_query_end: Coordinate,
    expected_reference_end: Coordinate,
    data_records: Iter<'a, AlignmentDataRecord>,
}

impl<'a> StepThrough<'a> {
    pub(crate) fn new(section: &'a AlignmentDataSection) -> Result<Self, Error> {
        let query_pointer = section
            .query_start_coordinate()
            .map_err(Error::InvalidCoordinate)?;

        let reference_pointer = section
            .reference_start_coordinate()
            .map_err(Error::InvalidCoordinate)?;

        let expected_query_end = section
            .query_end_coordinate()
            .map_err(Error::InvalidCoordinate)?;
        let expected_reference_end = section
            .reference_end_coordinate()
            .map_err(Error::InvalidCoordinate)?;

        Ok(Self {
            query_pointer,
            reference_pointer,
            expected_query_end,
            expected_reference_end,
            data_records: section.alignment_data_records().iter(),
        })
    }

    fn finish(&self) -> Result<(), Error> {
        if self.reference_pointer != self.expected_reference_end {
            return Err(Error::MisalignedDataSection);
        }

        if self.query_pointer != self.expected_query_end {
            return Err(Error::MisalignedDataSection);
        }

        Ok(())
    }
}

impl<'a> Iterator for StepThrough<'a> {
    type Item = Result<ContiguousIntervalPair, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let chunk = match self.data_records.next() {
            Some(c) => c,
            None => match self.finish() {
                Ok(_) => return None,
                Err(e) => return Some(Err(e)),
            },
        };

        // (1) Calculate the coordinates for the reference interval.
        // (a) Set the reference start to point to the current reference pointer.
        let reference_start = self.reference_pointer.clone();

        // (b) Push forward the reference pointer by the amount specified in the chunk.
        self.reference_pointer = match self.reference_pointer.clone().move_forward(chunk.size()) {
            Some(coordinate) => match coordinate {
                Ok(coordinate) => coordinate,
                Err(err) => return Some(Err(Error::InvalidCoordinate(err))),
            },
            None => {
                return Some(Err(Error::IntervalStepthroughOutOfBounds(
                    String::from("reference forward by chunk size"),
                    self.reference_pointer.clone(),
                    chunk.size(),
                )))
            }
        };

        // (c) Set the reference end to point to the current reference pointer.
        let reference_end = self.reference_pointer.clone();

        // (d) Create the reference interval.
        let reference = match Interval::try_new(reference_start, reference_end) {
            Ok(interval) => interval,
            Err(err) => return Some(Err(Error::InvalidInterval(err))),
        };

        // (2) Calculate the coordinates for the query interval.
        // (a) Set the query start to point to the current query pointer.
        let query_start = self.query_pointer.clone();

        // (b) Push forward the query pointer by the amount specified in the chunk.
        self.query_pointer = match self.query_pointer.clone().move_forward(chunk.size()) {
            Some(coordinate) => match coordinate {
                Ok(coordinate) => coordinate,
                Err(err) => return Some(Err(Error::InvalidCoordinate(err))),
            },
            None => {
                return Some(Err(Error::IntervalStepthroughOutOfBounds(
                    String::from("query forward by chunk size"),
                    self.query_pointer.clone(),
                    chunk.size(),
                )))
            }
        };

        // (c) Set the query end to point to the current query pointer.
        let query_end = self.query_pointer.clone();

        // (d) Create the query interval.
        let query = match Interval::try_new(query_start, query_end) {
            Ok(interval) => interval,
            Err(err) => return Some(Err(Error::InvalidInterval(err))),
        };

        // (3) Adjust the pointers to the current query and reference locations.
        // (a) If there is a query offset, add it to the query pointer.
        if let Some(dq) = chunk.dq() {
            self.query_pointer = match self.query_pointer.clone().move_forward(*dq) {
                Some(coordinate) => match coordinate {
                    Ok(coordinate) => coordinate,
                    Err(err) => return Some(Err(Error::InvalidCoordinate(err))),
                },
                None => {
                    return Some(Err(Error::IntervalStepthroughOutOfBounds(
                        String::from("query forward by dq"),
                        self.query_pointer.clone(),
                        *dq,
                    )))
                }
            }
        }

        // (b) If there is a reference offset, add it to the reference pointer.
        if let Some(dt) = chunk.dt() {
            self.reference_pointer = match self.reference_pointer.clone().move_forward(*dt) {
                Some(coordinate) => match coordinate {
                    Ok(coordinate) => coordinate,
                    Err(err) => return Some(Err(Error::InvalidCoordinate(err))),
                },
                None => {
                    return Some(Err(Error::IntervalStepthroughOutOfBounds(
                        String::from("reference forward by dt"),
                        self.reference_pointer.clone(),
                        *dt,
                    )))
                }
            }
        }

        // (4) Create the final contiguous interval pair.
        let pair = match ContiguousIntervalPair::try_new(reference, query) {
            Ok(p) => p,
            Err(err) => return Some(Err(Error::InvalidIntervalPair(err))),
        };

        Some(Ok(pair))
    }
}

#[cfg(test)]
mod tests {
    use crate::core::Coordinate;
    use crate::core::Interval;
    use crate::core::Strand;
    use crate::liftover::stepthrough::Error;

    #[test]
    fn test_it_correctly_steps_through_a_valid_chain_file() -> Result<(), Box<dyn std::error::Error>>
    {
        let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
        let mut reader = crate::Reader::new(&data[..]);

        let mut results = reader.sections().map(|x| x.unwrap()).collect::<Vec<_>>();
        assert_eq!(results.len(), 1);

        // SAFETY: we just checked that the length was 1.
        let section = results.pop().unwrap();
        let mut stepthrough = section.stepthrough()?;

        //===============================================//
        // Checks the first interval in the step through //
        //===============================================//

        let reference = Interval::try_new(
            Coordinate::try_new("seq0", 0, Strand::Positive)?,
            Coordinate::try_new("seq0", 3, Strand::Positive)?,
        )?;
        let query = Interval::try_new(
            Coordinate::try_new("seq0", 4, Strand::Negative)?,
            Coordinate::try_new("seq0", 1, Strand::Negative)?,
        )?;

        let result = stepthrough.next().unwrap().unwrap();
        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //================================================//
        // Checks the second interval in the step through //
        //================================================//

        let reference = Interval::try_new(
            Coordinate::try_new("seq0", 3, Strand::Positive)?,
            Coordinate::try_new("seq0", 4, Strand::Positive)?,
        )?;
        let query = Interval::try_new(
            Coordinate::try_new("seq0", 0, Strand::Negative)?,
            Coordinate::negative_bound("seq0"),
        )?;

        let result = stepthrough.next().unwrap().unwrap();
        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //==========================================//
        // Checks that the iterator is finished now //
        //==========================================//

        assert_eq!(stepthrough.next(), None);

        Ok(())
    }

    #[test]
    fn test_it_errors_when_there_is_an_alignment_issue() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 4 + 0 3 seq0 10 - 0 9 1\n2\t0\t1\n1";
        let mut reader = crate::Reader::new(&data[..]);

        let mut results = reader.sections().map(|x| x.unwrap()).collect::<Vec<_>>();
        assert_eq!(results.len(), 1);

        // SAFETY: we just checked that the length was 1.
        let section = results.pop().unwrap();
        let mut stepthrough = section.stepthrough()?;

        //===============================================//
        // Checks the first interval in the step through //
        //===============================================//

        let reference = Interval::try_new(
            Coordinate::try_new("seq0", 0, Strand::Positive)?,
            Coordinate::try_new("seq0", 2, Strand::Positive)?,
        )?;
        let query = Interval::try_new(
            Coordinate::try_new("seq0", 9, Strand::Negative)?,
            Coordinate::try_new("seq0", 7, Strand::Negative)?,
        )?;

        let result = stepthrough.next().unwrap().unwrap();
        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //================================================//
        // Checks the second interval in the step through //
        //================================================//

        let reference = Interval::try_new(
            Coordinate::try_new("seq0", 2, Strand::Positive)?,
            Coordinate::try_new("seq0", 3, Strand::Positive)?,
        )?;
        let query = Interval::try_new(
            Coordinate::try_new("seq0", 6, Strand::Negative)?,
            Coordinate::try_new("seq0", 5, Strand::Negative)?,
        )?;

        let result = stepthrough.next().unwrap().unwrap();
        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //============================================================//
        // Checks that we now receive a misaligned data section error //
        //============================================================//

        assert_eq!(stepthrough.next(), Some(Err(Error::MisalignedDataSection)));

        Ok(())
    }
}
