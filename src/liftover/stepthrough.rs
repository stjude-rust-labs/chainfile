//! Facilities for stepping through a contiguous liftover segment.
//!
//! A contiguous liftover segment is defined as a pair of intervals (one from
//! the reference and one from the query) that have the same distance and match
//! to each other contiguously.

use omics::coordinate;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;
use omics::coordinate::position::Number;
use thiserror::Error;
use tracing_log::log::warn;

use crate::alignment::Section;
use crate::alignment::section::data;
use crate::alignment::section::data::Record;
use crate::alignment::section::header::sequence;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;

pub mod interval_pair;

/// An error related to stepping through liftover segments.
#[derive(Debug, Error)]
pub enum Error {
    /// When stepping through the interval step through, moving the coordinates
    /// returned `None`, which indicates there was an out of bounds somewhere.
    #[error("interval stepthrough out of bounds: moving {0} from {1} by {2}")]
    IntervalStepthroughOutOfBounds(String, Coordinate, Number),

    /// An interval error.
    #[error("interval error: {0}")]
    Interval(coordinate::interval::Error),

    /// An error occurred when constructing a [`ContiguousIntervalPair`].
    #[error("invalid interval pair: {0}")]
    InvalidIntervalPair(interval_pair::Error),

    /// The data section does not cumulative add up the specified end
    /// coordinates. This indicates a malformed data section.
    #[error("misaligned data section, which indicates a malformed chain file")]
    MisalignedDataSection,

    /// A sequence error.
    #[error("sequence error: {0}")]
    Sequence(sequence::Error),
}

/// The core struct using for stepping through contiguous liftover segments.
pub struct StepThroughWithData {
    /// The pointer to the current reference coordinate.
    reference_pointer: Coordinate,

    /// The pointer to the ending reference coordinate.
    expected_reference_end: Coordinate,

    /// The pointer to the current query coordinate.
    query_pointer: Coordinate,

    /// The pointer to the ending query coordinate.
    expected_query_end: Coordinate,

    /// The data records.
    data: Box<dyn Iterator<Item = Record>>,

    /// Whether or not the step-through is finished.
    finished: bool,
}

impl std::fmt::Debug for StepThroughWithData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StepThroughWithData")
            .field("expected_query_end", &self.expected_query_end)
            .field("expected_reference_end", &self.expected_reference_end)
            .field("query_pointer", &self.query_pointer)
            .field("reference_pointer", &self.reference_pointer)
            .finish()
    }
}

impl StepThroughWithData {
    /// Creates a new [`StepThroughWithData`].
    pub(crate) fn new(section: &Section) -> Result<Self, Error> {
        let (reference_pointer, expected_reference_end) = section
            .header()
            .reference_sequence()
            .interval()
            .map_err(Error::Sequence)?
            .into_coordinates();

        let (query_pointer, expected_query_end) = section
            .header()
            .query_sequence()
            .interval()
            .map_err(Error::Sequence)?
            .into_coordinates();

        let data = Box::new(section.data().clone().into_iter());

        Ok(Self {
            query_pointer,
            reference_pointer,
            expected_query_end,
            expected_reference_end,
            data,
            finished: false,
        })
    }

    /// Ensure that the reference pointer and query pointer match the locations
    /// we expect when the stepthrough is finished. If they don't, an error is
    /// thrown.
    fn finish(&mut self) -> Result<(), Error> {
        if self.reference_pointer != self.expected_reference_end {
            return Err(Error::MisalignedDataSection);
        }

        if self.query_pointer != self.expected_query_end {
            return Err(Error::MisalignedDataSection);
        }

        self.finished = true;

        Ok(())
    }
}

impl Drop for StepThroughWithData {
    fn drop(&mut self) {
        if !self.finished && !std::thread::panicking() {
            warn!("a step-through has not been properly finished before dropping!")
        }
    }
}

impl Iterator for StepThroughWithData {
    type Item = Result<(ContiguousIntervalPair, data::Record), Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let chunk = match self.data.next() {
            Some(c) => c,
            None => match self.finish() {
                Ok(_) => return None,
                Err(e) => return Some(Err(e)),
            },
        };
        // (1) Capture reference start, then advance pointer by chunk size.
        let reference_start = self.reference_pointer.clone();

        if !self.reference_pointer.move_forward(chunk.size()) {
            return Some(Err(Error::IntervalStepthroughOutOfBounds(
                String::from("reference forward by chunk size"),
                self.reference_pointer.clone(),
                chunk.size(),
            )));
        }

        let reference_end = self.reference_pointer.clone();
        let reference = match Interval::try_new(reference_start, reference_end) {
            Ok(interval) => interval,
            Err(err) => return Some(Err(Error::Interval(err))),
        };

        // (2) Capture query start, then advance pointer by chunk size.
        let query_start = self.query_pointer.clone();

        if !self.query_pointer.move_forward(chunk.size()) {
            return Some(Err(Error::IntervalStepthroughOutOfBounds(
                String::from("query forward by chunk size"),
                self.query_pointer.clone(),
                chunk.size(),
            )));
        }

        let query_end = self.query_pointer.clone();
        let query = match Interval::try_new(query_start, query_end) {
            Ok(interval) => interval,
            Err(err) => return Some(Err(Error::Interval(err))),
        };

        // (3) Adjust pointers for gaps.
        if let Some(dq) = chunk.dq() {
            if !self.query_pointer.move_forward(dq) {
                return Some(Err(Error::IntervalStepthroughOutOfBounds(
                    String::from("query forward by dq"),
                    self.query_pointer.clone(),
                    dq,
                )));
            }
        }

        if let Some(dt) = chunk.dt() {
            if !self.reference_pointer.move_forward(dt) {
                return Some(Err(Error::IntervalStepthroughOutOfBounds(
                    String::from("reference forward by dt"),
                    self.reference_pointer.clone(),
                    dt,
                )));
            }
        }

        // (4) Create the final contiguous interval pair.
        let pair = match ContiguousIntervalPair::try_new(reference, query) {
            Ok(p) => p,
            Err(err) => return Some(Err(Error::InvalidIntervalPair(err))),
        };

        Some(Ok((pair, chunk)))
    }
}

/// A step-through.
#[derive(Debug)]
pub struct StepThrough(StepThroughWithData);

impl StepThrough {
    /// Creates a new [`StepThrough`].
    pub(crate) fn new(section: &Section) -> Result<Self, Error> {
        StepThroughWithData::new(section).map(Self)
    }
}

impl Iterator for StepThrough {
    type Item = Result<ContiguousIntervalPair, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.0.next()?.map(|(pair, _)| pair))
    }
}

#[cfg(test)]
mod tests {
    use omics::coordinate::Contig;
    use omics::coordinate::Strand;
    use omics::coordinate::interbase::Coordinate;
    use omics::coordinate::interval::interbase::Interval;

    #[test]
    fn test_it_correctly_steps_through_a_valid_chain_file() -> Result<(), Box<dyn std::error::Error>>
    {
        let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
        let mut reader = crate::Reader::new(&data[..]);

        let mut results = reader.sections().map(|x| x.unwrap()).collect::<Vec<_>>();
        assert_eq!(results.len(), 1);

        // SAFETY: we just checked that the length was 1. Thus, this will always
        // unwrap.
        let section = results.pop().unwrap();
        let mut stepthrough = section.stepthrough()?;

        //===============================================//
        // Checks the first interval in the step through //
        //===============================================//

        let result = stepthrough.next().unwrap().unwrap();

        let reference = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 0u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 3u64),
        )?;

        let query = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 5u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 2u64),
        )?;

        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //================================================//
        // Checks the second interval in the step through //
        //================================================//

        let result = stepthrough.next().unwrap().unwrap();

        let reference = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 3u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 4u64),
        )?;

        let query = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 1u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 0u64),
        )?;

        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //==========================================//
        // Checks that the iterator is finished now //
        //==========================================//

        assert!(stepthrough.next().is_none());

        Ok(())
    }

    #[test]
    fn test_it_errors_when_there_is_an_alignment_issue() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 4 + 0 3 seq0 9 - 0 9 1\n2\t0\t1\n1";
        let mut reader = crate::Reader::new(&data[..]);

        let mut results = reader.sections().map(|x| x.unwrap()).collect::<Vec<_>>();
        assert_eq!(results.len(), 1);

        // SAFETY: we just checked that the length was 1. Thus, this will always
        // unwrap.
        let section = results.pop().unwrap();
        let mut stepthrough = section.stepthrough()?;

        //===============================================//
        // Checks the first interval in the step through //
        //===============================================//

        let result = stepthrough.next().unwrap().unwrap();

        let reference = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 0u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 2u64),
        )?;

        let query = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 9u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 7u64),
        )?;

        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //================================================//
        // Checks the second interval in the step through //
        //================================================//

        let result = stepthrough.next().unwrap().unwrap();

        let reference = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 2u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Positive, 3u64),
        )?;

        let query = Interval::try_new(
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 6u64),
            Coordinate::new(Contig::new_unchecked("seq0"), Strand::Negative, 5u64),
        )?;

        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //============================================================//
        // Checks that we now receive a misaligned data section error //
        //============================================================//

        let err = stepthrough.next().unwrap().unwrap_err();
        assert_eq!(
            err.to_string(),
            String::from("misaligned data section, which indicates a malformed chain file")
        );

        Ok(())
    }
}
