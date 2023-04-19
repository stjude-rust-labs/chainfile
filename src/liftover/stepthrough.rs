//! Utilities for stepping through a contiguous liftover segment.
//!
//! A contiguous liftover segment is defined as a pair of intervals (one from
//! the reference and one from the query) that have the same distance and match
//! to each other contiguously.

use std::slice::Iter;

use crate::core::interval;
use crate::core::Coordinate;
use crate::core::Interval;
use crate::core::Strand;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
use crate::reader::AlignmentDataSection;
use crate::record::AlignmentDataRecord;
use crate::record::HeaderRecord;

pub mod interval_pair;

/// An error related to stepping through liftover segments.
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// An error occurred when parsing an [`Interval`].
    InvalidInterval(interval::Error),
    /// An error occurred when constructing a [`ContiguousIntervalPair`].
    InvalidIntervalPair(interval_pair::Error),
    /// The data section does not cumulative add up the specified end
    /// coordinates. This indicates a malformed data section.
    MisalignedDataSection,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidInterval(err) => write!(f, "error when parsing interval: {}", err),
            Error::InvalidIntervalPair(err) => {
                write!(f, "error when constructing interval pair: {}", err)
            }
            Error::MisalignedDataSection => write!(
                f,
                "misaligned data section, which indicates a malformed chain file"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// The core struct using for stepping through contiguous liftover segments.
pub struct StepThrough<'a> {
    query: usize,
    reference: usize,
    header: &'a HeaderRecord,
    data_records: Iter<'a, AlignmentDataRecord>,
}

impl<'a> StepThrough<'a> {
    pub(crate) fn new(section: &'a AlignmentDataSection) -> Self {
        let query = section.header().query_sequence().alignment_start();
        let reference = section.header().reference_sequence().alignment_start();

        Self {
            query,
            reference,
            header: section.header(),
            data_records: section.alignment_data_records().iter(),
        }
    }

    fn finish(&self) -> Result<(), Error> {
        if self.query != self.header.query_sequence().alignment_end() {
            return Err(Error::MisalignedDataSection);
        }

        if self.reference != self.header.reference_sequence().alignment_end() {
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

        // Note that this check against the strand is required because, as
        // stated in the last line of the file format spec, "When the strand
        // value is '-', position coordinates are listed in terms of the
        // reverse-complemented sequence."
        let reference_start_position = match self.header.reference_sequence().strand() {
            Strand::Positive => self.reference,
            Strand::Negative => self.header.reference_sequence().chromosome_size() - self.reference,
        };

        let reference_start = Coordinate::new(
            self.header.reference_sequence().chromosome_name(),
            reference_start_position,
            self.header.reference_sequence().strand().clone(),
        );

        // Note that this check against the strand is required because, as
        // stated in the last line of the file format spec, "When the strand
        // value is '-', position coordinates are listed in terms of the
        // reverse-complemented sequence."
        let reference_end_position = match self.header.reference_sequence().strand() {
            Strand::Positive => self.reference + chunk.size(),
            Strand::Negative => {
                self.header.reference_sequence().chromosome_size() - (self.reference + chunk.size())
            }
        };

        let reference_end = Coordinate::new(
            self.header.reference_sequence().chromosome_name(),
            reference_end_position,
            self.header.reference_sequence().strand().clone(),
        );

        let reference = match Interval::try_new(reference_start, reference_end) {
            Ok(interval) => interval,
            Err(err) => return Some(Err(Error::InvalidInterval(err))),
        };

        // (2) Calculate the coordinates for the query interval.

        // Note that this check against the strand is required because, as
        // stated in the last line of the file format spec, "When the strand
        // value is '-', position coordinates are listed in terms of the
        // reverse-complemented sequence."
        let query_start_position = match self.header.query_sequence().strand() {
            Strand::Positive => self.query,
            Strand::Negative => self.header.query_sequence().chromosome_size() - self.query,
        };

        let query_start = Coordinate::new(
            self.header.query_sequence().chromosome_name(),
            query_start_position,
            self.header.query_sequence().strand().clone(),
        );

        // Note that this check against the strand is required because, as
        // stated in the last line of the file format spec, "When the strand
        // value is '-', position coordinates are listed in terms of the
        // reverse-complemented sequence."
        let query_end_position = match self.header.query_sequence().strand() {
            Strand::Positive => self.query + chunk.size(),
            Strand::Negative => {
                self.header.query_sequence().chromosome_size() - (self.query + chunk.size())
            }
        };

        let query_end = Coordinate::new(
            self.header.query_sequence().chromosome_name(),
            query_end_position,
            self.header.query_sequence().strand().clone(),
        );

        let query = match Interval::try_new(query_start, query_end) {
            Ok(interval) => interval,
            Err(err) => return Some(Err(Error::InvalidInterval(err))),
        };

        // (3) Adjust the pointers to the current query and reference locations.
        // (a) Always add the chunk length to both query and reference.
        self.query += chunk.size();
        self.reference += chunk.size();

        // (b) If there is a query offset, add it to the query pointer.
        if let Some(dq) = chunk.dq() {
            self.query += *dq;
        }

        // (c) If there is a reference offset, add it to the reference pointer.
        if let Some(dt) = chunk.dt() {
            self.reference += *dt;
        }

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
        let mut stepthrough = section.stepthrough();

        //===============================================//
        // Checks the first interval in the step through //
        //===============================================//

        let reference = Interval::try_new(
            Coordinate::new("seq0", 0, Strand::Positive),
            Coordinate::new("seq0", 3, Strand::Positive),
        )?;
        let query = Interval::try_new(
            Coordinate::new("seq0", 5, Strand::Negative),
            Coordinate::new("seq0", 2, Strand::Negative),
        )?;

        let result = stepthrough.next().unwrap().unwrap();
        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //================================================//
        // Checks the second interval in the step through //
        //================================================//

        let reference = Interval::try_new(
            Coordinate::new("seq0", 3, Strand::Positive),
            Coordinate::new("seq0", 4, Strand::Positive),
        )?;
        let query = Interval::try_new(
            Coordinate::new("seq0", 1, Strand::Negative),
            Coordinate::new("seq0", 0, Strand::Negative),
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
        let mut stepthrough = section.stepthrough();

        //===============================================//
        // Checks the first interval in the step through //
        //===============================================//

        let reference = Interval::try_new(
            Coordinate::new("seq0", 0, Strand::Positive),
            Coordinate::new("seq0", 2, Strand::Positive),
        )?;
        let query = Interval::try_new(
            Coordinate::new("seq0", 10, Strand::Negative),
            Coordinate::new("seq0", 8, Strand::Negative),
        )?;

        let result = stepthrough.next().unwrap().unwrap();
        assert_eq!(result.reference(), &reference);
        assert_eq!(result.query(), &query);

        //================================================//
        // Checks the second interval in the step through //
        //================================================//

        let reference = Interval::try_new(
            Coordinate::new("seq0", 2, Strand::Positive),
            Coordinate::new("seq0", 3, Strand::Positive),
        )?;
        let query = Interval::try_new(
            Coordinate::new("seq0", 7, Strand::Negative),
            Coordinate::new("seq0", 6, Strand::Negative),
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
