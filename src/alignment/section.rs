//! Sections of alignment data within a chain file.

use nonempty::NonEmpty;

use crate::alignment::section::header::Sequence;
use crate::liftover::StepThroughWithData;
use crate::liftover::stepthrough;
use crate::liftover::stepthrough::StepThrough;

mod builder;
pub mod data;
pub mod header;
pub mod sections;

pub use builder::Builder;
pub use sections::Sections;

/// An alignment section.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Section {
    /// The header record.
    header: header::Record,

    /// The data records.
    data: NonEmpty<data::Record>,
}

impl Section {
    /// Gets the header record for the [`Section`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io::BufRead;
    /// use std::io::{self};
    ///
    /// use chainfile::alignment::section::header::Record;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chainfile::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// assert_eq!(
    ///     section.header(),
    ///     &"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1".parse::<Record>()?
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn header(&self) -> &header::Record {
        &self.header
    }

    /// Gets the data records for the [`Section`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    /// use std::io::BufRead;
    /// use std::io::{self};
    ///
    /// use chainfile::alignment::section::data::Record;
    /// use chainfile::alignment::section::data::record::Kind;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chainfile::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// let mut records = section
    ///     .data()
    ///     .into_iter()
    ///     .cloned()
    ///     .collect::<VecDeque<Record>>();
    /// assert_eq!(records.len(), 2);
    ///
    /// // SAFETY: we just checked that the length was two, so this first pop will succeed.
    /// let record = records.pop_front().unwrap();
    /// assert_eq!(
    ///     record,
    ///     Record::try_new(3, Some(0), Some(1), Kind::NonTerminating)?
    /// );
    ///
    /// // SAFETY: we just checked that the length was two, so this second pop will succeed.
    /// let record = records.pop_front().unwrap();
    /// assert_eq!(record, Record::try_new(1, None, None, Kind::Terminating)?);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn data(&self) -> &NonEmpty<data::Record> {
        &self.data
    }

    /// Gets a liftover [`StepThrough`] for this alignment [`Section`].
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    /// use std::io::BufRead;
    /// use std::io::{self};
    ///
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::Strand;
    /// use omics::coordinate::interval::interbase::Interval;
    /// use omics::coordinate::position::interbase::Position;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chainfile::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    ///
    /// let mut stepthrough = section.stepthrough()?.collect::<VecDeque<_>>();
    /// assert_eq!(stepthrough.len(), 2);
    ///
    /// // SAFETY: we just checked that the length was two, so this first pop will succeed.
    /// let result = stepthrough.pop_front().unwrap()?;
    /// let expected = ContiguousIntervalPair::try_new(
    ///     "seq0:+:0-3".parse::<Interval>()?,
    ///     "seq0:-:5-2".parse::<Interval>()?,
    /// )?;
    /// assert_eq!(result, expected);
    ///
    /// // SAFETY: we just checked that the length was two, so this second pop will succeed.
    /// let result = stepthrough.pop_front().unwrap()?;
    /// let expected = ContiguousIntervalPair::try_new(
    ///     "seq0:+:3-4".parse::<Interval>()?,
    ///     "seq0:-:1-0".parse::<Interval>()?,
    /// )?;
    /// assert_eq!(result, expected);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn stepthrough(&self) -> Result<StepThrough, stepthrough::Error> {
        StepThrough::new(self)
    }

    /// Gets a liftover [`StepThrough`] for this alignment [`Section`].
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    /// use std::io::BufRead;
    /// use std::io::{self};
    ///
    /// use chainfile::alignment::section::data::record::Kind;
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::Strand;
    /// use omics::coordinate::interval::interbase::Interval;
    /// use omics::coordinate::position::interbase::Position;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chainfile::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    ///
    /// let mut stepthrough = section.stepthrough_with_data()?.collect::<VecDeque<_>>();
    /// assert_eq!(stepthrough.len(), 2);
    ///
    /// // SAFETY: we just checked that the length was two, so this first pop will succeed.
    /// let (result, record) = stepthrough.pop_front().unwrap()?;
    /// let expected = ContiguousIntervalPair::try_new(
    ///     "seq0:+:0-3".parse::<Interval>()?,
    ///     "seq0:-:5-2".parse::<Interval>()?,
    /// )?;
    /// assert_eq!(result, expected);
    /// assert_eq!(record.size(), 3);
    /// assert_eq!(record.dt(), Some(0));
    /// assert_eq!(record.dq(), Some(1));
    /// assert_eq!(record.kind(), Kind::NonTerminating);
    ///
    /// // SAFETY: we just checked that the length was two, so this second pop will succeed.
    /// let (result, record) = stepthrough.pop_front().unwrap()?;
    /// let expected = ContiguousIntervalPair::try_new(
    ///     "seq0:+:3-4".parse::<Interval>()?,
    ///     "seq0:-:1-0".parse::<Interval>()?,
    /// )?;
    /// assert_eq!(result, expected);
    /// assert_eq!(record.size(), 1);
    /// assert_eq!(record.dt(), None);
    /// assert_eq!(record.dq(), None);
    /// assert_eq!(record.kind(), Kind::Terminating);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn stepthrough_with_data(&self) -> Result<StepThroughWithData, stepthrough::Error> {
        StepThroughWithData::new(self)
    }

    /// Gets the reference sequence.
    pub fn reference_sequence(&self) -> &Sequence {
        self.header.reference_sequence()
    }

    /// Gets the query sequence.
    pub fn query_sequence(&self) -> &Sequence {
        self.header.query_sequence()
    }
}

#[cfg(test)]
mod tests {
    use std::io;

    use crate::Reader;

    #[test]
    fn test_valid_sections() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let results = reader
            .sections()
            .map(|result| result.unwrap())
            .collect::<Vec<_>>();

        assert_eq!(results.len(), 1);
        assert_eq!(results.first().unwrap().data().len(), 2);

        Ok(())
    }

    #[test]
    fn test_invalid_sections_alignment_data_between_sections()
    -> Result<(), Box<dyn std::error::Error>> {
        let data = b"2\t0\t0";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let error = reader.sections().next().unwrap().unwrap_err();
        assert_eq!(
            error.to_string(),
            "parse error: found alignment data between sections: record: 2\t0\t0"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_sections_abrupt_end_in_section() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1\t0\t0";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let error = reader.sections().next().unwrap().unwrap_err();
        assert_eq!(
            error.to_string(),
            "parse error: the file abruptly ended in the middle of an alignment section"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_sections_blank_line_in_section() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1\t0\t0\n\n0";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let error = reader.sections().next().unwrap().unwrap_err();
        assert_eq!(
            error.to_string(),
            "parse error: found blank line in alignment section: line 4"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_sections_header_in_section() -> Result<(), Box<dyn std::error::Error>> {
        let data =
            b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1\t0\t0\nchain 0 seq0 2 + 0 2 seq0 2 - 0 2 1\n0";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let error = reader.sections().next().unwrap().unwrap_err();
        assert_eq!(
            error.to_string(),
            "parse error: found header in alignment section: record: chain 0 seq0 2 + 0 2 seq0 2 \
             - 0 2 1"
        );

        Ok(())
    }
}
