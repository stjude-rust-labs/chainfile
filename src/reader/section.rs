//! Utilities related to alignment sections within a chain file.

use std::io::BufRead;

use crate::core::coordinate;
use crate::core::Coordinate;
use crate::core::Strand;
use crate::liftover::stepthrough;
use crate::liftover::StepThrough;
use crate::line::Line;
use crate::reader;
use crate::record::alignment_data::AlignmentDataRecordType;
use crate::record::header::Sequence;
use crate::record::AlignmentDataRecord;
use crate::record::HeaderRecord;
use crate::Reader;

/// An alignment data section.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentDataSection {
    header: HeaderRecord,
    alignment_data_records: Vec<AlignmentDataRecord>,
}

impl AlignmentDataSection {
    /// Creates a new [`AlignmentDataSection`].
    fn new(header: HeaderRecord) -> Self {
        Self {
            header,
            alignment_data_records: Default::default(),
        }
    }

    /// Adds an alignment record to the [`AlignmentDataSection`].
    fn add_alignment_data_record(&mut self, record: AlignmentDataRecord) {
        self.alignment_data_records.push(record);
    }

    /// Gets the header for the alignment data section.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::record::HeaderRecord;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// assert_eq!(section.header(), &"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1".parse::<HeaderRecord>()?);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn header(&self) -> &HeaderRecord {
        &self.header
    }

    /// Gets the alignment data records for the alignment data section.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    ///
    /// use chainfile as chain;
    /// use chain::record::AlignmentDataRecord;
    /// use chain::record::alignment_data::AlignmentDataRecordType;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// let mut records = section.alignment_data_records()
    ///     .into_iter()
    ///     .cloned()
    ///     .collect::<VecDeque<AlignmentDataRecord>>();
    /// assert_eq!(records.len(), 2);
    ///
    /// // SAFETY: we just checked that the length was two, so this first pop will succeed.
    /// let record = records.pop_front().unwrap();
    /// assert_eq!(
    ///     record,
    ///     AlignmentDataRecord::try_new(3, Some(0), Some(1), AlignmentDataRecordType::NonTerminating)?
    /// );
    ///
    /// // SAFETY: we just checked that the length was two, so this second pop will succeed.
    /// let record = records.pop_front().unwrap();
    /// assert_eq!(
    ///     record,
    ///     AlignmentDataRecord::try_new(1, None, None, AlignmentDataRecordType::Terminating)?
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alignment_data_records(&self) -> &Vec<AlignmentDataRecord> {
        &self.alignment_data_records
    }

    /// Gets a liftover step through for this alignment data section.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::collections::VecDeque;
    ///
    /// use chainfile as chain;
    /// use chain::core::Interval;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    /// use chain::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
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
    /// let expected = ContiguousIntervalPair::try_new("seq0:0-3".parse::<Interval>()?, "seq0:4-1".parse::<Interval>()?)?;
    /// assert_eq!(result, expected);
    ///
    /// // SAFETY: we just checked that the length was two, so this second pop will succeed.
    /// let result = stepthrough.pop_front().unwrap()?;
    /// let expected = ContiguousIntervalPair::try_new("seq0:3-4".parse::<Interval>()?, "seq0:0-$".parse::<Interval>()?)?;
    /// assert_eq!(result, expected);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn stepthrough(&self) -> Result<StepThrough<'_>, stepthrough::Error> {
        StepThrough::new(self)
    }

    /// Gets the starting coordinate for the reference interval specified in the
    /// section header.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// let start = section.reference_start_coordinate()?;
    /// assert_eq!(start.contig(), &String::from("seq0"));
    /// assert_eq!(start.position(), &Position::new(0));
    /// assert_eq!(start.strand(), &Strand::Positive);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_start_coordinate(&self) -> Result<Coordinate, coordinate::Error> {
        start_coordinate(self.header.reference_sequence())
    }

    /// Gets the ending coordinate for the reference interval specified in the
    /// section header.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// let start = section.reference_end_coordinate()?;
    /// assert_eq!(start.contig(), &String::from("seq0"));
    /// assert_eq!(start.position(), &Position::new(4));
    /// assert_eq!(start.strand(), &Strand::Positive);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_end_coordinate(&self) -> Result<Coordinate, coordinate::Error> {
        end_coordinate(self.header.reference_sequence())
    }

    /// Gets the starting coordinate for the query interval specified in the
    /// section header.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// let start = section.query_start_coordinate()?;
    /// assert_eq!(start.contig(), &String::from("seq0"));
    /// assert_eq!(start.position(), &Position::new(4));
    /// assert_eq!(start.strand(), &Strand::Negative);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query_start_coordinate(&self) -> Result<Coordinate, coordinate::Error> {
        start_coordinate(self.header.query_sequence())
    }

    /// Gets the ending coordinate for the query interval specified in the
    /// section header.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    /// use std::io::{self, BufRead};
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    /// let mut reader = chain::Reader::new(cursor);
    ///
    /// let mut sections = reader.sections().collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// // SAFETY: we just checked that the length was one.
    /// let section = sections.pop().unwrap()?;
    /// let start = section.query_end_coordinate()?;
    /// assert_eq!(start.contig(), &String::from("seq0"));
    /// assert_eq!(start.position(), &Position::negative_bound());
    /// assert_eq!(start.strand(), &Strand::Negative);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query_end_coordinate(&self) -> Result<Coordinate, coordinate::Error> {
        end_coordinate(self.header.query_sequence())
    }
}

fn start_coordinate(sequence: &Sequence) -> Result<Coordinate, coordinate::Error> {
    let coordinate = Coordinate::try_new(
        sequence.chromosome_name(),
        sequence.alignment_start(),
        sequence.strand().clone(),
    )?;

    match sequence.strand() {
        Strand::Positive => Ok(coordinate),
        Strand::Negative => coordinate.compliment_position(sequence.chromosome_size()),
    }
}

fn end_coordinate(sequence: &Sequence) -> Result<Coordinate, coordinate::Error> {
    let coordinate = Coordinate::try_new(
        sequence.chromosome_name(),
        sequence.alignment_end(),
        sequence.strand().clone(),
    )?;

    match sequence.strand() {
        Strand::Positive => Ok(coordinate),
        Strand::Negative => coordinate.compliment_position(sequence.chromosome_size()),
    }
}
#[derive(Debug)]
enum State {
    InBetweenSections,
    ReadingSection,
}

/// An error related to the parsing of an alignment data section.
#[derive(Debug)]
pub enum ParseError {
    /// Alignment data was found in between alignment data sections.
    AlignmentDataBetweenSections(AlignmentDataRecord),
    /// The file abruptly ended before closing an alignment data section.
    AbruptEndInSection,
    /// There was a blank line within an alignment data section.
    BlankLineInSection(usize),
    /// There was a header record within an alignment data section.
    HeaderInSection(HeaderRecord),
    /// There was an issue reading from the underlying reader.
    Reader(reader::ParseError),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::AlignmentDataBetweenSections(record) => write!(
                f,
                "found alignment data between sections\n\nrecord: {}",
                record
            ),
            ParseError::AbruptEndInSection => {
                write!(
                    f,
                    "the file abruptly ended in the middle of an alignment data section"
                )
            }
            ParseError::BlankLineInSection(line_no) => {
                write!(
                    f,
                    "found blank line in alignment data section: line {}",
                    line_no
                )
            }
            ParseError::HeaderInSection(record) => write!(
                f,
                "found header in alignment data section\n\nrecord: {}",
                record
            ),
            ParseError::Reader(err) => write!(f, "reader error: {}", err),
        }
    }
}

impl std::error::Error for ParseError {}

/// An iterator struct the traverses the alignment data sections while keeping
/// track of the current state of the reader.
pub struct AlignmentDataSections<'a, T>
where
    T: BufRead,
{
    reader: &'a mut Reader<T>,
    state: State,
    line_no: usize,
}

impl<'a, T> AlignmentDataSections<'a, T>
where
    T: BufRead,
{
    pub(super) fn new(reader: &'a mut Reader<T>) -> Self {
        Self {
            reader,
            state: State::InBetweenSections,
            line_no: 0usize,
        }
    }
}

impl<'a, T> Iterator for AlignmentDataSections<'a, T>
where
    T: BufRead,
{
    type Item = Result<AlignmentDataSection, Box<ParseError>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut section: Option<AlignmentDataSection> = None;

        loop {
            // (1) Reads the current line and returns an error if a parsing
            // error occurs within the reader.
            let line = match self.reader.read_line() {
                Ok(l) => l,
                Err(err) => return Some(Err(Box::new(ParseError::Reader(err)))),
            };

            self.line_no += 1;

            // (2) Checks to see if the line read returned a result. If the line
            // read _did not_ return a result, we either exit gracefully (in the
            // case that we are between sections) or error out (in the case we
            // are in the middle of reading a section).
            let line = match line {
                Some(l) => l,
                None => match self.state {
                    State::InBetweenSections => return None,
                    State::ReadingSection => {
                        return Some(Err(Box::new(ParseError::AbruptEndInSection)))
                    }
                },
            };

            // (3) Gets the current state and errors out if we encounter a parse error.
            self.state = match get_state(&self.state, &line, self.line_no) {
                Ok(s) => s,
                Err(err) => return Some(Err(err)),
            };

            // (4) If the section is instantiated and if we have an alignment
            // record line, add the line to the section.
            if let Some(ref mut s) = section {
                match line {
                    Line::AlignmentData(record) => s.add_alignment_data_record(record),
                    _ => unreachable!(),
                }
            } else {
                match line {
                    Line::AlignmentData(_) => unreachable!(),
                    Line::Empty => {}
                    Line::Header(header) => section = Some(AlignmentDataSection::new(header)),
                }
            }

            // (5) Performs the associated action given our current state.
            match self.state {
                State::InBetweenSections => {
                    if let Some(section) = section {
                        return Some(Ok(section));
                    }
                }
                State::ReadingSection => {}
            }
        }
    }
}

fn get_state(last: &State, line: &Line, line_no: usize) -> Result<State, Box<ParseError>> {
    match (last, line) {
        (State::InBetweenSections, Line::Empty) => {
            // We're in between sections and we hit an empty line, we can just
            // ignore this and move on to the next line.
            Ok(State::InBetweenSections)
        }
        (State::InBetweenSections, Line::Header(_)) => {
            // We were in-between sections, but now we hit a header line which
            // indicates the start of a new section.
            Ok(State::ReadingSection)
        }
        (State::InBetweenSections, Line::AlignmentData(record)) => {
            // We found alignment data in between sections, this is not allowed
            // in the format.
            Err(Box::new(ParseError::AlignmentDataBetweenSections(
                record.clone(),
            )))
        }
        (State::ReadingSection, Line::Empty) => {
            // We found a blank line within an alignment data section, this is
            // not allowed in the format.
            Err(Box::new(ParseError::BlankLineInSection(line_no)))
        }
        (State::ReadingSection, Line::Header(record)) => {
            // We found a header within a alignment data section, this is not
            // allowed in the format.
            Err(Box::new(ParseError::HeaderInSection(record.clone())))
        }
        (State::ReadingSection, Line::AlignmentData(record)) => {
            // We found alignment data in a section. We need to examine the
            // alignment data record that was parsed to determine if the record
            // is terminating for the section or not.
            match record.record_type() {
                AlignmentDataRecordType::NonTerminating => Ok(State::ReadingSection),
                AlignmentDataRecordType::Terminating => Ok(State::InBetweenSections),
            }
        }
    }
}

#[cfg(test)]
pub mod tests {
    use std::io;

    use super::*;

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
        assert_eq!(results.get(0).unwrap().alignment_data_records().len(), 2);
        Ok(())
    }

    #[test]
    fn test_invalid_sections_alignment_data_between_sections(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let data = b"2\t0\t0";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let error = reader.sections().next().unwrap().unwrap_err();
        assert_eq!(
            error.to_string(),
            "found alignment data between sections\n\nrecord: 2\t0\t0"
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
            "the file abruptly ended in the middle of an alignment data section"
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
            "found blank line in alignment data section: line 4"
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
            "found header in alignment data section\n\nrecord: chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1"
        );

        Ok(())
    }
}
