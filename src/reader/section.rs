//! Utilities related to alignment sections within a chain file.

use std::io::BufRead;

use crate::line::Line;
use crate::reader;
use crate::record::alignment_data::AlignmentDataRecordType;
use crate::record::AlignmentDataRecord;
use crate::record::HeaderRecord;
use crate::Reader;

/// An alignment data section.
#[derive(Clone, Debug)]
pub struct AlignmentDataSection {
    header: HeaderRecord,
    alignment_data_records: Vec<AlignmentDataRecord>,
}

impl AlignmentDataSection {
    fn new(header: HeaderRecord) -> Self {
        Self {
            header,
            alignment_data_records: Default::default(),
        }
    }

    fn add_alignment_data_record(&mut self, record: AlignmentDataRecord) {
        self.alignment_data_records.push(record);
    }

    /// Gets the header for the alignment data section.
    pub fn header(&self) -> &HeaderRecord {
        &self.header
    }

    /// Gets the alignment data records for the alignment data section.
    pub fn alignment_data_records(&self) -> &Vec<AlignmentDataRecord> {
        &self.alignment_data_records
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
    type Item = Result<AlignmentDataSection, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut section: Option<AlignmentDataSection> = None;

        loop {
            // (1) Reads the current line and returns an error if a parsing
            // error occurs within the reader.
            let line = match self.reader.read_line() {
                Ok(l) => l,
                Err(err) => return Some(Err(ParseError::Reader(err))),
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
                    State::ReadingSection => return Some(Err(ParseError::AbruptEndInSection)),
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

fn get_state(last: &State, line: &Line, line_no: usize) -> Result<State, ParseError> {
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
            Err(ParseError::AlignmentDataBetweenSections(record.clone()))
        }
        (State::ReadingSection, Line::Empty) => {
            // We found a blank line within an alignment data section, this is
            // not allowed in the format.
            Err(ParseError::BlankLineInSection(line_no))
        }
        (State::ReadingSection, Line::Header(record)) => {
            // We found a header within a alignment data section, this is not
            // allowed in the format.
            Err(ParseError::HeaderInSection(record.clone()))
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
        let data = b"chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1\n2\t0\t0\n0";
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
        let data = b"chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1\n2\t0\t0";
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
        let data = b"chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1\n2\t0\t0\n\n0";
        let cursor = io::Cursor::new(data);

        let mut reader = Reader::new(cursor);
        let error = reader.sections().next().unwrap().unwrap_err();
        assert_eq!(
            error.to_string(),
            "found blank line in alignment data section: line 3"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_sections_header_in_section() -> Result<(), Box<dyn std::error::Error>> {
        let data =
            b"chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1\n2\t0\t0\nchain 0 seq0 2 + 0 2 seq0 2 - 0 2 1\n0";
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