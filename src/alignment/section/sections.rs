//! An iterator over a set of [alignment sections](crate::alignment::Section)s.

use std::io::BufRead;

use crate::Line;
use crate::Reader;
use crate::alignment::Section;
use crate::alignment::section::Builder;
use crate::alignment::section::builder;
use crate::alignment::section::data;
use crate::alignment::section::data::record::Kind;
use crate::alignment::section::header;
use crate::reader;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// An error related to the parsing of an alignment section.
#[derive(Debug)]
pub enum ParseError {
    /// The file abruptly ended before closing an alignment section.
    AbruptEndInSection,

    /// There was a blank line within an alignment section.
    BlankLineInSection(usize),

    /// Alignment data was found in between alignment sections.
    DataBetweenSections(data::Record),

    /// There was a header record within an alignment section.
    HeaderInSection(header::Record),

    /// There was an issue reading from the underlying reader.
    Reader(reader::Error),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::AbruptEndInSection => {
                write!(
                    f,
                    "the file abruptly ended in the middle of an alignment section"
                )
            }
            ParseError::BlankLineInSection(line_no) => {
                write!(f, "found blank line in alignment section: line {line_no}")
            }
            ParseError::DataBetweenSections(record) => {
                write!(f, "found alignment data between sections: record: {record}")
            }
            ParseError::HeaderInSection(record) => {
                write!(f, "found header in alignment section: record: {record}")
            }
            ParseError::Reader(err) => write!(f, "reader error: {err}"),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to [`Sections`].
#[derive(Debug)]
pub enum Error {
    /// A builder error.
    Builder(builder::Error),

    /// A parse error.
    Parse(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Builder(err) => write!(f, "builder error: {err}"),
            Error::Parse(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// Sections
////////////////////////////////////////////////////////////////////////////////////////

/// The state of the iterator.
#[derive(Debug)]
enum State {
    /// The reader is in-between alignment sections.
    InBetweenSections,

    /// The reader is in the middle of an alignment section.
    ReadingSection,
}

/// An iterator struct the traverses the alignment sections while keeping
/// track of the current state of the reader.
#[derive(Debug)]
pub struct Sections<'a, T>
where
    T: BufRead,
{
    /// The inner reader.
    reader: &'a mut Reader<T>,

    /// The state of the iterator.
    state: State,

    /// The line number.
    line_no: usize,
}

impl<'a, T> Sections<'a, T>
where
    T: BufRead,
{
    /// Creates a new [`Sections`].
    pub(crate) fn new(reader: &'a mut Reader<T>) -> Self {
        Self {
            reader,
            state: State::InBetweenSections,
            line_no: 0usize,
        }
    }
}

impl<T> Iterator for Sections<'_, T>
where
    T: BufRead,
{
    type Item = Result<Section>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut builder: Option<Builder> = None;
        let mut buffer = String::new();

        loop {
            // (1) Reads the current line and returns an error if a parsing
            // error occurs within the reader.
            let line = match self.reader.read_line(&mut buffer) {
                Ok(l) => l,
                Err(err) => return Some(Err(Error::Parse(ParseError::Reader(err)))),
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
                        return Some(Err(Error::Parse(ParseError::AbruptEndInSection)));
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
            builder = match builder {
                Some(inner_builder) => {
                    match line {
                        Line::AlignmentData(record) => Some(inner_builder.push_data(record)),
                        Line::Empty => Some(inner_builder),
                        // The `get_state()` method will throw an error in the step above
                        // before this code is even reached. As such, this is unreachable.
                        Line::Header(_) => unreachable!(),
                    }
                }
                None => {
                    match line {
                        // The `get_state()` method will throw an error in the step above
                        // before this code is even reached. As such, this is unreachable.
                        Line::AlignmentData(_) => unreachable!(),
                        Line::Empty => None,
                        Line::Header(record) => {
                            // SAFETY: since this is the first time `header()`
                            // is called (right after the creation of the
                            // builder), it will never error because of multiple
                            // headers being provided.
                            Some(Builder::default().header(record).unwrap())
                        }
                    }
                }
            };

            // (5) Performs the associated action given our current state.
            match self.state {
                State::InBetweenSections => {
                    if let Some(builder) = builder {
                        match builder.try_build() {
                            Ok(section) => return Some(Ok(section)),
                            Err(err) => return Some(Err(Error::Builder(err))),
                        }
                    }
                }
                State::ReadingSection => {}
            }
        }
    }
}

/// Gets the current state given the previous state and the line that was just
/// read in from the [`Reader`].
#[allow(clippy::result_large_err)]
fn get_state(last: &State, line: &Line, line_no: usize) -> Result<State> {
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
            Err(Error::Parse(ParseError::DataBetweenSections(
                record.clone(),
            )))
        }
        (State::ReadingSection, Line::Empty) => {
            // We found a blank line within an alignment section, this is
            // not allowed in the format.
            Err(Error::Parse(ParseError::BlankLineInSection(line_no)))
        }
        (State::ReadingSection, Line::Header(record)) => {
            // We found a header within a alignment section, this is not
            // allowed in the format.
            Err(Error::Parse(ParseError::HeaderInSection(record.clone())))
        }
        (State::ReadingSection, Line::AlignmentData(record)) => {
            // We found alignment data in a section. We need to examine the
            // alignment data record that was parsed to determine if the record
            // is terminating for the section or not.
            match record.kind() {
                Kind::NonTerminating => Ok(State::ReadingSection),
                Kind::Terminating => Ok(State::InBetweenSections),
            }
        }
    }
}
