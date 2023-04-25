//! A chain file reader.

use std::io::{self, BufRead};
use std::iter;

use crate::line;
use crate::line::Line;

pub mod section;

pub use section::AlignmentDataSection;
pub use section::AlignmentDataSections;

/// An error related to the parsing of a chain file.
#[derive(Debug)]
pub enum ParseError {
    /// An I/O error.
    InputOutput(io::Error),
    /// A line parsing error.
    LineParsing(line::ParseError),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InputOutput(err) => write!(f, "i/o error: {}", err),
            ParseError::LineParsing(err) => write!(f, "line parsing error: {}", err),
        }
    }
}

impl std::error::Error for ParseError {}

/// A chain file reader.
#[derive(Clone, Debug)]
pub struct Reader<T>
where
    T: BufRead,
{
    inner: T,
}

impl<T> Reader<T>
where
    T: BufRead,
{
    /// Creates a chain file reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let reader = chain::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: T) -> Self {
        Self::from(inner)
    }

    /// Gets a reference to the inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use std::io;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    ///
    /// let reader = chain::Reader::new(cursor);
    /// assert_eq!(reader.inner().position(), 0);
    /// ```
    pub fn inner(&self) -> &T {
        &self.inner
    }

    /// Gets a mutable reference to the inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use std::io::Read;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chain::Reader::new(&data[..]);
    /// let mut buffer = vec![0; data.len()];
    ///
    /// reader.inner_mut().read_exact(&mut buffer).unwrap();
    /// assert_eq!(buffer, data[..]);
    /// ```
    pub fn inner_mut(&mut self) -> &mut T {
        &mut self.inner
    }

    /// Consumes self and returns the inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use std::io::BufRead;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let reader = chain::Reader::new(&data[..]);
    /// let mut lines = reader.into_inner().lines().map(|line| line.unwrap());
    ///
    /// assert_eq!(
    ///     lines.next(),
    ///     Some(String::from("chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1"))
    /// );
    /// assert_eq!(lines.next(), Some(String::from("3\t0\t1")));
    /// assert_eq!(lines.next(), Some(String::from("1")));
    /// assert_eq!(lines.next(), None);
    /// ```
    pub fn into_inner(self) -> T {
        self.inner
    }

    /// Reads a raw, textual line from the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use std::io;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chain::Reader::new(&data[..]);
    ///
    /// let mut buffer = String::new();
    ///
    /// assert_eq!(reader.read_line_raw(&mut buffer)?, 36);
    /// assert_eq!(buffer, "chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1");
    ///
    /// assert_eq!(reader.read_line_raw(&mut buffer)?, 6);
    /// assert_eq!(buffer, "3\t0\t1");
    ///
    /// assert_eq!(reader.read_line_raw(&mut buffer)?, 1);
    /// assert_eq!(buffer, "1");
    ///
    /// assert_eq!(reader.read_line_raw(&mut buffer)?, 0);
    ///
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_line_raw(&mut self, buffer: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buffer)
    }

    /// Attempts to read a single `Line` from the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use std::io;
    /// use chain::line::Line;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chain::Reader::new(&data[..]);
    ///
    /// assert!(matches!(reader.read_line()?, Some(Line::Header(_))));
    /// assert!(matches!(reader.read_line()?, Some(Line::AlignmentData(_))));
    /// assert!(matches!(reader.read_line()?, Some(Line::AlignmentData(_))));
    /// assert!(matches!(reader.read_line()?, None));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_line(&mut self) -> Result<Option<Line>, ParseError> {
        let mut buffer = String::new();
        let read = self
            .read_line_raw(&mut buffer)
            .map_err(ParseError::InputOutput)?;

        match read {
            0 => Ok(None),
            _ => {
                let line = buffer.parse::<Line>().map_err(ParseError::LineParsing)?;
                Ok(Some(line))
            }
        }
    }

    /// Returns an iterator over the `Line`s in the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use std::io::BufRead;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chain::Reader::new(&data[..]);
    ///
    /// let lines = reader.lines().collect::<Vec<_>>();
    /// assert_eq!(lines.len(), 3);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn lines(&mut self) -> impl Iterator<Item = io::Result<Line>> + '_ {
        let mut buffer = String::new();

        iter::from_fn(move || {
            buffer.clear();

            match self.read_line_raw(&mut buffer) {
                Ok(0) => None,
                Ok(_) => Some(
                    buffer
                        .parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
                ),
                Err(e) => Some(Err(e)),
            }
        })
    }

    /// Returns an iterator over the alignment data sections in the underlying
    /// reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chain::Reader::new(&data[..]);
    ///
    /// let sections = reader.sections().map(|result| result.unwrap()).collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn sections(&mut self) -> AlignmentDataSections<'_, T> {
        AlignmentDataSections::new(self)
    }
}

impl<T> From<T> for Reader<T>
where
    T: BufRead,
{
    fn from(inner: T) -> Self {
        Self { inner }
    }
}

/// Reads a line from a buffered reader.
///
/// This method is copied almost directly from noodles-gtf. I repurposed it
/// because it captures pretty much exactly what I need to do for this reader.
fn read_line<T>(reader: &mut T, buffer: &mut String) -> io::Result<usize>
where
    T: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    buffer.clear();

    match reader.read_line(buffer) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buffer.ends_with(LINE_FEED) {
                buffer.pop();

                if buffer.ends_with(CARRIAGE_RETURN) {
                    buffer.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    #[test]
    fn test_read_line() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"hello\r\nworld!";
        let mut cursor = io::Cursor::new(data);

        let mut buffer = String::new();
        let len = read_line(&mut cursor, &mut buffer)?;
        assert_eq!(buffer, "hello");
        assert_eq!(len, 7);

        let len = read_line(&mut cursor, &mut buffer)?;
        assert_eq!(buffer, "world!");
        assert_eq!(len, 6);

        Ok(())
    }
}
