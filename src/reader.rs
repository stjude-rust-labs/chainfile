//! A chain file reader.

use std::io::BufRead;
use std::io::{self};
use std::iter;

use crate::Line;
use crate::alignment::section::Sections;
use crate::line;

/// The new line character.
const NEW_LINE: char = '\n';

/// The carriage return character.
const CARRIAGE_RETURN: char = '\r';

/// An error related to a [`Reader`].
#[derive(Debug)]
pub enum Error {
    /// An I/O error.
    Io(io::Error),

    /// A line error.
    Line(line::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(err) => write!(f, "i/o error: {err}"),
            Error::Line(err) => write!(f, "line error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A chain file reader.
#[derive(Clone, Debug)]
pub struct Reader<T>(T)
where
    T: BufRead;

impl<T> Reader<T>
where
    T: BufRead,
{
    /// Creates a chain file reader.
    ///
    /// # Examples
    ///
    /// ```
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let reader = chainfile::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: T) -> Self {
        Self::from(inner)
    }

    /// Gets a reference to the inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let cursor = io::Cursor::new(data);
    ///
    /// let reader = chainfile::Reader::new(cursor);
    /// assert_eq!(reader.inner().position(), 0);
    /// ```
    pub fn inner(&self) -> &T {
        &self.0
    }

    /// Gets a mutable reference to the inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io::Read;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chainfile::Reader::new(&data[..]);
    /// let mut buffer = vec![0; data.len()];
    ///
    /// reader.inner_mut().read_exact(&mut buffer).unwrap();
    /// assert_eq!(buffer, data[..]);
    /// ```
    pub fn inner_mut(&mut self) -> &mut T {
        &mut self.0
    }

    /// Consumes self and returns the inner reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io::BufRead;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let reader = chainfile::Reader::new(&data[..]);
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
        self.0
    }

    /// Reads a raw, textual line from the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chainfile::Reader::new(&data[..]);
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
        read_line(self.inner_mut(), buffer)
    }

    /// Attempts to read a [`Line`] from the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io;
    ///
    /// use chainfile::Line;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chainfile::Reader::new(&data[..]);
    ///
    /// let mut buffer = String::new();
    /// assert!(matches!(
    ///     reader.read_line(&mut buffer)?,
    ///     Some(Line::Header(_))
    /// ));
    /// assert!(matches!(
    ///     reader.read_line(&mut buffer)?,
    ///     Some(Line::AlignmentData(_))
    /// ));
    /// assert!(matches!(
    ///     reader.read_line(&mut buffer)?,
    ///     Some(Line::AlignmentData(_))
    /// ));
    /// assert!(matches!(reader.read_line(&mut buffer)?, None));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_line(&mut self, buffer: &mut String) -> Result<Option<Line>, Error> {
        let read = self.read_line_raw(buffer).map_err(Error::Io)?;

        match read {
            0 => Ok(None),
            _ => {
                let line = buffer.parse::<Line>().map_err(Error::Line)?;
                Ok(Some(line))
            }
        }
    }

    /// Returns an iterator over the `Line`s in the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::io::BufRead;
    ///
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chainfile::Reader::new(&data[..]);
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

    /// Returns an iterator over the alignment sections in the underlying
    /// reader.
    ///
    /// # Examples
    ///
    /// ```
    /// let data = b"chain 0 seq0 4 + 0 4 seq0 5 - 0 5 1\n3\t0\t1\n1";
    /// let mut reader = chainfile::Reader::new(&data[..]);
    ///
    /// let sections = reader
    ///     .sections()
    ///     .map(|result| result.unwrap())
    ///     .collect::<Vec<_>>();
    /// assert_eq!(sections.len(), 1);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn sections(&mut self) -> Sections<'_, T> {
        Sections::new(self)
    }
}

impl<T> From<T> for Reader<T>
where
    T: BufRead,
{
    fn from(inner: T) -> Self {
        Self(inner)
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
    buffer.clear();

    match reader.read_line(buffer) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buffer.ends_with(NEW_LINE) {
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
    use std::io;

    use super::*;

    #[test]
    fn test_read_line() {
        let data = b"hello\r\nworld!";
        let mut cursor = io::Cursor::new(data);

        let mut buffer = String::new();
        let len = read_line(&mut cursor, &mut buffer).unwrap();
        assert_eq!(buffer, "hello");
        assert_eq!(len, 7);

        let len = read_line(&mut cursor, &mut buffer).unwrap();
        assert_eq!(buffer, "world!");
        assert_eq!(len, 6);
    }
}
