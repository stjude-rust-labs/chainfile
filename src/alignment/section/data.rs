//! An alignment data record.

use std::num::ParseIntError;
use std::str::FromStr;

use crate::alignment::section::data::record::Kind;

pub mod record;

/// The delimiter for an alignment data record.
const ALIGNMENT_DATA_DELIMITER: char = '\t';

/// The number of expected fields in a non-terminating alignment data record.
pub const NUM_ALIGNMENT_DATA_FIELDS_NONTERMINATING: usize = 3;

/// The number of expected fields in a terminating alignment data record.
pub const NUM_ALIGNMENT_DATA_FIELDS_TERMINATING: usize = 1;

/// An error related to the parsing of an alignment data record.
#[derive(Debug)]
pub enum ParseError {
    /// An incorrect number of fields in the alignment data line.
    IncorrectNumberOfFields(usize),

    /// An invalid size.
    InvalidSize(ParseIntError),

    /// An invalid dt.
    InvalidDt(ParseIntError),

    /// An invalid dq.
    InvalidDq(ParseIntError),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::IncorrectNumberOfFields(n) => write!(
                f,
                "invalid number of fields in alignment data: expected {} (non-terminating) or {} \
                 (terminating) fields, found {} fields",
                NUM_ALIGNMENT_DATA_FIELDS_NONTERMINATING, NUM_ALIGNMENT_DATA_FIELDS_TERMINATING, n
            ),
            ParseError::InvalidSize(err) => write!(f, "invalid size: {}", err),
            ParseError::InvalidDt(err) => write!(f, "invalid dt: {}", err),
            ParseError::InvalidDq(err) => write!(f, "invalid dq: {}", err),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a [`Record`].
#[derive(Debug)]
pub enum Error {
    /// An invalid dt value for a non-terminating alignment data line.
    InvalidNonTerminatingDt,

    /// An invalid dq value for a non-terminating alignment data line.
    InvalidNonTerminatingDq,

    /// An invalid dt value for a terminating alignment data line.
    InvalidTerminatingDt,

    /// An invalid dq value for a terminating alignment data line.
    InvalidTerminatingDq,

    /// A parse error.
    Parse(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidNonTerminatingDt => write!(
                f,
                "expected value for dt in non-terminating alignment data line, found no value"
            ),
            Error::InvalidNonTerminatingDq => write!(
                f,
                "expected value for dq in non-terminating alignment data line, found no value"
            ),
            Error::InvalidTerminatingDt => write!(
                f,
                "expected no value for dt in terminating alignment data line, found value"
            ),
            Error::InvalidTerminatingDq => write!(
                f,
                "expected no value for dq in terminating alignment data line, found value"
            ),
            Error::Parse(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

/// An alignment data record within a chain file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record(usize, Option<usize>, Option<usize>, Kind);

impl Record {
    /// Attempts to create a new [`Record`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::data::Record;
    /// use chainfile::alignment::section::data::record::Kind;
    ///
    /// let record = Record::try_new(10, Some(0), Some(1), Kind::NonTerminating)?;
    ///
    /// assert_eq!(record.size(), 10);
    /// assert_eq!(record.dt(), &Some(0));
    /// assert_eq!(record.dq(), &Some(1));
    /// assert_eq!(record.record_type(), &Kind::NonTerminating);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        size: usize,
        dt: Option<usize>,
        dq: Option<usize>,
        record_type: Kind,
    ) -> Result<Self> {
        match record_type {
            Kind::NonTerminating => {
                if dt.is_none() {
                    return Err(Error::InvalidNonTerminatingDt);
                }

                if dq.is_none() {
                    return Err(Error::InvalidNonTerminatingDq);
                }
            }
            Kind::Terminating => {
                if dt.is_some() {
                    return Err(Error::InvalidTerminatingDt);
                }

                if dq.is_some() {
                    return Err(Error::InvalidTerminatingDq);
                }
            }
        }

        Ok(Record(size, dt, dq, record_type))
    }

    /// Retuns the size of the ungapped alignment.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::data;
    ///
    /// let alignment: data::Record = "9\t1\t0".parse()?;
    /// assert_eq!(alignment.size(), 9);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn size(&self) -> usize {
        self.0
    }

    /// Returns the difference between this block and the next block for the
    /// reference/target sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::data;
    ///
    /// let alignment: data::Record = "9\t1\t0".parse()?;
    /// assert_eq!(*alignment.dt(), Some(1));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn dt(&self) -> &Option<usize> {
        &self.1
    }

    /// Returns the difference between this block and the next block for the
    /// query sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::data;
    ///
    /// let alignment: data::Record = "9\t1\t0".parse()?;
    /// assert_eq!(*alignment.dq(), Some(0));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn dq(&self) -> &Option<usize> {
        &self.2
    }

    /// Returns the record type.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::data;
    /// use chainfile::alignment::section::data::record::Kind;
    ///
    /// let alignment: data::Record = "9\t1\t0".parse()?;
    /// assert_eq!(*alignment.record_type(), Kind::NonTerminating);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn record_type(&self) -> &Kind {
        &self.3
    }
}

impl FromStr for Record {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts = s.split(ALIGNMENT_DATA_DELIMITER).collect::<Vec<_>>();
        let record_type = match parts.len() {
            NUM_ALIGNMENT_DATA_FIELDS_TERMINATING => Kind::Terminating,
            NUM_ALIGNMENT_DATA_FIELDS_NONTERMINATING => Kind::NonTerminating,
            _ => {
                return Err(Error::Parse(ParseError::IncorrectNumberOfFields(
                    parts.len(),
                )));
            }
        };

        let size = parts[0]
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidSize(err)))?;
        let (dt, dq) = match record_type {
            Kind::NonTerminating => {
                let dt = parts[1]
                    .parse()
                    .map_err(|err| Error::Parse(ParseError::InvalidDt(err)))?;
                let dq = parts[2]
                    .parse()
                    .map_err(|err| Error::Parse(ParseError::InvalidDq(err)))?;
                (Some(dt), Some(dq))
            }
            Kind::Terminating => (None, None),
        };

        Record::try_new(size, dt, dq, record_type)
    }
}

impl std::fmt::Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.size())?;

        match self.record_type() {
            Kind::NonTerminating => {
                write!(
                    f,
                    "{}{}{}{}",
                    ALIGNMENT_DATA_DELIMITER,
                    self.dt().expect("non-terminating record must have a dt"),
                    ALIGNMENT_DATA_DELIMITER,
                    self.dq().expect("non-terminating record must have a dq"),
                )?;
            }
            Kind::Terminating => {}
        }

        Ok(())
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_nonterminating_alignment_data() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let record = "9\t1\t0".parse::<Record>()?;

        assert_eq!(record.size(), 9);
        assert_eq!(*record.dt(), Some(1));
        assert_eq!(*record.dq(), Some(0));
        assert_eq!(*record.record_type(), Kind::NonTerminating);

        Ok(())
    }

    #[test]
    fn test_terminating_alignment_data() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let record = "9".parse::<Record>()?;

        assert_eq!(record.size(), 9);
        assert_eq!(*record.dt(), None);
        assert_eq!(*record.dq(), None);
        assert_eq!(*record.record_type(), Kind::Terminating);

        Ok(())
    }

    #[test]
    fn test_invalid_number_of_fields() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = "9\t0".parse::<Record>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid number of fields in alignment data: expected 3 \
             (non-terminating) or 1 (terminating) fields, found 2 fields"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_size() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = "?\t0\t1".parse::<Record>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid size: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_dt() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = "9\t?\t1".parse::<Record>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid dt: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_dq() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = "9\t0\t?".parse::<Record>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid dq: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_nonterminating_dt() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Record::try_new(9, None, Some(1), Kind::NonTerminating).unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected value for dt in non-terminating alignment data line, found no value"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_nonterminating_dq() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Record::try_new(9, Some(0), None, Kind::NonTerminating).unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected value for dq in non-terminating alignment data line, found no value"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_terminating_dt() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Record::try_new(9, Some(0), None, Kind::Terminating).unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected no value for dt in terminating alignment data line, found value"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_terminating_dq() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Record::try_new(9, None, Some(1), Kind::Terminating).unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected no value for dq in terminating alignment data line, found value"
        );

        Ok(())
    }

    #[test]
    fn test_nonterminating_alignment_data_display()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let alignment = Record::try_new(9, Some(1), Some(0), Kind::NonTerminating)?;

        assert_eq!(alignment.to_string(), "9\t1\t0");

        Ok(())
    }

    #[test]
    fn test_terminating_alignment_data_display()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let alignment = Record::try_new(9, None, None, Kind::Terminating)?;

        assert_eq!(alignment.to_string(), "9");

        Ok(())
    }
}