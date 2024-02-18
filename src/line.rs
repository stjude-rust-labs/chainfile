//! A line within a chain file.

use crate::alignment;
use crate::alignment::section::data;
use crate::alignment::section::header;
use crate::alignment::section::header::PREFIX;

/// An error associated with parsing the chain file.
#[derive(Debug)]
pub enum Error {
    /// An invalid header record.
    InvalidHeaderRecord(header::Error, String),

    /// An invalid alignment data record.
    InvalidAlignmentDataRecord(data::Error, String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidHeaderRecord(err, line) => {
                write!(f, "invalid header record: {}: line: {}", err, line)
            }
            Error::InvalidAlignmentDataRecord(err, line) => {
                write!(f, "invalid alignment data record: {}: line: {}", err, line)
            }
        }
    }
}

impl std::error::Error for Error {}

/// A line within a chain file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Line {
    /// An empty line.
    Empty,

    /// A header line.
    Header(header::Record),

    /// An alignment data line.
    AlignmentData(alignment::section::data::Record),
}

impl std::fmt::Display for Line {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Line::Empty => write!(f, ""),
            Line::Header(record) => write!(f, "{}", record),
            Line::AlignmentData(record) => write!(f, "{}", record),
        }
    }
}

impl std::str::FromStr for Line {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Ok(Self::Empty)
        } else if s.starts_with(PREFIX) {
            s.parse::<header::Record>()
                .map(Line::Header)
                .map_err(|e| Error::InvalidHeaderRecord(e, s.into()))
        } else {
            s.parse::<alignment::section::data::Record>()
                .map(Line::AlignmentData)
                .map_err(|e| Error::InvalidAlignmentDataRecord(e, s.into()))
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::alignment::section::data::record::Kind;

    #[test]
    pub fn test_valid_header_line() -> Result<(), Box<dyn std::error::Error>> {
        let line = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<Line>()?;
        assert!(matches!(line, Line::Header(_)));
        Ok(())
    }

    #[test]
    pub fn test_valid_nonterminating_alignment_data_line() -> Result<(), Box<dyn std::error::Error>>
    {
        let line = "9\t0\t1".parse::<Line>()?;
        assert!(matches!(line, Line::AlignmentData(_)));
        if let Line::AlignmentData(record) = line {
            assert_eq!(*record.record_type(), Kind::NonTerminating);
        }
        Ok(())
    }

    #[test]
    pub fn test_valid_terminating_alignment_data_line() -> Result<(), Box<dyn std::error::Error>> {
        let line = "9".parse::<Line>()?;
        assert!(matches!(line, Line::AlignmentData(_)));
        if let Line::AlignmentData(record) = line {
            assert_eq!(*record.record_type(), Kind::Terminating);
        }
        Ok(())
    }

    #[test]
    pub fn test_invalid_header_line() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?"
            .parse::<Line>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid header record: parse error: invalid id: invalid digit found in string: line: \
             chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?"
        );
        Ok(())
    }

    #[test]
    pub fn test_invalid_alignment_data_line() -> Result<(), Box<dyn std::error::Error>> {
        let err = "9\t1".parse::<Line>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid alignment data record: parse error: invalid number of fields in alignment \
             data: expected 3 (non-terminating) or 1 (terminating) fields, found 2 fields: line: \
             9\t1"
        );
        Ok(())
    }
}
