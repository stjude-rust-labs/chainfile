//! A line within a chain file.

use std::str::FromStr;

use crate::record::alignment_data;
use crate::record::alignment_data::AlignmentDataRecord;
use crate::record::header;
use crate::record::header::HeaderRecord;
use crate::record::header::HEADER_PREFIX;

/// An error associated with parsing the chain file.
#[derive(Debug)]
pub enum ParseError {
    /// An invalid header record.
    InvalidHeaderRecord(header::ParseError, String),
    /// An invalid alignment data record.
    InvalidAlignmentDataRecord(alignment_data::ParseError, String),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InvalidHeaderRecord(err, line) => {
                write!(f, "invalid header record: {}\n\nline: {}", err, line)
            }
            ParseError::InvalidAlignmentDataRecord(err, line) => {
                write!(
                    f,
                    "invalid alignment data record: {}\n\nline: {}",
                    err, line
                )
            }
        }
    }
}

impl std::error::Error for ParseError {}

/// A line within a chain file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Line {
    /// An empty line.
    Empty,
    /// A header line.
    Header(HeaderRecord),
    /// An alignment data line.
    AlignmentData(AlignmentDataRecord),
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

impl FromStr for Line {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Ok(Self::Empty)
        } else if s.starts_with(HEADER_PREFIX) {
            s.parse::<HeaderRecord>()
                .map(Line::Header)
                .map_err(|e| ParseError::InvalidHeaderRecord(e, s.into()))
        } else {
            s.parse::<AlignmentDataRecord>()
                .map(Line::AlignmentData)
                .map_err(|e| ParseError::InvalidAlignmentDataRecord(e, s.into()))
        }
    }
}

#[cfg(test)]
pub mod tests {
    use crate::record::alignment_data::AlignmentDataRecordType;

    use super::*;

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
            assert_eq!(
                *record.record_type(),
                AlignmentDataRecordType::NonTerminating
            );
        }
        Ok(())
    }

    #[test]
    pub fn test_valid_terminating_alignment_data_line() -> Result<(), Box<dyn std::error::Error>> {
        let line = "9".parse::<Line>()?;
        assert!(matches!(line, Line::AlignmentData(_)));
        if let Line::AlignmentData(record) = line {
            assert_eq!(*record.record_type(), AlignmentDataRecordType::Terminating);
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
            "invalid header record: invalid id: invalid digit \
             found in string\n\nline: chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?"
        );
        Ok(())
    }

    #[test]
    pub fn test_invalid_alignment_data_line() -> Result<(), Box<dyn std::error::Error>> {
        let err = "9\t1".parse::<Line>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid alignment data record: invalid number of fields in alignment data: \
            expected 3 (non-terminating) or 1 (terminating) fields, found 2 fields\n\n\
            line: 9\t1"
        );
        Ok(())
    }
}
