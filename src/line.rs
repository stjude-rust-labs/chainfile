//! A line within a chain file.

use crate::alignment::section::data;
use crate::alignment::section::data::Record as AlignmentDataRecord;
use crate::alignment::section::header;
use crate::alignment::section::header::HEADER_PREFIX;
use crate::alignment::section::header::Record as HeaderRecord;

/// An error associated with parsing the chain file.
#[derive(Debug)]
pub enum Error {
    /// An invalid header record.
    InvalidHeaderRecord {
        /// The inner error.
        inner: header::Error,

        /// The literal line.
        line: String,
    },

    /// An invalid alignment data record.
    InvalidAlignmentDataRecord {
        /// The inner error.
        inner: data::Error,

        /// The literal line.
        line: String,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidHeaderRecord { inner, line } => {
                write!(f, "invalid header record: {}\n\nline: `{}`", inner, line)
            }
            Error::InvalidAlignmentDataRecord { inner, line } => {
                write!(
                    f,
                    "invalid alignment data record: {}\n\nline: `{}`",
                    inner, line
                )
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
    Header(HeaderRecord),

    /// An alignment data line.
    AlignmentData(AlignmentDataRecord),
}

impl Line {
    /// Attempts to return a reference to the inner header record.
    pub fn as_header(&self) -> Option<&HeaderRecord> {
        match self {
            Line::Header(record) => Some(record),
            _ => None,
        }
    }

    /// Consumes `self` and attempts to return a reference to the inner header
    /// record.
    pub fn into_header(self) -> Option<HeaderRecord> {
        match self {
            Line::Header(record) => Some(record),
            _ => None,
        }
    }

    /// Attempts to return a reference to the inner alignment data record.
    pub fn as_alignment_data(&self) -> Option<&AlignmentDataRecord> {
        match self {
            Line::AlignmentData(record) => Some(record),
            _ => None,
        }
    }

    /// Consumes `self` and attempts to return a reference to the inner
    /// alignment data record.
    pub fn into_alignment_data_record(self) -> Option<AlignmentDataRecord> {
        match self {
            Line::AlignmentData(record) => Some(record),
            _ => None,
        }
    }
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
        } else if s.starts_with(HEADER_PREFIX) {
            s.parse::<HeaderRecord>()
                .map(Line::Header)
                .map_err(|err| Error::InvalidHeaderRecord {
                    inner: err,
                    line: s.into(),
                })
        } else {
            s.parse::<AlignmentDataRecord>()
                .map(Line::AlignmentData)
                .map_err(|err| Error::InvalidAlignmentDataRecord {
                    inner: err,
                    line: s.into(),
                })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::section::data::record::Kind;

    #[test]
    pub fn valid_header_line() {
        let line = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<Line>()
            .unwrap();

        let record = line.into_header().unwrap();

        assert_eq!(record.score(), 0);
        assert_eq!(record.id(), 1);
    }

    #[test]
    pub fn valid_nonterminating_alignment_data_line() {
        let line = "9\t0\t1".parse::<Line>().unwrap();
        let record = line.into_alignment_data_record().unwrap();

        assert_eq!(record.size(), 9);
        assert_eq!(record.dt().unwrap(), 0);
        assert_eq!(record.dq().unwrap(), 1);
        assert_eq!(record.kind(), Kind::NonTerminating);
    }

    #[test]
    pub fn valid_terminating_alignment_data_line() {
        let line = "9".parse::<Line>().unwrap();
        let record = line.into_alignment_data_record().unwrap();

        assert_eq!(record.size(), 9);
        assert!(record.dt().is_none());
        assert!(record.dq().is_none());
        assert_eq!(record.kind(), Kind::Terminating);
    }

    #[test]
    pub fn invalid_header_line() {
        let err = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?"
            .parse::<Line>()
            .unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid header record: parse error: invalid id: invalid digit found in \
             string\n\nline: `chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?`"
        );
    }

    #[test]
    pub fn invalid_alignment_data_line() {
        let err = "9\t1".parse::<Line>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid alignment data record: parse error: invalid number of fields in alignment \
             data: expected 3 (non-terminating) or 1 (terminating) fields, found 2 \
             fields\n\nline: `9\t1`"
        );
    }
}
