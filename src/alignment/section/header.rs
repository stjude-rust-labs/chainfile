//! A header record.

pub mod sequence;

use std::num::ParseIntError;
use std::str::FromStr;

use omics::coordinate::position::Number;
pub use sequence::Sequence;

/// The prefix for a header record.
pub const HEADER_PREFIX: &str = "chain";

/// The delimiter for a header record.
pub const DELIMITER: char = ' ';

/// The number of expected fields in a header record.
pub const NUM_HEADER_FIELDS: usize = 13;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

/// An error associated with parsing a header record.
#[derive(Debug)]
pub enum ParseError {
    /// An incorrect number of fields in the header line.
    IncorrectNumberOfFields(usize),

    /// An invalid prefix.
    InvalidPrefix(String),

    /// An invalid score.
    InvalidScore(ParseIntError),

    /// An invalid reference sequence.
    InvalidReferenceSequence(sequence::Error),

    /// An invalid query sequence.
    InvalidQuerySequence(sequence::Error),

    /// An invalid id.
    InvalidId(ParseIntError),

    /// The end position exceeds the size of the chromosome.
    EndPositionExceedsSize(String, Number, Number),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::IncorrectNumberOfFields(fields) => write!(
                f,
                "invalid number of fields in header: expected {} fields, found {} fields",
                NUM_HEADER_FIELDS, fields
            ),
            ParseError::InvalidPrefix(prefix) => {
                write!(
                    f,
                    "invalid prefix: expected \"{}\", found \"{}\"",
                    HEADER_PREFIX, prefix
                )
            }
            ParseError::InvalidScore(err) => write!(f, "invalid score: {}", err),
            ParseError::InvalidReferenceSequence(err) => {
                write!(f, "invalid reference sequence: {}", err)
            }
            ParseError::InvalidQuerySequence(err) => write!(f, "invalid query sequence: {}", err),
            ParseError::InvalidId(err) => write!(f, "invalid id: {}", err),
            ParseError::EndPositionExceedsSize(chrom, pos, size) => write!(
                f,
                "the end position ({}) exceeds the size of the chromosome `{}` ({})",
                pos, chrom, size
            ),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a [`Record`].
#[derive(Debug)]
pub enum Error {
    /// A parse error.
    Parse(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Parse(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// Record
////////////////////////////////////////////////////////////////////////////////////////

/// A header record within a chain file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    /// The chain score.
    score: usize,

    /// The reference sequence.
    reference_sequence: Sequence,

    /// The query sequence.
    query_sequence: Sequence,

    /// The chain id.
    id: usize,
}

impl Record {
    /// Gets the score.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<header::Record>()?;
    ///
    /// assert_eq!(header.score(), 0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn score(&self) -> usize {
        self.score
    }

    /// Gets the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header;
    /// use omics::coordinate::Strand;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<header::Record>()?;
    ///
    /// assert_eq!(header.reference_sequence().chromosome_name(), "seq0");
    /// assert_eq!(header.reference_sequence().chromosome_size(), 2);
    /// assert_eq!(header.reference_sequence().strand(), Strand::Positive);
    /// assert_eq!(header.reference_sequence().alignment_start(), 0);
    /// assert_eq!(header.reference_sequence().alignment_end(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequence(&self) -> &Sequence {
        &self.reference_sequence
    }

    /// Gets the query sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header;
    /// use omics::coordinate::Strand;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<header::Record>()?;
    ///
    /// assert_eq!(header.query_sequence().chromosome_name(), "seq0");
    /// assert_eq!(header.query_sequence().chromosome_size(), 2);
    /// assert_eq!(header.query_sequence().strand(), Strand::Negative);
    /// assert_eq!(header.query_sequence().alignment_start(), 0);
    /// assert_eq!(header.query_sequence().alignment_end(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query_sequence(&self) -> &Sequence {
        &self.query_sequence
    }

    /// Gets the id.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header;
    /// use omics::coordinate::Strand;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<header::Record>()?;
    ///
    /// assert_eq!(header.id(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn id(&self) -> usize {
        self.id
    }
}

impl FromStr for Record {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts = s.split(DELIMITER).collect::<Vec<_>>();
        if parts.len() != NUM_HEADER_FIELDS {
            return Err(Error::Parse(ParseError::IncorrectNumberOfFields(
                parts.len(),
            )));
        }

        let chain = parts[0];
        if chain != HEADER_PREFIX {
            return Err(Error::Parse(ParseError::InvalidPrefix(chain.into())));
        }

        let score = parts[1]
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidScore(err)))?;
        let reference_sequence =
            Sequence::try_from_str_parts(parts[2], parts[3], parts[4], parts[5], parts[6])
                .map_err(|err| Error::Parse(ParseError::InvalidReferenceSequence(err)))?;
        let query_sequence =
            Sequence::try_from_str_parts(parts[7], parts[8], parts[9], parts[10], parts[11])
                .map_err(|err| Error::Parse(ParseError::InvalidQuerySequence(err)))?;
        let id = parts[12]
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidId(err)))?;

        if reference_sequence.chromosome_size() < reference_sequence.alignment_end() {
            return Err(Error::Parse(ParseError::EndPositionExceedsSize(
                reference_sequence.chromosome_name().to_string(),
                reference_sequence.alignment_end(),
                reference_sequence.chromosome_size(),
            )));
        }

        if query_sequence.chromosome_size() < query_sequence.alignment_end() {
            return Err(Error::Parse(ParseError::EndPositionExceedsSize(
                query_sequence.chromosome_name().to_string(),
                query_sequence.alignment_end(),
                query_sequence.chromosome_size(),
            )));
        }

        Ok(Record {
            score,
            reference_sequence,
            query_sequence,
            id,
        })
    }
}

impl std::fmt::Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}{}{}{}{}{}{}{}",
            HEADER_PREFIX,
            DELIMITER,
            self.score,
            DELIMITER,
            self.reference_sequence,
            DELIMITER,
            self.query_sequence,
            DELIMITER,
            self.id
        )
    }
}

#[cfg(test)]
mod tests {
    use omics::coordinate::Strand;

    use super::*;

    #[test]
    pub fn parse() {
        let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<Record>()
            .unwrap();

        assert_eq!(header.score(), 0);

        assert_eq!(header.reference_sequence().chromosome_name(), "seq0");
        assert_eq!(header.reference_sequence().chromosome_size(), 2);
        assert_eq!(header.reference_sequence().strand(), Strand::Positive);
        assert_eq!(header.reference_sequence().alignment_start(), 0);
        assert_eq!(header.reference_sequence().alignment_end(), 2);

        assert_eq!(header.query_sequence().chromosome_name(), "seq0");
        assert_eq!(header.query_sequence().chromosome_size(), 2);
        assert_eq!(header.query_sequence().strand(), Strand::Negative);
        assert_eq!(header.query_sequence().alignment_start(), 0);
        assert_eq!(header.query_sequence().alignment_end(), 2);

        assert_eq!(header.id(), 1);
    }

    #[test]
    fn incorrect_number_of_fields() {
        let err = "foo 0 seq0 2 + 0 2 seq0 2 - 0 2"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::IncorrectNumberOfFields(_))
        ));

        assert_eq!(
            err.to_string(),
            "parse error: invalid number of fields in header: expected 13 fields, found 12 fields"
        );
    }

    #[test]
    fn invalid_prefix() {
        let err = "foo 0 seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(err, Error::Parse(ParseError::InvalidPrefix(_))));
        assert_eq!(
            err.to_string(),
            "parse error: invalid prefix: expected \"chain\", found \"foo\""
        );
    }

    #[test]
    fn invalid_score() {
        let err = "chain ? seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(err, Error::Parse(ParseError::InvalidScore(_))));
        assert_eq!(
            err.to_string(),
            "parse error: invalid score: invalid digit found in string"
        );
    }

    #[test]
    fn invalid_reference_sequence() {
        let err = "chain 0 seq0 ? + 0 2 seq0 2 - 0 2 1"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::InvalidReferenceSequence(_))
        ));

        assert_eq!(
            err.to_string(),
            "parse error: invalid reference sequence: parse error: invalid chromosome size: \
             invalid digit found in string"
        );
    }

    #[test]
    fn invalid_query_sequence() {
        let err = "chain 0 seq0 2 + 0 2 seq0 ? - 0 2 1"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::InvalidQuerySequence(_))
        ));

        assert_eq!(
            err.to_string(),
            "parse error: invalid query sequence: parse error: invalid chromosome size: invalid \
             digit found in string"
        );
    }

    #[test]
    fn invalid_id() {
        let err = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(err, Error::Parse(ParseError::InvalidId(_))));
        assert_eq!(
            err.to_string(),
            "parse error: invalid id: invalid digit found in string"
        );
    }

    #[test]
    fn end_is_greater_than_size_reference() {
        let err = "chain 0 seq0 2 + 0 3 seq0 2 - 0 1 1"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::EndPositionExceedsSize(_, _, _))
        ));

        assert_eq!(
            err.to_string(),
            "parse error: the end position (3) exceeds the size of the chromosome `seq0` (2)"
        );
    }

    #[test]
    fn end_is_greater_than_size_query() {
        let err = "chain 0 seq0 2 + 0 1 seq0 2 - 0 3 1"
            .parse::<Record>()
            .unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::EndPositionExceedsSize(_, _, _))
        ));

        assert_eq!(
            err.to_string(),
            "parse error: the end position (3) exceeds the size of the chromosome `seq0` (2)"
        );
    }

    #[test]
    pub fn display() {
        let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<Record>()
            .unwrap();

        assert_eq!(header.to_string(), "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1");
    }
}
