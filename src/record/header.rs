//! A header record.

pub mod sequence;

use std::num::ParseIntError;
use std::str::FromStr;

pub use sequence::Sequence;

/// The prefix for a header record.
pub const HEADER_PREFIX: &str = "chain";

/// The delimiter for a header record.
pub const HEADER_DELIMITER: char = ' ';

/// The number of expected fields in a header record.
pub const NUM_HEADER_FIELDS: usize = 13;

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
    InvalidReferenceSequence(sequence::ParseError),
    /// An invalid query sequence.
    InvalidQuerySequence(sequence::ParseError),
    /// An invalid id.
    InvalidId(ParseIntError),
    /// The end position exceeds the size of the chromosome.
    EndPositionExceedsSize(String, usize, usize),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::IncorrectNumberOfFields(n) => write!(
                f,
                "invalid number of fields in header: expected {} fields, found {} fields",
                NUM_HEADER_FIELDS, n
            ),
            ParseError::InvalidPrefix(s) => {
                write!(
                    f,
                    "invalid prefix: expected \"{}\", found \"{}\"",
                    HEADER_PREFIX, s
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
                "the end position ({}) exceeds the size of the chromosome {} ({})",
                pos, chrom, size
            ),
        }
    }
}

impl std::error::Error for ParseError {}

/// A header record within a chain file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HeaderRecord {
    /// The chain score.
    score: usize,
    /// The reference sequence.
    reference_sequence: Sequence,
    /// The query sequence.
    query_sequence: Sequence,
    /// The chain id.
    id: usize,
}

impl HeaderRecord {
    /// Returns the score for the chain header line.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::record::HeaderRecord;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<HeaderRecord>()?;
    ///
    /// assert_eq!(header.score(), 0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn score(&self) -> usize {
        self.score
    }

    /// Returns the reference sequence for the chain header line.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Strand;
    /// use chain::record::HeaderRecord;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<HeaderRecord>()?;
    ///
    /// assert_eq!(header.reference_sequence().chromosome_name(), "seq0");
    /// assert_eq!(header.reference_sequence().chromosome_size(), 2);
    /// assert_eq!(header.reference_sequence().strand(), &Strand::Positive);
    /// assert_eq!(header.reference_sequence().alignment_start(), 0);
    /// assert_eq!(header.reference_sequence().alignment_end(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequence(&self) -> &Sequence {
        &self.reference_sequence
    }

    /// Returns the query sequence for the chain header line.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Strand;
    /// use chain::record::HeaderRecord;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<HeaderRecord>()?;
    ///
    /// assert_eq!(header.query_sequence().chromosome_name(), "seq0");
    /// assert_eq!(header.query_sequence().chromosome_size(), 2);
    /// assert_eq!(header.query_sequence().strand(), &Strand::Negative);
    /// assert_eq!(header.query_sequence().alignment_start(), 0);
    /// assert_eq!(header.query_sequence().alignment_end(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query_sequence(&self) -> &Sequence {
        &self.query_sequence
    }

    /// Returns the id for the chain header line.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Strand;
    /// use chain::record::HeaderRecord;
    ///
    /// let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<HeaderRecord>()?;
    ///
    /// assert_eq!(header.id(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn id(&self) -> usize {
        self.id
    }
}

impl FromStr for HeaderRecord {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(HEADER_DELIMITER).collect::<Vec<_>>();
        if parts.len() != NUM_HEADER_FIELDS {
            return Err(ParseError::IncorrectNumberOfFields(parts.len()));
        }

        let chain = parts[0];
        if chain != HEADER_PREFIX {
            return Err(ParseError::InvalidPrefix(chain.into()));
        }

        let score = parts[1].parse().map_err(ParseError::InvalidScore)?;
        let reference_sequence = Sequence::new(parts[2], parts[3], parts[4], parts[5], parts[6])
            .map_err(ParseError::InvalidReferenceSequence)?;
        let query_sequence = Sequence::new(parts[7], parts[8], parts[9], parts[10], parts[11])
            .map_err(ParseError::InvalidQuerySequence)?;
        let id = parts[12].parse().map_err(ParseError::InvalidId)?;

        if reference_sequence.chromosome_size() < reference_sequence.alignment_end() {
            return Err(ParseError::EndPositionExceedsSize(
                reference_sequence.chromosome_name().clone(),
                reference_sequence.alignment_end(),
                reference_sequence.chromosome_size(),
            ));
        }

        if query_sequence.chromosome_size() < query_sequence.alignment_end() {
            return Err(ParseError::EndPositionExceedsSize(
                query_sequence.chromosome_name().clone(),
                query_sequence.alignment_end(),
                query_sequence.chromosome_size(),
            ));
        }

        Ok(HeaderRecord {
            score,
            reference_sequence,
            query_sequence,
            id,
        })
    }
}

impl std::fmt::Display for HeaderRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let score = self.score.to_string();
        let reference_sequence = self.reference_sequence.to_string();
        let query_sequence = self.query_sequence.to_string();
        let id = self.id.to_string();

        let parts = [
            HEADER_PREFIX,
            score.as_str(),
            reference_sequence.as_str(),
            query_sequence.as_str(),
            id.as_str(),
        ];

        write!(f, "{}", parts.join(HEADER_DELIMITER.to_string().as_str()))
    }
}

#[cfg(test)]
pub mod tests {
    use crate::core::Strand;

    use super::*;

    #[test]
    pub fn test_parsing_header_record() -> Result<(), Box<dyn std::error::Error>> {
        let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<HeaderRecord>()?;

        assert_eq!(header.score(), 0);

        assert_eq!(header.reference_sequence().chromosome_name(), "seq0");
        assert_eq!(header.reference_sequence().chromosome_size(), 2);
        assert_eq!(header.reference_sequence().strand(), &Strand::Positive);
        assert_eq!(header.reference_sequence().alignment_start(), 0);
        assert_eq!(header.reference_sequence().alignment_end(), 2);

        assert_eq!(header.query_sequence().chromosome_name(), "seq0");
        assert_eq!(header.query_sequence().chromosome_size(), 2);
        assert_eq!(header.query_sequence().strand(), &Strand::Negative);
        assert_eq!(header.query_sequence().alignment_start(), 0);
        assert_eq!(header.query_sequence().alignment_end(), 2);

        assert_eq!(header.id(), 1);

        Ok(())
    }

    #[test]
    fn test_invalid_number_of_fields() -> Result<(), Box<dyn std::error::Error>> {
        let err = "foo 0 seq0 2 + 0 2 seq0 2 - 0 2"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid number of fields in header: expected 13 fields, found 12 fields"
        );
        Ok(())
    }

    #[test]
    fn test_invalid_prefix() -> Result<(), Box<dyn std::error::Error>> {
        let err = "foo 0 seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid prefix: expected \"chain\", found \"foo\""
        );
        Ok(())
    }

    #[test]
    fn test_invalid_score() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain ? seq0 2 + 0 2 seq0 2 - 0 2 1"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid score: invalid digit found in string"
        );
        Ok(())
    }

    #[test]
    fn test_invalid_reference_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain 0 seq0 ? + 0 2 seq0 2 - 0 2 1"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid reference sequence: invalid chromosome size: invalid digit found in string"
        );
        Ok(())
    }

    #[test]
    fn test_invalid_query_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain 0 seq0 2 + 0 2 seq0 ? - 0 2 1"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "invalid query sequence: invalid chromosome size: invalid digit found in string"
        );
        Ok(())
    }

    #[test]
    fn test_invalid_id() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 ?"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(err.to_string(), "invalid id: invalid digit found in string");
        Ok(())
    }

    #[test]
    fn test_end_is_greater_than_size_reference() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain 0 seq0 2 + 0 3 seq0 2 - 0 1 1"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "the end position (3) exceeds the size of the chromosome seq0 (2)"
        );
        Ok(())
    }

    #[test]
    fn test_end_is_greater_than_size_query() -> Result<(), Box<dyn std::error::Error>> {
        let err = "chain 0 seq0 2 + 0 1 seq0 2 - 0 3 1"
            .parse::<HeaderRecord>()
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "the end position (3) exceeds the size of the chromosome seq0 (2)"
        );
        Ok(())
    }

    #[test]
    pub fn test_header_record_display() -> Result<(), Box<dyn std::error::Error>> {
        let header = "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1".parse::<HeaderRecord>()?;
        assert_eq!(header.to_string(), "chain 0 seq0 2 + 0 2 seq0 2 - 0 2 1");
        Ok(())
    }
}
