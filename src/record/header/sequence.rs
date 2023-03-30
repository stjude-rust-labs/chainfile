//! A sequence within a header record.

use std::num::ParseIntError;

use crate::record::header::strand::ParseStrandError;
use crate::record::header::strand::Strand;
use crate::record::header::HEADER_DELIMITER;

/// Errors associated with parsing a sequence.
#[derive(Debug)]
pub enum ParseError {
    /// An invalid chromosome size.
    InvalidChromosomeSize(ParseIntError),
    /// An invalid strand.
    InvalidStrand(ParseStrandError),
    /// An invalid alignment start.
    InvalidAlignmentStart(ParseIntError),
    /// An invalid alignment end.
    InvalidAlignmentEnd(ParseIntError),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::InvalidChromosomeSize(err) => write!(f, "invalid chromosome size: {}", err),
            ParseError::InvalidStrand(err) => write!(f, "invalid strand: {}", err),
            ParseError::InvalidAlignmentStart(err) => write!(f, "invalid alignment start: {}", err),
            ParseError::InvalidAlignmentEnd(err) => write!(f, "invalid alignment end: {}", err),
        }
    }
}

impl std::error::Error for ParseError {}

/// The sequence portion(s) of a header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sequence {
    /// The chromosome name.
    chromosome_name: String,
    /// The chromosome size.
    chromosome_size: usize,
    /// The strand.
    strand: Strand,
    /// The start of the alignment.
    alignment_start: usize,
    /// The end of the alignment.
    alignment_end: usize,
}

impl Sequence {
    /// Creates a new sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let sequence = chain::record::header::Sequence::new("seq0", "2", "+", "0", "2");
    /// ```
    pub fn new(
        chromosome_name: &str,
        chromosome_size: &str,
        strand: &str,
        alignment_start: &str,
        alignment_end: &str,
    ) -> Result<Self, ParseError> {
        Ok(Self {
            chromosome_name: chromosome_name.into(),
            chromosome_size: chromosome_size
                .parse()
                .map_err(ParseError::InvalidChromosomeSize)?,
            strand: strand.parse().map_err(ParseError::InvalidStrand)?,
            alignment_start: alignment_start
                .parse()
                .map_err(ParseError::InvalidAlignmentStart)?,
            alignment_end: alignment_end
                .parse()
                .map_err(ParseError::InvalidAlignmentEnd)?,
        })
    }

    /// Returns the chromosome name.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let sequence = chain::record::header::Sequence::new("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.chromosome_name(), "seq0");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn chromosome_name(&self) -> &String {
        &self.chromosome_name
    }

    /// Returns the chromosome size.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let sequence = chain::record::header::Sequence::new("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.chromosome_size(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn chromosome_size(&self) -> usize {
        self.chromosome_size
    }

    /// Returns the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::record::header::strand::Strand;
    ///
    /// let sequence = chain::record::header::Sequence::new("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.strand(), &Strand::Positive);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> &Strand {
        &self.strand
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let sequence = chain::record::header::Sequence::new("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.alignment_start(), 0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alignment_start(&self) -> usize {
        self.alignment_start
    }

    /// Returns the alignment end.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    ///
    /// let sequence = chain::record::header::Sequence::new("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.alignment_end(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alignment_end(&self) -> usize {
        self.alignment_end
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let chromosome_name = self.chromosome_name.to_string();
        let chromosome_size = self.chromosome_size.to_string();
        let strand = self.strand.to_string();
        let alignment_start = self.alignment_start.to_string();
        let alignment_end = self.alignment_end.to_string();

        let parts = [
            chromosome_name.as_str(),
            chromosome_size.as_str(),
            strand.as_str(),
            alignment_start.as_str(),
            alignment_end.as_str(),
        ];

        write!(f, "{}", parts.join(HEADER_DELIMITER.to_string().as_str()))
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_sequence_creation() -> Result<(), Box<dyn std::error::Error>> {
        let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
        assert_eq!(sequence.chromosome_name(), "seq0");
        assert_eq!(sequence.chromosome_size(), 2);
        assert_eq!(sequence.strand(), &Strand::Positive);
        assert_eq!(sequence.alignment_start(), 0);
        assert_eq!(sequence.alignment_end(), 2);
        Ok(())
    }

    #[test]
    fn test_invalid_chromosome_size() -> Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "A", "+", "0", "2").unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid chromosome size: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_strand() -> Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "2", "?", "0", "2").unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid strand: parse strand error: ? is not a valid strand"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_alignment_start() -> Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "2", "+", "?", "2").unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid alignment start: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_alignment_end() -> Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "2", "+", "0", "?").unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid alignment end: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_sequence_display() -> Result<(), Box<dyn std::error::Error>> {
        let sequence = Sequence::new("seq0", "2", "+", "0", "1")?;
        assert_eq!(sequence.to_string(), "seq0 2 + 0 1");
        Ok(())
    }
}
