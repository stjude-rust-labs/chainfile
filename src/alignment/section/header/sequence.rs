//! A sequence within a header record.

use std::num::ParseIntError;

use omics::coordinate;
use omics::coordinate::Strand;
use omics::coordinate::interval;
use omics::coordinate::interval::r#trait::Interval as _;
use omics::coordinate::interval::zero::Interval;
use omics::coordinate::strand;
use omics::coordinate::system::Zero;
use omics::coordinate::zero::Coordinate;

use crate::alignment::section::header::DELIMITER;

/// Errors associated with parsing a sequence.
#[derive(Debug)]
pub enum ParseError {
    /// An invalid chromosome size.
    InvalidChromosomeSize(ParseIntError),

    /// An invalid strand.
    InvalidStrand(strand::Error),

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

/// An error related to a [`Sequence`].
#[derive(Debug)]
pub enum Error {
    /// A coordinate error.
    Coordinate(coordinate::Error),

    /// An interval error.
    Interval(interval::Error),

    /// A parse error.
    Parse(ParseError),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Coordinate(err) => write!(f, "coordinate error: {err}"),
            Error::Interval(err) => write!(f, "interval error: {err}"),
            Error::Parse(err) => write!(f, "parse error: {err}"),
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

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
    /// use chainfile::alignment::section::header::Sequence;
    /// use omics::coordinate::Strand;
    ///
    /// let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.chromosome_name(), "seq0");
    /// assert_eq!(sequence.chromosome_size(), 2);
    /// assert_eq!(sequence.strand(), &Strand::Positive);
    /// assert_eq!(sequence.alignment_start(), 0);
    /// assert_eq!(sequence.alignment_end(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(
        chromosome_name: &str,
        chromosome_size: &str,
        strand: &str,
        alignment_start: &str,
        alignment_end: &str,
    ) -> Result<Self> {
        Ok(Self {
            chromosome_name: chromosome_name.into(),
            chromosome_size: chromosome_size
                .parse()
                .map_err(|err| Error::Parse(ParseError::InvalidChromosomeSize(err)))?,
            strand: strand
                .parse()
                .map_err(|err| Error::Parse(ParseError::InvalidStrand(err)))?,
            alignment_start: alignment_start
                .parse()
                .map_err(|err| Error::Parse(ParseError::InvalidAlignmentStart(err)))?,
            alignment_end: alignment_end
                .parse()
                .map_err(|err| Error::Parse(ParseError::InvalidAlignmentEnd(err)))?,
        })
    }

    /// Returns the chromosome name.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
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
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
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
    /// use chainfile::alignment::section::header::Sequence;
    /// use omics::coordinate::Strand;
    ///
    /// let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
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
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
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
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
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

        write!(f, "{}", parts.join(DELIMITER.to_string().as_str()))
    }
}

impl TryFrom<Sequence> for omics::coordinate::Interval<Zero> {
    type Error = Error;

    fn try_from(value: Sequence) -> Result<Self> {
        let start = Coordinate::try_new(
            value.chromosome_name().clone(),
            Strand::Positive,
            value.alignment_start(),
        )
        .map_err(Error::Coordinate)?;

        let end = Coordinate::try_new(
            value.chromosome_name().clone(),
            Strand::Positive,
            value.alignment_end(),
        )
        .map_err(Error::Coordinate)?;

        let interval = Interval::try_new(start, end).map_err(Error::Interval)?;

        match value.strand() {
            Strand::Positive => Ok(interval),
            Strand::Negative => {
                // SAFETY: for all of the sequences read from chainfiles, this
                // will always unwrap.
                Ok(interval.complement().map_err(Error::Interval)?.unwrap())
            }
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_sequence_creation() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let sequence = Sequence::new("seq0", "2", "+", "0", "2")?;
        assert_eq!(sequence.chromosome_name(), "seq0");
        assert_eq!(sequence.chromosome_size(), 2);
        assert_eq!(sequence.strand(), &Strand::Positive);
        assert_eq!(sequence.alignment_start(), 0);
        assert_eq!(sequence.alignment_end(), 2);
        Ok(())
    }

    #[test]
    fn test_invalid_chromosome_size() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "A", "+", "0", "2").unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid chromosome size: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_strand() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "2", "?", "0", "2").unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid strand: parse error: invalid value for strand: ?"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_alignment_start() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "2", "+", "?", "2").unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid alignment start: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_alignment_end() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let err = Sequence::new("seq0", "2", "+", "0", "?").unwrap_err();

        assert_eq!(
            err.to_string(),
            "parse error: invalid alignment end: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_sequence_display() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let sequence = Sequence::new("seq0", "2", "+", "0", "1")?;
        assert_eq!(sequence.to_string(), "seq0 2 + 0 1");
        Ok(())
    }
}
