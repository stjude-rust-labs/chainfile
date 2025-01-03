//! A sequence within a header record.

use std::num::ParseIntError;

use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval;
use omics::coordinate::interval::interbase::Interval;
use omics::coordinate::position::Number;
use omics::coordinate::strand;

use crate::alignment::section::header::DELIMITER;

////////////////////////////////////////////////////////////////////////////////////////
// Errors
////////////////////////////////////////////////////////////////////////////////////////

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
            ParseError::InvalidChromosomeSize(err) => write!(f, "invalid chromosome size: {err}"),
            ParseError::InvalidStrand(err) => write!(f, "invalid strand: {err}"),
            ParseError::InvalidAlignmentStart(err) => write!(f, "invalid alignment start: {err}"),
            ParseError::InvalidAlignmentEnd(err) => write!(f, "invalid alignment end: {err}"),
        }
    }
}

impl std::error::Error for ParseError {}

/// An error related to a sequence.
#[derive(Debug)]
pub enum Error {
    /// An interval error.
    Interval(interval::Error),

    /// A parse error.
    Parse(ParseError),

    /// The start position is greater than the end position.
    StartPositionGreaterThanEndPosition,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Interval(err) => write!(f, "interval error: {err}"),
            Error::Parse(err) => write!(f, "parse error: {err}"),
            Error::StartPositionGreaterThanEndPosition => {
                write!(f, "start position greater than end position")
            }
        }
    }
}

impl std::error::Error for Error {}

/// A [`Result`](std::result::Result) with an [`Error`].
type Result<T> = std::result::Result<T, Error>;

////////////////////////////////////////////////////////////////////////////////////////
// Sequence
////////////////////////////////////////////////////////////////////////////////////////

/// The sequence portion of a header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Sequence {
    /// The chromosome name.
    chromosome_name: Contig,

    /// The chromosome size.
    chromosome_size: Number,

    /// The strand.
    strand: Strand,

    /// The start of the alignment.
    alignment_start: Number,

    /// The end of the alignment.
    alignment_end: Number,
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
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    ///
    /// assert_eq!(sequence.chromosome_name(), "seq0");
    /// assert_eq!(sequence.chromosome_size(), 2);
    /// assert_eq!(sequence.strand(), Strand::Positive);
    /// assert_eq!(sequence.alignment_start(), 0);
    /// assert_eq!(sequence.alignment_end(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_from_str_parts(
        chromosome_name: &str,
        chromosome_size: &str,
        strand: &str,
        alignment_start: &str,
        alignment_end: &str,
    ) -> Result<Self> {
        let chromosome_name = chromosome_name.into();

        let chromosome_size = chromosome_size
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidChromosomeSize(err)))?;

        let strand = strand
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidStrand(err)))?;

        let alignment_start = alignment_start
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidAlignmentStart(err)))?;

        let alignment_end = alignment_end
            .parse()
            .map_err(|err| Error::Parse(ParseError::InvalidAlignmentEnd(err)))?;

        if alignment_start > alignment_end {
            return Err(Error::StartPositionGreaterThanEndPosition);
        }

        Ok(Self {
            chromosome_name,
            chromosome_size,
            strand,
            alignment_start,
            alignment_end,
        })
    }

    /// Returns the chromosome name.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    /// assert_eq!(sequence.chromosome_name(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn chromosome_name(&self) -> &str {
        self.chromosome_name.as_ref()
    }

    /// Returns the chromosome size.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    /// assert_eq!(sequence.chromosome_size(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn chromosome_size(&self) -> Number {
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
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    /// assert_eq!(sequence.strand(), Strand::Positive);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    /// assert_eq!(sequence.alignment_start(), 0);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alignment_start(&self) -> Number {
        self.alignment_start
    }

    /// Returns the alignment end.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header::Sequence;
    ///
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    /// assert_eq!(sequence.alignment_end(), 2);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn alignment_end(&self) -> Number {
        self.alignment_end
    }

    /// Gets the sequence as an interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::alignment::section::header::Sequence;
    /// use omics::coordinate::Strand;
    ///
    /// // Positive strand.
    ///
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2")?;
    /// let interval = sequence.interval()?;
    ///
    /// assert_eq!(interval.start().contig().as_str(), "seq0");
    /// assert_eq!(interval.start().strand(), Strand::Positive);
    /// assert_eq!(interval.start().position().get(), 0);
    ///
    /// assert_eq!(interval.end().contig().as_str(), "seq0");
    /// assert_eq!(interval.end().strand(), Strand::Positive);
    /// assert_eq!(interval.end().position().get(), 2);
    ///
    /// // Negative strand.
    ///
    /// let sequence = Sequence::try_from_str_parts("seq0", "2", "-", "0", "2")?;
    /// let interval = sequence.interval()?;
    ///
    /// assert_eq!(interval.start().contig().as_str(), "seq0");
    /// assert_eq!(interval.start().strand(), Strand::Negative);
    /// assert_eq!(interval.start().position().get(), 2);
    ///
    /// assert_eq!(interval.end().contig().as_str(), "seq0");
    /// assert_eq!(interval.end().strand(), Strand::Negative);
    /// assert_eq!(interval.end().position().get(), 0);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn interval(&self) -> Result<Interval> {
        // NOTE: this is necessary because chainfiles store the positions of
        // sequences on the negative strand as their reverse complement.
        let (start_pos, end_pos) = match self.strand {
            Strand::Positive => (
                self.alignment_start() as Number,
                self.alignment_end() as Number,
            ),
            Strand::Negative => (
                // NOTE: coordinates on the negative strand are stored as the
                // reverse complement of the real sequence.
                self.chromosome_size - self.alignment_start() as Number,
                self.chromosome_size - self.alignment_end() as Number,
            ),
        };

        let start = Coordinate::new(self.chromosome_name(), self.strand(), start_pos);
        let end = Coordinate::new(self.chromosome_name(), self.strand(), end_pos);

        Interval::try_new(start, end).map_err(Error::Interval)
    }
}

impl std::fmt::Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}{}{}{}{}{}{}{}",
            self.chromosome_name.as_str(),
            DELIMITER,
            self.chromosome_size,
            DELIMITER,
            self.strand,
            DELIMITER,
            self.alignment_start,
            DELIMITER,
            self.alignment_end
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_new() {
        let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "2").unwrap();

        assert_eq!(sequence.chromosome_name(), "seq0");
        assert_eq!(sequence.chromosome_size(), 2);
        assert_eq!(sequence.strand(), Strand::Positive);
        assert_eq!(sequence.alignment_start(), 0);
        assert_eq!(sequence.alignment_end(), 2);
    }

    #[test]
    fn invalid_chromosome_size() {
        let err = Sequence::try_from_str_parts("seq0", "A", "+", "0", "2").unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::InvalidChromosomeSize(_))
        ));
        assert_eq!(
            err.to_string(),
            "parse error: invalid chromosome size: invalid digit found in string"
        );
    }

    #[test]
    fn invalid_strand() {
        let err = Sequence::try_from_str_parts("seq0", "2", "?", "0", "2").unwrap_err();

        assert!(matches!(err, Error::Parse(ParseError::InvalidStrand(_))));
        assert_eq!(
            err.to_string(),
            "parse error: invalid strand: parse error: invalid strand: ?"
        );
    }

    #[test]
    fn invalid_alignment_start() {
        let err = Sequence::try_from_str_parts("seq0", "2", "+", "?", "2").unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::InvalidAlignmentStart(_))
        ));
        assert_eq!(
            err.to_string(),
            "parse error: invalid alignment start: invalid digit found in string"
        );
    }

    #[test]
    fn invalid_alignment_end() {
        let err = Sequence::try_from_str_parts("seq0", "2", "+", "0", "?").unwrap_err();

        assert!(matches!(
            err,
            Error::Parse(ParseError::InvalidAlignmentEnd(_))
        ));
        assert_eq!(
            err.to_string(),
            "parse error: invalid alignment end: invalid digit found in string"
        );
    }

    #[test]
    fn display() {
        let sequence = Sequence::try_from_str_parts("seq0", "2", "+", "0", "1").unwrap();
        assert_eq!(sequence.to_string(), "seq0 2 + 0 1");
    }
}
