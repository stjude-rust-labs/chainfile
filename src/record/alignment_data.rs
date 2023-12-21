//! An alignment data record.

use std::num::ParseIntError;
use std::str::FromStr;

/// The delimiter for an alignment data record.
const ALIGNMENT_DATA_DELIMITER: char = '\t';

/// The number of expected fields in a non-terminating alignment data record.
pub const NUM_ALIGNMENT_DATA_FIELDS_NONTERMINATING: usize = 3;

/// The number of expected fields in a terminating alignment data record.
pub const NUM_ALIGNMENT_DATA_FIELDS_TERMINATING: usize = 1;

/// A type of alignment data record (whether the alignment data record is
/// terminating or not).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AlignmentDataRecordType {
    /// Every line except the last line in an alignment section.
    NonTerminating,
    /// The last line in an alignment section.
    Terminating,
}

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
    /// An invalid dt value for a non-terminating alignment data line.
    InvalidNonTerminatingDt,
    /// An invalid dq value for a non-terminating alignment data line.
    InvalidNonTerminatingDq,
    /// An invalid dt value for a terminating alignment data line.
    InvalidTerminatingDt,
    /// An invalid dq value for a terminating alignment data line.
    InvalidTerminatingDq,
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
            ParseError::InvalidNonTerminatingDt => write!(
                f,
                "expected value for dt in non-terminating alignment data line, found no value"
            ),
            ParseError::InvalidNonTerminatingDq => write!(
                f,
                "expected value for dq in non-terminating alignment data line, found no value"
            ),
            ParseError::InvalidTerminatingDt => write!(
                f,
                "expected no value for dt in terminating alignment data line, found value"
            ),
            ParseError::InvalidTerminatingDq => write!(
                f,
                "expected no value for dq in terminating alignment data line, found value"
            ),
        }
    }
}

impl std::error::Error for ParseError {}

/// An alignment data record within a chain file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentDataRecord(usize, Option<usize>, Option<usize>, AlignmentDataRecordType);

impl AlignmentDataRecord {
    /// Attempts to create a new [`AlignmentDataRecord`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::record::alignment_data::AlignmentDataRecordType;
    /// use chain::record::AlignmentDataRecord;
    /// use chainfile as chain;
    ///
    /// let record = AlignmentDataRecord::try_new(
    ///     10,
    ///     Some(0),
    ///     Some(1),
    ///     AlignmentDataRecordType::NonTerminating,
    /// )?;
    ///
    /// assert_eq!(record.size(), 10);
    /// assert_eq!(record.dt(), &Some(0));
    /// assert_eq!(record.dq(), &Some(1));
    /// assert_eq!(
    ///     record.record_type(),
    ///     &AlignmentDataRecordType::NonTerminating
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        size: usize,
        dt: Option<usize>,
        dq: Option<usize>,
        record_type: AlignmentDataRecordType,
    ) -> Result<Self, ParseError> {
        match record_type {
            AlignmentDataRecordType::NonTerminating => {
                if dt.is_none() {
                    return Err(ParseError::InvalidNonTerminatingDt);
                }

                if dq.is_none() {
                    return Err(ParseError::InvalidNonTerminatingDq);
                }
            }
            AlignmentDataRecordType::Terminating => {
                if dt.is_some() {
                    return Err(ParseError::InvalidTerminatingDt);
                }

                if dq.is_some() {
                    return Err(ParseError::InvalidTerminatingDq);
                }
            }
        }

        Ok(AlignmentDataRecord(size, dt, dq, record_type))
    }

    /// Retuns the size of the ungapped alignment.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::record::AlignmentDataRecord;
    /// use chainfile as chain;
    ///
    /// let alignment: AlignmentDataRecord = "9\t1\t0".parse()?;
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
    /// use chain::record::AlignmentDataRecord;
    /// use chainfile as chain;
    ///
    /// let alignment: AlignmentDataRecord = "9\t1\t0".parse()?;
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
    /// use chain::record::AlignmentDataRecord;
    /// use chainfile as chain;
    ///
    /// let alignment: AlignmentDataRecord = "9\t1\t0".parse()?;
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
    /// use chain::record::alignment_data::AlignmentDataRecordType;
    /// use chain::record::AlignmentDataRecord;
    /// use chainfile as chain;
    ///
    /// let alignment: AlignmentDataRecord = "9\t1\t0".parse()?;
    /// assert_eq!(
    ///     *alignment.record_type(),
    ///     AlignmentDataRecordType::NonTerminating
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn record_type(&self) -> &AlignmentDataRecordType {
        &self.3
    }
}

impl FromStr for AlignmentDataRecord {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split(ALIGNMENT_DATA_DELIMITER).collect::<Vec<_>>();
        let record_type = match parts.len() {
            NUM_ALIGNMENT_DATA_FIELDS_TERMINATING => AlignmentDataRecordType::Terminating,
            NUM_ALIGNMENT_DATA_FIELDS_NONTERMINATING => AlignmentDataRecordType::NonTerminating,
            _ => return Err(ParseError::IncorrectNumberOfFields(parts.len())),
        };

        let size = parts[0].parse().map_err(ParseError::InvalidSize)?;
        let (dt, dq) = match record_type {
            AlignmentDataRecordType::NonTerminating => {
                let dt = parts[1].parse().map_err(ParseError::InvalidDt)?;
                let dq = parts[2].parse().map_err(ParseError::InvalidDq)?;
                (Some(dt), Some(dq))
            }
            AlignmentDataRecordType::Terminating => (None, None),
        };

        AlignmentDataRecord::try_new(size, dt, dq, record_type)
    }
}

impl std::fmt::Display for AlignmentDataRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.size())?;

        match self.record_type() {
            AlignmentDataRecordType::NonTerminating => {
                write!(
                    f,
                    "{}{}{}{}",
                    ALIGNMENT_DATA_DELIMITER,
                    self.dt().expect("non-terminating record must have a dt"),
                    ALIGNMENT_DATA_DELIMITER,
                    self.dq().expect("non-terminating record must have a dq"),
                )?;
            }
            AlignmentDataRecordType::Terminating => {}
        }

        Ok(())
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_nonterminating_alignment_data() -> Result<(), Box<dyn std::error::Error>> {
        let record = "9\t1\t0".parse::<AlignmentDataRecord>()?;

        assert_eq!(record.size(), 9);
        assert_eq!(*record.dt(), Some(1));
        assert_eq!(*record.dq(), Some(0));
        assert_eq!(
            *record.record_type(),
            AlignmentDataRecordType::NonTerminating
        );

        Ok(())
    }

    #[test]
    fn test_terminating_alignment_data() -> Result<(), Box<dyn std::error::Error>> {
        let record = "9".parse::<AlignmentDataRecord>()?;

        assert_eq!(record.size(), 9);
        assert_eq!(*record.dt(), None);
        assert_eq!(*record.dq(), None);
        assert_eq!(*record.record_type(), AlignmentDataRecordType::Terminating);

        Ok(())
    }

    #[test]
    fn test_invalid_number_of_fields() -> Result<(), Box<dyn std::error::Error>> {
        let err = "9\t0".parse::<AlignmentDataRecord>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid number of fields in alignment data: expected 3 (non-terminating) or 1 \
             (terminating) fields, found 2 fields"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_size() -> Result<(), Box<dyn std::error::Error>> {
        let err = "?\t0\t1".parse::<AlignmentDataRecord>().unwrap_err();

        assert_eq!(
            err.to_string(),
            "invalid size: invalid digit found in string"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_dt() -> Result<(), Box<dyn std::error::Error>> {
        let err = "9\t?\t1".parse::<AlignmentDataRecord>().unwrap_err();

        assert_eq!(err.to_string(), "invalid dt: invalid digit found in string");

        Ok(())
    }

    #[test]
    fn test_invalid_dq() -> Result<(), Box<dyn std::error::Error>> {
        let err = "9\t0\t?".parse::<AlignmentDataRecord>().unwrap_err();

        assert_eq!(err.to_string(), "invalid dq: invalid digit found in string");

        Ok(())
    }

    #[test]
    fn test_invalid_nonterminating_dt() -> Result<(), Box<dyn std::error::Error>> {
        let err =
            AlignmentDataRecord::try_new(9, None, Some(1), AlignmentDataRecordType::NonTerminating)
                .unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected value for dt in non-terminating alignment data line, found no value"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_nonterminating_dq() -> Result<(), Box<dyn std::error::Error>> {
        let err =
            AlignmentDataRecord::try_new(9, Some(0), None, AlignmentDataRecordType::NonTerminating)
                .unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected value for dq in non-terminating alignment data line, found no value"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_terminating_dt() -> Result<(), Box<dyn std::error::Error>> {
        let err =
            AlignmentDataRecord::try_new(9, Some(0), None, AlignmentDataRecordType::Terminating)
                .unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected no value for dt in terminating alignment data line, found value"
        );

        Ok(())
    }

    #[test]
    fn test_invalid_terminating_dq() -> Result<(), Box<dyn std::error::Error>> {
        let err =
            AlignmentDataRecord::try_new(9, None, Some(1), AlignmentDataRecordType::Terminating)
                .unwrap_err();

        assert_eq!(
            err.to_string(),
            "expected no value for dq in terminating alignment data line, found value"
        );

        Ok(())
    }

    #[test]
    fn test_nonterminating_alignment_data_display() -> Result<(), Box<dyn std::error::Error>> {
        let alignment = AlignmentDataRecord::try_new(
            9,
            Some(1),
            Some(0),
            AlignmentDataRecordType::NonTerminating,
        )?;

        assert_eq!(alignment.to_string(), "9\t1\t0");

        Ok(())
    }

    #[test]
    fn test_terminating_alignment_data_display() -> Result<(), Box<dyn std::error::Error>> {
        let alignment =
            AlignmentDataRecord::try_new(9, None, None, AlignmentDataRecordType::Terminating)?;

        assert_eq!(alignment.to_string(), "9");

        Ok(())
    }
}
