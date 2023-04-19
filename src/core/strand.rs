//! The strand upon which a coordinate is located.

use std::io;
use std::str::FromStr;

/// An error related to the parsing of a strand.
#[derive(Debug)]
pub struct ParseStrandError(io::Error);

impl std::fmt::Display for ParseStrandError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "parse strand error: {}", self.0)
    }
}

impl std::error::Error for ParseStrandError {}

/// The strand of a header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Strand {
    /// The positive strand (`+`).
    Positive,
    /// The negative strand (`-`).
    Negative,
}

impl FromStr for Strand {
    type Err = ParseStrandError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Self::Positive),
            "-" => Ok(Self::Negative),
            c => Err(ParseStrandError(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{} is not a valid strand", c),
            ))),
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Positive => write!(f, "+"),
            Strand::Negative => write!(f, "-"),
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_strand_from_str() -> Result<(), Box<dyn std::error::Error>> {
        let strand: Strand = "+".parse()?;
        assert_eq!(strand, Strand::Positive);

        let strand: Strand = "-".parse()?;
        assert_eq!(strand, Strand::Negative);

        let err = "?".parse::<Strand>().unwrap_err();
        assert_eq!(
            err.to_string(),
            "parse strand error: ? is not a valid strand"
        );

        Ok(())
    }

    #[test]
    fn test_strand_display() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(Strand::Positive.to_string(), "+");
        assert_eq!(Strand::Negative.to_string(), "-");
        Ok(())
    }
}
