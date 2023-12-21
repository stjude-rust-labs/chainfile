//! A 0-based, half-open interval consisting of a start and stop coordinate.
//!
//! Intervals can be either on the positive strand or the negative strand of a
//! contig.
//!
//! ```text
//! ================ seq0 ===============
//!
//! | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
//! -------------------------------------
//! |   |   | X | X | X | X | O |   |   |  <= seq0:3-7
//! |   |   |   |   | O | X | X | X | X |  <= seq0:9-5
//! ```
//!
//! - The first interval above (`seq0:3-7`) is a forward-stranded interval from
//!   3 up until (but not including) 7.
//! - The second interval above (`seq0:9-5`) is a reverse-stranded interval from
//!   9 down until (but not including) 5.
//!
//! ## Parsing Intervals
//!
//! Two types of intervals can be parsed from a string: singular positions or
//! ranges of positions.
//!
//! - Singular positions are provided in the form `<contig>:<position>` (e.g.,
//!   `seq0:1`). Singular positions are always parsed as an interval containing
//!   only the position specified. These are always placed on the positive
//!   strand, primarily because there is no coherent way to derive the
//!   strandedness from the interval short of explicitly specifying the strand
//!   in the interval. The authors did not want to support this, so we chose to
//!   always interpret these as being on the positive strand. Furthermore, when
//!   an interval is a single position, directionality is largely
//!   inconsequential for the purposes of this library.
//! - Range of positions are provided in the form `<contig>:<start>-<end>`
//!   (e.g., `seq0:0-1000`). Importantly, the start and the end positions are
//!   used to determine which strand the interval falls on. If the end position
//!   is greater than the start position, then the interval is interpretted as
//!   being [`Strand::Positive`]. If the start position is greater than the end
//!   position, then the interval is interpretted as being [`Strand::Negative`].
//!   Note that, if the start position and end position are equal, this would be
//!   a zero-sized interval, which would trigger an
//!   [`Error::ZeroSizedInterval`].

use std::collections::VecDeque;
use std::str::FromStr;

use crate::core::coordinate;
use crate::core::coordinate::Position;
use crate::core::Contig;
use crate::core::Coordinate;
use crate::core::Strand;

/// An error related to an interval.
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// Attempted to create an invalid coordinate.
    InvalidCoordinate(coordinate::Error),
    /// Could not create an interval for coordinates on different contigs.
    NonsensicalContigs,
    /// Could not create an interval for coordinates on different strands.
    NonsensicalStrands,
    /// The start position cannot equal the end position, which would result in
    /// a zero-size interval.
    ZeroSizedInterval,
    /// The start position is greater than the end position for a positive
    /// stranded interval, which is not allowed.
    StartGreaterThanEndForPositiveStrand,
    /// The end position is greater than the start position for a negative
    /// stranded interval, which is not allowed.
    EndGreaterThanStartForNegativeStrand,
    /// Attempted to clamp an interval, but the two intervals had mismatching
    /// contigs.
    MismatchedContigDuringClamp(String, String),
    /// Attempted to clamp an interval, but the two intervals had mismatching
    /// strands.
    MismatchedStrandDuringClamp(Strand, Strand),
    /// Could not parse an interval from the given value.
    ParseError(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::InvalidCoordinate(err) => write!(f, "invalid coordinate: {}", err),
            Error::NonsensicalContigs => {
                write!(f, "interval coordinates must be on the same contig")
            }
            Error::NonsensicalStrands => {
                write!(f, "interval coordinates must be on the same strand")
            }
            Error::ZeroSizedInterval => write!(
                f,
                "start position equals end position, which is a zero-sized interval"
            ),
            Error::StartGreaterThanEndForPositiveStrand => {
                write!(
                    f,
                    "start position cannot be greater than the end position for a positive \
                     stranded interval"
                )
            }
            Error::EndGreaterThanStartForNegativeStrand => {
                write!(
                    f,
                    "end position cannot be greater than the start position for a negative \
                     stranded interval"
                )
            }
            Error::MismatchedContigDuringClamp(a, b) => {
                write!(f, "mismatched contig while clamping: {} and {}", a, b)
            }
            Error::MismatchedStrandDuringClamp(a, b) => {
                write!(f, "mismatched strand while clamping: {} and {}", a, b)
            }
            Error::ParseError(val) => write!(f, "could not parse interval from the value: {}", val),
        }
    }
}

impl std::error::Error for Error {}

/// A half-open interval consisting of a start and end coordinate.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Interval(Coordinate, Coordinate);

impl Interval {
    /// Attempts to create a new [`Interval`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// // Positive-stranded interval
    ///
    /// let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// // Negative-stranded interval
    ///
    /// let start = Coordinate::try_new("seq0", 1000, Strand::Negative)?;
    /// let end = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    ///
    /// Interval::try_new(start, end)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(start: Coordinate, end: Coordinate) -> Result<Interval, Error> {
        // (1) Check that the coordinates are from the same contig.
        if start.contig() != end.contig() {
            return Err(Error::NonsensicalContigs);
        }

        // (2) Check that the coordinates are from the same strand.
        if start.strand() != end.strand() {
            return Err(Error::NonsensicalStrands);
        }

        // (3) Ensure that the start position does not equal the end position,
        // which would result in a zero-sized interval.
        if start.position() == end.position() {
            return Err(Error::ZeroSizedInterval);
        }

        // (4) Ensure that the magnitudes of the positions are oriented
        // correctly based on the strand value.
        match start.strand() {
            Strand::Positive => {
                if start.position() > end.position() {
                    return Err(Error::StartGreaterThanEndForPositiveStrand);
                }
            }
            Strand::Negative => {
                if end.position() > start.position() {
                    return Err(Error::EndGreaterThanStartForNegativeStrand);
                }
            }
        }

        Ok(Interval(start, end))
    }

    /// Gets the start position of the interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;
    ///
    /// let interval = Interval::try_new(start.clone(), end)?;
    /// assert_eq!(interval.start(), &start);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn start(&self) -> &Coordinate {
        &self.0
    }

    /// Gets the end position of the interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;
    ///
    /// let interval = Interval::try_new(start, end.clone())?;
    /// assert_eq!(interval.end(), &end);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn end(&self) -> &Coordinate {
        &self.1
    }

    /// Gets the contig for this interval. Note that the contig for the start
    /// and end coordinate are guaranteed to be the same by the checks done at
    /// [`Interval`] creation time.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;
    ///
    /// let interval = Interval::try_new(start, end)?;
    /// assert_eq!(interval.contig(), &String::from("seq0"));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        self.0.contig()
    }

    /// Gets the strand for this interval. Note that the strand for the start
    /// and end coordinate are guaranteed to be the same by the checks done at
    /// [`Interval`] creation time.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;
    ///
    /// let interval = Interval::try_new(start, end)?;
    /// assert_eq!(interval.strand(), &Strand::Positive);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> &Strand {
        self.0.strand()
    }

    /// Gets the distance of the interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// // Positive-stranded interval
    ///
    /// let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;
    ///
    /// let interval = Interval::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// // Negative-stranded interval
    ///
    /// let start = Coordinate::try_new("seq0", 1000, Strand::Negative)?;
    /// let end = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    ///
    /// let interval = Interval::try_new(start, end)?;
    /// assert_eq!(interval.distance(), 1000);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn distance(&self) -> usize {
        match self.strand() {
            Strand::Positive => self.end().position() - self.start().position(),
            Strand::Negative => self.start().position() - self.end().position(),
        }
    }

    /// Indicates whether a coordinate falls within an interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    ///
    /// assert!(interval.contains(&Coordinate::try_new("seq0", 0, Strand::Positive)?));
    /// assert!(interval.contains(&Coordinate::try_new("seq0", 999, Strand::Positive)?));
    ///
    /// assert!(!interval.contains(&Coordinate::try_new("seq1", 0, Strand::Positive)?));
    /// assert!(!interval.contains(&Coordinate::try_new("seq0", 50, Strand::Negative)?));
    /// assert!(!interval.contains(&Coordinate::try_new("seq0", 1000, Strand::Positive)?));
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:1000-0".parse::<Interval>()?;
    ///
    /// assert!(interval.contains(&Coordinate::try_new("seq0", 1000, Strand::Negative)?));
    /// assert!(interval.contains(&Coordinate::try_new("seq0", 1, Strand::Negative)?));
    ///
    /// assert!(!interval.contains(&Coordinate::try_new("seq1", 1000, Strand::Negative)?));
    /// assert!(!interval.contains(&Coordinate::try_new("seq0", 1000, Strand::Positive)?));
    /// assert!(!interval.contains(&Coordinate::try_new("seq0", 0, Strand::Negative)?));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contains(&self, coordinate: &Coordinate) -> bool {
        if self.contig() != coordinate.contig() {
            return false;
        }

        if self.strand() != coordinate.strand() {
            return false;
        }

        match self.strand() {
            Strand::Positive => {
                if self.start().position() > coordinate.position() {
                    return false;
                }

                if self.end().position() <= coordinate.position() {
                    return false;
                }
            }
            Strand::Negative => {
                if self.start().position() < coordinate.position() {
                    return false;
                }

                if self.end().position() >= coordinate.position() {
                    return false;
                }
            }
        }

        true
    }

    /// Consumes self to clamp an [`Interval`] based on the value of another
    /// interval.
    ///
    /// This method can, at most, shorten the consumed interval to be within the
    /// `interval` argument's bounds. If either end of the `interval` argument
    /// falls within the consumed interval's bounds, then the returned interval
    /// will be shortened to match the `interval` argument's start/end position.
    /// If the `interval` argument completely encapsulates the consumed
    /// interval, then the original interval will be returned unmodified.
    ///
    /// Note that this function returns `Result<Interval, Error>` instead of an
    /// `Option` because its expected that, if you are clamping an interval,
    /// you're probably doing so with the expectation that the intervals are on
    /// the same contig and strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::interval::Error;
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let clamp = "seq0:250-2000".parse::<Interval>()?;
    /// let result = interval.clamp(&clamp)?;
    /// assert_eq!(&result, &"seq0:250-1000".parse::<Interval>()?);
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:2000-1000".parse::<Interval>()?;
    /// let clamp = "seq0:3000-1250".parse::<Interval>()?;
    /// let result = interval.clamp(&clamp)?;
    /// assert_eq!(&result, &"seq0:2000-1250".parse::<Interval>()?);
    ///
    /// // Differing contigs
    ///
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let clamp = "seq1:250-2000".parse::<Interval>()?;
    /// let result = interval.clamp(&clamp);
    /// assert!(matches!(
    ///     result,
    ///     Err(Error::MismatchedContigDuringClamp(_, _))
    /// ));
    ///
    /// // Differing strands
    ///
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let clamp = "seq0:2000-250".parse::<Interval>()?;
    /// let result = interval.clamp(&clamp);
    /// assert!(matches!(
    ///     result,
    ///     Err(Error::MismatchedStrandDuringClamp(_, _))
    /// ));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn clamp(self, interval: &Interval) -> Result<Interval, Error> {
        if self.contig() != interval.contig() {
            return Err(Error::MismatchedContigDuringClamp(
                self.contig().clone(),
                interval.contig().clone(),
            ));
        }

        if self.strand() != interval.strand() {
            return Err(Error::MismatchedStrandDuringClamp(
                self.strand().clone(),
                interval.strand().clone(),
            ));
        }

        let start_position = match self.strand() {
            Strand::Positive => std::cmp::max(self.start().position(), interval.start().position()),
            Strand::Negative => std::cmp::min(self.start().position(), interval.start().position()),
        };

        let end_position = match self.strand() {
            Strand::Positive => std::cmp::min(self.end().position(), interval.end().position()),
            Strand::Negative => std::cmp::max(self.end().position(), interval.end().position()),
        };

        let start =
            Coordinate::try_new(self.contig(), start_position.clone(), self.strand().clone())
                .map_err(Error::InvalidCoordinate)?;
        let end = Coordinate::try_new(self.contig(), end_position.clone(), self.strand().clone())
            .map_err(Error::InvalidCoordinate)?;

        Interval::try_new(start, end)
    }

    /// Gets the offset of the provided coordinate from the starting position of
    /// the interval.
    ///
    /// - If the provided coordinate _does not_ fall within the interval, then
    ///   `None` is returned.
    /// - If the provided coordinate _does_ fall within the interval, then the
    ///   distance from the starting position of the interval is returned
    ///   wrapped in `Some`.
    ///
    /// Note that for [`Strand::Positive`] intervals, the magnitude should be
    /// interpreted as a _positive_ magnitude, whereas a [`Strand::Negative`]
    /// offset will be a _negative_ magnitude (though it will not be reflected
    /// in the size of the returned result).
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    ///
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 5, Strand::Positive)?),
    ///     Some(5)
    /// );
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 999, Strand::Positive)?),
    ///     Some(999)
    /// );
    ///
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq1", 5, Strand::Positive)?),
    ///     None
    /// );
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 5, Strand::Negative)?),
    ///     None
    /// );
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 1000, Strand::Positive)?),
    ///     None
    /// );
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:1000-0".parse::<Interval>()?;
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 995, Strand::Negative)?),
    ///     Some(5)
    /// );
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 1, Strand::Negative)?),
    ///     Some(999)
    /// );
    ///
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq1", 995, Strand::Negative)?),
    ///     None
    /// );
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 995, Strand::Positive)?),
    ///     None
    /// );
    /// assert_eq!(
    ///     interval.offset_from_start(&Coordinate::try_new("seq0", 0, Strand::Negative)?),
    ///     None
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn offset_from_start(&self, coordinate: &Coordinate) -> Option<usize> {
        if !self.contains(coordinate) {
            return None;
        }

        match self.strand() {
            Strand::Positive => Some(coordinate.position() - self.start().position()),
            Strand::Negative => Some(self.start().position() - coordinate.position()),
        }
    }

    /// Attempts to translate an offset from the start of the interval to a
    /// coordinate.
    ///
    /// - If the translated coordinate _does not_ fall within the interval, then
    ///   `None` is returned.
    /// - If the translated coordinate _does_ fall within the interval, then the
    ///   coordinate is returned wrapped in `Some`.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chainfile as chain;
    ///
    /// // Positive-stranded interval
    ///
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    ///
    /// assert_eq!(
    ///     interval.translate_offset_from_start(5),
    ///     Some(Ok(Coordinate::try_new("seq0", 5, Strand::Positive)?))
    /// );
    /// assert_eq!(
    ///     interval.translate_offset_from_start(999),
    ///     Some(Ok(Coordinate::try_new("seq0", 999, Strand::Positive)?))
    /// );
    /// assert_eq!(interval.translate_offset_from_start(1000), None);
    ///
    /// // Negative-stranded interval
    ///
    /// let interval = "seq0:1000-$".parse::<Interval>()?;
    ///
    /// assert_eq!(
    ///     interval.translate_offset_from_start(5),
    ///     Some(Ok(Coordinate::try_new("seq0", 995, Strand::Negative)?))
    /// );
    /// assert_eq!(
    ///     interval.translate_offset_from_start(1000),
    ///     Some(Ok(Coordinate::try_new("seq0", 0, Strand::Negative)?))
    /// );
    /// assert_eq!(interval.translate_offset_from_start(1001), None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn translate_offset_from_start(&self, offset: usize) -> Option<Result<Coordinate, Error>> {
        let coordinate = self
            .start()
            .clone()
            .move_forward(offset)?
            .map_err(Error::InvalidCoordinate);

        match coordinate {
            Ok(coordinate) => match self.contains(&coordinate) {
                true => Some(Ok(coordinate)),
                false => None,
            },
            Err(err) => Some(Err(err)),
        }
    }
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} ({})",
            self.contig(),
            self.start().position(),
            self.end().position(),
            self.strand()
        )
    }
}

impl FromStr for Interval {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut parts = s.split(':').collect::<VecDeque<_>>();

        // (1) Check to ensure that there is one and only one colon, which
        // separates the contig from the position range.
        if parts.len() != 2 {
            return Err(Error::ParseError(s.to_string()));
        }

        // SAFETY: we checked that there are two parts above.
        let contig = parts.pop_front().unwrap();
        let mut position_parts = parts
            .pop_front()
            .unwrap()
            .split('-')
            .collect::<VecDeque<_>>();

        // (2) Parse the single position or range.
        match position_parts.len() {
            1 => {
                let start_position = parse_position(position_parts.pop_front().unwrap())?;

                let start = Coordinate::try_new(contig, start_position, Strand::Positive)
                    .map_err(Error::InvalidCoordinate)?;
                let end = start
                    .clone()
                    .move_forward(1)
                    .unwrap()
                    .map_err(Error::InvalidCoordinate)?;

                Interval::try_new(start, end)
            }
            2 => {
                // SAFETY: we just checked if there are at least two items.
                let start_position = parse_position(position_parts.pop_front().unwrap())?;
                let end_position = parse_position(position_parts.pop_front().unwrap())?;

                let strand = if end_position >= start_position {
                    Strand::Positive
                } else {
                    Strand::Negative
                };

                Interval::try_new(
                    Coordinate::try_new(contig, start_position, strand.clone())
                        .map_err(Error::InvalidCoordinate)?,
                    Coordinate::try_new(contig, end_position, strand)
                        .map_err(Error::InvalidCoordinate)?,
                )
            }
            _ => Err(Error::ParseError(s.to_string())),
        }
    }
}

fn parse_position(position: &str) -> Result<Position, Error> {
    match position {
        "$" => Ok(Position::negative_bound()),
        other => match other.parse::<usize>() {
            Ok(position) => Ok(Position::new(position)),
            Err(_) => Err(Error::ParseError(position.to_string())),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Interval;
    use crate::core::Strand;

    #[test]
    fn test_it_creates_a_valid_interval() -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
        let end = Coordinate::try_new("seq0", 10, Strand::Positive)?;

        let interval = Interval::try_new(start, end)?;
        assert_eq!(interval.distance(), 10);

        Ok(())
    }

    #[test]
    fn test_it_errors_when_contigs_differ() -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
        let end = Coordinate::try_new("chr2", 10, Strand::Positive)?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::NonsensicalContigs));

        Ok(())
    }

    #[test]
    fn test_it_errors_when_strands_differ() -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
        let end = Coordinate::try_new("seq0", 10, Strand::Negative)?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::NonsensicalStrands));

        Ok(())
    }

    #[test]
    fn test_it_does_not_create_a_zero_sized_interval() -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", 0, Strand::Positive)?;
        let end = Coordinate::try_new("seq0", 0, Strand::Positive)?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::ZeroSizedInterval));

        Ok(())
    }

    #[test]
    fn test_it_does_not_allow_start_to_be_greater_than_end_for_positive_stranded_interval()
    -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", 10, Strand::Positive)?;
        let end = Coordinate::try_new("seq0", 0, Strand::Positive)?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::StartGreaterThanEndForPositiveStrand));

        Ok(())
    }

    #[test]
    fn test_it_does_not_allow_end_to_be_greater_than_start_for_negative_stranded_interval()
    -> Result<(), Box<dyn std::error::Error>> {
        let start = Coordinate::try_new("seq0", 0, Strand::Negative)?;
        let end = Coordinate::try_new("seq0", 10, Strand::Negative)?;

        let err = Interval::try_new(start, end).unwrap_err();
        assert!(matches!(err, Error::EndGreaterThanStartForNegativeStrand));

        Ok(())
    }

    #[test]
    fn test_it_clamps_correctly_for_positive_stranded_intervals()
    -> Result<(), Box<dyn std::error::Error>> {
        let interval = "seq0:1000-2000".parse::<Interval>()?;
        assert!(matches!(
            interval.clone().clamp(&"seq1:0-3000".parse::<Interval>()?),
            Err(Error::MismatchedContigDuringClamp(_, _))
        ));
        assert!(matches!(
            interval.clone().clamp(&"seq0:3000-0".parse::<Interval>()?),
            Err(Error::MismatchedStrandDuringClamp(_, _))
        ));
        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:0-3000".parse::<Interval>()?)?,
            "seq0:1000-2000".parse::<Interval>()?
        );
        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:1250-3000".parse::<Interval>()?)?,
            "seq0:1250-2000".parse::<Interval>()?
        );
        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:0-1750".parse::<Interval>()?)?,
            "seq0:1000-1750".parse::<Interval>()?
        );
        assert_eq!(
            interval.clamp(&"seq0:1250-1750".parse::<Interval>()?)?,
            "seq0:1250-1750".parse::<Interval>()?
        );

        Ok(())
    }

    #[test]
    fn test_it_clamps_correctly_for_negative_stranded_intervals()
    -> Result<(), Box<dyn std::error::Error>> {
        let interval = "seq0:2000-1000".parse::<Interval>()?;
        assert!(matches!(
            interval.clone().clamp(&"seq1:3000-0".parse::<Interval>()?),
            Err(Error::MismatchedContigDuringClamp(_, _))
        ));
        assert!(matches!(
            interval.clone().clamp(&"seq0:0-3000".parse::<Interval>()?),
            Err(Error::MismatchedStrandDuringClamp(_, _))
        ));
        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:3000-0".parse::<Interval>()?)?,
            "seq0:2000-1000".parse::<Interval>()?
        );
        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:3000-1250".parse::<Interval>()?)?,
            "seq0:2000-1250".parse::<Interval>()?
        );
        assert_eq!(
            interval
                .clone()
                .clamp(&"seq0:1750-0".parse::<Interval>()?)?,
            "seq0:1750-1000".parse::<Interval>()?
        );
        assert_eq!(
            interval.clamp(&"seq0:1750-1250".parse::<Interval>()?)?,
            "seq0:1750-1250".parse::<Interval>()?
        );

        Ok(())
    }

    #[test]
    fn test_parsing_intervals_works_for_valid_single_position()
    -> Result<(), Box<dyn std::error::Error>> {
        let interval = "seq0:1".parse::<Interval>()?;
        let start = Coordinate::try_new("seq0", 1, Strand::Positive)?;
        let end = Coordinate::try_new("seq0", 2, Strand::Positive)?;

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);
        Ok(())
    }

    #[test]
    fn test_parsing_intervals_works_for_valid_position_range()
    -> Result<(), Box<dyn std::error::Error>> {
        // Testing positive stranded interval
        let interval = "seq0:1-1000".parse::<Interval>()?;
        let start = Coordinate::try_new("seq0", 1, Strand::Positive)?;
        let end = Coordinate::try_new("seq0", 1000, Strand::Positive)?;

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);

        // Testing negative stranded interval
        let interval = "seq0:1000-1".parse::<Interval>()?;
        let start = Coordinate::try_new("seq0", 1000, Strand::Negative)?;
        let end = Coordinate::try_new("seq0", 1, Strand::Negative)?;

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);

        // Testing running up until negative bound
        let interval = "seq0:1000-$".parse::<Interval>()?;
        let start = Coordinate::try_new("seq0", 1000, Strand::Negative)?;
        let end = Coordinate::negative_bound("seq0");

        assert_eq!(*interval.start(), start);
        assert_eq!(*interval.end(), end);
        Ok(())
    }

    #[test]
    fn test_various_invalid_intervals() -> Result<(), Box<dyn std::error::Error>> {
        let err = "1".parse::<Interval>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "1-1000".parse::<Interval>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "seq0:".parse::<Interval>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "seq0:1-".parse::<Interval>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        let err = "seq0:1-10000:".parse::<Interval>().unwrap_err();
        assert!(matches!(err, Error::ParseError(_)));

        Ok(())
    }

    #[test]
    fn test_interval_to_string() -> Result<(), Box<dyn std::error::Error>> {
        // Positive stranded interval
        let start = Coordinate::try_new("seq0", Position::new(0), Strand::Positive)?;
        let end = Coordinate::try_new("seq0", Position::new(10), Strand::Positive)?;
        let interval = Interval::try_new(start, end)?;

        assert_eq!(interval.to_string(), "seq0:0-10 (+)");

        // Negative stranded interval
        let start = Coordinate::try_new("seq0", Position::new(10), Strand::Negative)?;
        let end = Coordinate::try_new("seq0", Position::new(0), Strand::Negative)?;
        let interval = Interval::try_new(start, end)?;

        assert_eq!(interval.to_string(), "seq0:10-0 (-)");

        // Negative stranded interval with negative bound
        let start = Coordinate::try_new("seq0", Position::new(10), Strand::Negative)?;
        let end = Coordinate::negative_bound("seq0");
        let interval = Interval::try_new(start, end)?;

        assert_eq!(interval.to_string(), "seq0:10-$ (-)");

        Ok(())
    }
}
