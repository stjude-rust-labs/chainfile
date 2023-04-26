//! A 0-based location within a genome consisting of a contig, a position, and a strand.

use crate::core::Interval;
use crate::core::Strand;

pub mod position;

/// A contiguous molecule upon which a coordinate is located.
pub type Contig = String;
pub use position::Position;

/// An error related to a [`Coordinate`].
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// Attempted to create a Negative-bound coordinate that was not on the
    /// negative strand.
    NegativeBoundOnNonNegativeStrand,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::NegativeBoundOnNonNegativeStrand => write!(
                f,
                "cannot place negative bound coordinate on non-negative strand"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// A zero-based single location, situated upon a contig at a particular
/// position on a given strand.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Coordinate(Contig, Position, Strand);

impl Coordinate {
    /// Creates a new coordinate from the relevant parts.
    ///
    /// Note that a negative bounded coordinate can only sit on the
    /// [`Strand::Negative`], so trying to create a negative bound on the
    /// postive strand will result in an
    /// [`Error::NegativeBoundOnNonNegativeStrand`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    /// use chain::core::coordinate::Error;
    ///
    /// // Positive-stranded
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// assert_eq!(coordinate.contig(), &String::from("seq0"));
    /// assert_eq!(coordinate.position(), &Position::new(0));
    /// assert_eq!(coordinate.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    /// assert_eq!(coordinate.contig(), &String::from("seq0"));
    /// assert_eq!(coordinate.position(), &Position::new(0));
    /// assert_eq!(coordinate.strand(), &Strand::Negative);
    ///
    /// // Negative-bound
    ///
    /// let coordinate = Coordinate::try_new("seq0", Position::negative_bound(), Strand::Negative)?;
    /// assert_eq!(coordinate.contig(), &String::from("seq0"));
    /// assert_eq!(coordinate.position(), &Position::negative_bound());
    /// assert_eq!(coordinate.strand(), &Strand::Negative);
    ///
    /// // Attempting to create negative-bound on positive strand
    ///
    /// let err = Coordinate::try_new("seq0", Position::negative_bound(), Strand::Positive).unwrap_err();
    /// assert_eq!(err, Error::NegativeBoundOnNonNegativeStrand);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(
        contig: impl Into<Contig>,
        position: impl Into<Position>,
        strand: impl Into<Strand>,
    ) -> Result<Coordinate, Error> {
        let contig = contig.into();
        let position = position.into();
        let strand = strand.into();

        if position == Position::NegativeBound && strand != Strand::Negative {
            return Err(Error::NegativeBoundOnNonNegativeStrand);
        }

        Ok(Coordinate(contig, position, strand))
    }

    /// Creates a new [`Coordinate`] with [`Position::NegativeBound`] for the
    /// position from the contig provided.
    ///
    /// Note that a negative bounded coordinate can only sit on the
    /// [`Strand::Negative`], so trying to create a negative bound on the
    /// postive strand will result in an
    /// [`Error::NegativeBoundOnNonNegativeStrand`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    /// assert_eq!(coordinate.position(), &Position::negative_bound());
    /// assert_eq!(coordinate.strand(), &Strand::Negative);
    /// ```
    pub fn negative_bound(contig: impl Into<Contig>) -> Coordinate {
        // SAFETY: this will never fail because we have hardcoded it to live on
        // the negative strand here.
        Coordinate::try_new(contig, Position::NegativeBound, Strand::Negative).unwrap()
    }

    /// Gets the associated [`Contig`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Strand;
    ///
    /// // Zero-based coordinate
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    ///
    /// // Negative-bound coordinate
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// assert_eq!(coordinate.contig().as_str(), "seq0");
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contig(&self) -> &Contig {
        &self.0
    }

    /// Gets the associated [`Position`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// // Zero-based coordinate
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// assert_eq!(coordinate.position(), &Position::new(0));
    ///
    /// // Negative-bound coordinate
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// assert_eq!(coordinate.position(), &Position::negative_bound());
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn position(&self) -> &Position {
        &self.1
    }

    /// Gets the associated [`Strand`] by reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Strand;
    ///
    /// // Zero-based coordinate
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// assert_eq!(coordinate.strand(), &Strand::Positive);
    ///
    /// // Negative-bound coordinate
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// assert_eq!(coordinate.strand(), &Strand::Negative);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn strand(&self) -> &Strand {
        &self.2
    }

    /// Consumes self to attempt to move the coordinate forward by a specified
    /// magnitude.
    ///
    /// A checked add (for positive strand) or subtract (for negative strand) is
    /// performed to ensure we don't overflow/underflow. Note that, though the
    /// position is checked for usize overflow/underflow, we don't do any bounds
    /// checking to make sure that the coordinates fall within any given
    /// interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let result = coordinate.move_forward(10).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(10));
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded
    ///
    /// let coordinate = Coordinate::try_new("seq0", 1000, Strand::Negative)?;
    /// let result = coordinate.move_forward(10).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(990));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded, but with magnitude zero
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let result = coordinate.move_forward(0).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(0));
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded, but with magnitude zero
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    /// let result = coordinate.move_forward(0).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(0));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Negative-stranded to negative-bound
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    /// let result = coordinate.move_forward(1).unwrap()?;
    ///
    /// assert_eq!(result, Coordinate::negative_bound("seq0"));
    ///
    /// // Negative-stranded underflow
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    /// let result = coordinate.move_forward(10);
    ///
    /// assert_eq!(result, None);
    ///
    /// // Negative-bound underflow
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// let result = coordinate.move_forward(10);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_forward(self, magnitude: usize) -> Option<Result<Coordinate, Error>> {
        if magnitude == 0 {
            return Some(Ok(self));
        }

        let position = match (self.position(), self.strand()) {
            (Position::ZeroBased(position), Strand::Positive) => {
                usize::checked_add(*position, magnitude)
            }
            (Position::ZeroBased(position), Strand::Negative) => {
                match usize::checked_sub(*position, magnitude - 1) {
                    Some(0) => return Some(Ok(Coordinate::negative_bound(self.contig()))),
                    Some(result_plus_one) => usize::checked_sub(result_plus_one, 1),
                    None => None,
                }
            }
            // NegativeBound should always be on the negative strand.
            // If, by some chance, it isn't, then the user needs to be
            // alerted of this.
            (Position::NegativeBound, Strand::Positive) => {
                unreachable!("negative bound cannot be on positive strand")
            }
            // If the negative bound _is_ on the negative strand, then moving
            // forward by any magnitude would move the value into the
            // non-allowed negative space.
            (Position::NegativeBound, Strand::Negative) => match magnitude {
                0 => return Some(Ok(self)),
                _ => None,
            },
        }?;

        Some(Coordinate::try_new(
            self.contig(),
            position,
            self.strand().clone(),
        ))
    }

    /// Consumes self to attempt to move the coordinate backward by a specified
    /// magnitude.
    ///
    /// A checked subtract (for positive strand) or add (for negative strand) is
    /// performed to ensure we don't underflow/overflow. Note that, though the
    /// position is checked for usize underflow/overflow, we don't do any bounds
    /// checking to make sure that the coordinates fall within any given
    /// interval.
    ///
    /// Note that we don't handle the negative bound in this method. This is
    /// because it would only be possible to reach the negative bound by moving
    /// backwards on the _positive_ strand, but the negative bound is not
    /// allowed on the positive strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded
    ///
    /// let coordinate = Coordinate::try_new("seq0", 500, Strand::Positive)?;
    /// let result = coordinate.move_backward(10).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(490));
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded
    ///
    /// let coordinate = Coordinate::try_new("seq0", 500, Strand::Negative)?;
    /// let result = coordinate.move_backward(10).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(510));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded, but with magnitude zero
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let result = coordinate.move_backward(0).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(0));
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded, but with magnitude zero
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Negative)?;
    /// let result = coordinate.move_backward(0).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(0));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded underflow
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let result = coordinate.move_backward(10);
    ///
    /// assert_eq!(result, None);
    ///
    /// // Negative-bound
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// let result = coordinate.move_backward(10).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(9));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_backward(self, magnitude: usize) -> Option<Result<Coordinate, Error>> {
        if magnitude == 0 {
            return Some(Ok(self));
        }

        let position = match (self.position(), self.strand()) {
            (Position::ZeroBased(position), Strand::Positive) => {
                usize::checked_sub(*position, magnitude)
            }
            (Position::ZeroBased(position), Strand::Negative) => {
                usize::checked_add(*position, magnitude)
            }
            (Position::NegativeBound, Strand::Positive) => {
                unreachable!()
            }
            (Position::NegativeBound, Strand::Negative) => match magnitude {
                0 => return Some(Ok(self)),
                magnitude => Some(magnitude - 1),
            },
        }?;

        Some(Coordinate::try_new(
            self.contig(),
            position,
            self.strand().clone(),
        ))
    }

    /// Consumes self to attempt to move the coordinate forward by a specified
    /// magnitude while also perform a bounds check within an interval.
    ///
    /// First, a checked add (for positive strand) or subtract (for negative
    /// strand) is performed to ensure we don't overflow/underflow. Next, the
    /// calculated result is checked to ensure it falls within the specified
    /// interval—this is to ensure that, although the `usize` limits are not
    /// broken, the limits of the interval that contains this coordinate are
    /// also not broken.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(10, &interval).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(10));
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::try_new("seq0", 1000, Strand::Negative)?;
    /// let interval = "seq0:1000-0".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(10, &interval).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(990));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded that _does not_ fall within interval
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let interval = "seq0:0-10".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(10, &interval);
    ///
    /// assert_eq!(result, None);
    ///
    /// // Negative-bound that _does not_ fall within interval
    /// // (and also would not move forward due to underflow).
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// let interval = "seq0:10-$".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(12, &interval);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_forward_checked_bounds(
        self,
        magnitude: usize,
        interval: &Interval,
    ) -> Option<Result<Coordinate, Error>> {
        match self.move_forward(magnitude)? {
            Ok(result) => match interval.contains(&result) {
                true => Some(Ok(result)),
                false => None,
            },
            Err(err) => Some(Err(err)),
        }
    }

    /// Consumes self to attempt to move the coordinate backward by a specified
    /// magnitude while also performing a bounds check within an interval.
    ///
    /// First, a checked subtract (for positive strand) or add (for negative
    /// strand) is performed to ensure we don't underflow/overflow. Next, the
    /// calculated result is checked to ensure it falls within the specified
    /// interval—this is to ensure that, although the `usize` limits are not
    /// broken, the limits of the interval that contains this coordinate are
    /// also not broken.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::try_new("seq0", 500, Strand::Positive)?;
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(10, &interval).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(490));
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::try_new("seq0", 500, Strand::Negative)?;
    /// let interval = "seq0:1000-0".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(10, &interval).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(510));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded that _does not_ fall within interval
    ///
    /// let coordinate = Coordinate::try_new("seq0", 0, Strand::Positive)?;
    /// let interval = "seq0:0-10".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(10, &interval);
    ///
    /// assert_eq!(result, None);
    ///
    /// // Negative-bound that _does_ fall within interval
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// let interval = "seq0:10-$".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(11, &interval).unwrap()?;
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), &Position::new(10));
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Negative-bound that _does not_ fall within interval
    ///
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// let interval = "seq0:10-$".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(12, &interval);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_backward_checked_bounds(
        self,
        magnitude: usize,
        interval: &Interval,
    ) -> Option<Result<Coordinate, Error>> {
        match self.move_backward(magnitude)? {
            Ok(result) => match interval.contains(&result) {
                true => Some(Ok(result)),
                false => None,
            },
            Err(err) => Some(Err(err)),
        }
    }

    /// Complements the position with respect to some chromosome size.
    ///
    /// This method _does_ handle the case when the negative bound should be
    /// returned. Note that the strand is _not_ complemented, only the position
    /// is. This is useful when working with chain files, as the chain file will
    /// store the position of the reverse compliment of a sequence without
    /// complementing the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Position;
    /// use chain::core::Strand;
    ///
    /// // Positive-strand
    /// let coordinate = Coordinate::try_new("seq0", 5, Strand::Positive)?;
    /// let inverted = coordinate.complement_position(16)?;
    /// assert_eq!(inverted, Coordinate::try_new("seq0", 10, Strand::Positive)?);
    ///
    /// // Negative-strand
    /// let coordinate = Coordinate::try_new("seq0", 5, Strand::Negative)?;
    /// let inverted = coordinate.complement_position(16)?;
    /// assert_eq!(inverted, Coordinate::try_new("seq0", 10, Strand::Negative)?);
    ///
    /// // Zero-based to negative-bound
    /// let coordinate = Coordinate::try_new("seq0", 16, Strand::Negative)?;
    /// let inverted = coordinate.complement_position(16)?;
    /// assert_eq!(inverted, Coordinate::negative_bound("seq0"));
    ///
    /// // Negative-bound to zero-based
    /// let coordinate = Coordinate::negative_bound("seq0");
    /// let inverted = coordinate.complement_position(16)?;
    /// assert_eq!(inverted, Coordinate::try_new("seq0", 16, Strand::Negative)?);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn complement_position(self, contig_size: usize) -> Result<Coordinate, Error> {
        let position = match self.position() {
            Position::ZeroBased(position) => position,
            Position::NegativeBound => {
                return Coordinate::try_new(self.contig(), contig_size, self.strand().clone())
            }
        };

        let new_position_plus_one = contig_size - position;
        if new_position_plus_one == 0 {
            match self.strand() {
                Strand::Positive => unreachable!(
                    "it should never be possible to want to return the negative bound \
                    if the source coordinate is on the positive strand"
                ),
                Strand::Negative => Ok(Coordinate::negative_bound(self.contig())),
            }
        } else {
            Coordinate::try_new(
                self.contig(),
                new_position_plus_one - 1,
                self.strand().clone(),
            )
        }
    }
}

impl std::fmt::Display for Coordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{{}, {}, {}}}",
            self.contig(),
            self.position(),
            self.strand()
        )
    }
}
