//! A 0-based location within a genome consisting of a contig, a position, and a strand.

use crate::core::Interval;
use crate::core::Strand;

/// A contiguous molecule upon which a coordinate is located.
pub type Contig = String;
/// The exact, 0-based position upon a contig which a coordinate is located.
pub type Position = usize;

/// A zero-based single location, situated upon a contig at a particular
/// position on a given strand.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Coordinate(Contig, Position, Strand);

impl Coordinate {
    /// Creates a new [`Coordinate`] from the relevant parts.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Strand;
    ///
    /// let coordinate = Coordinate::new("chr1", 0, Strand::Positive);
    /// ```
    pub fn new(contig: impl Into<Contig>, position: Position, strand: Strand) -> Coordinate {
        Coordinate(contig.into(), position, strand)
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
    /// let coordinate = Coordinate::new("chr1", 0, Strand::Positive);
    /// assert_eq!(coordinate.contig().as_str(), "chr1");
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
    /// use chain::core::Strand;
    ///
    /// let coordinate = Coordinate::new("chr1", 0, Strand::Positive);
    /// assert_eq!(coordinate.position(), 0);
    /// ```
    pub fn position(&self) -> Position {
        self.1
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
    /// let coordinate = Coordinate::new("chr1", 0, Strand::Positive);
    /// assert_eq!(coordinate.strand(), &Strand::Positive);
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
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded
    ///
    /// let coordinate = Coordinate::new("seq0", 0, Strand::Positive);
    /// let result = coordinate.move_forward(10).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 10);
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded
    ///
    /// let coordinate = Coordinate::new("seq0", 1000, Strand::Negative);
    /// let result = coordinate.move_forward(10).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 990);
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Negative-stranded underflow
    ///
    /// let coordinate = Coordinate::new("seq0", 0, Strand::Negative);
    /// let result = coordinate.move_forward(10);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_forward(self, magnitude: usize) -> Option<Coordinate> {
        let position = match self.strand() {
            Strand::Positive => usize::checked_add(self.position(), magnitude)?,
            Strand::Negative => usize::checked_sub(self.position(), magnitude)?,
        };

        Some(Coordinate::new(
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
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded
    ///
    /// let coordinate = Coordinate::new("seq0", 500, Strand::Positive);
    /// let result = coordinate.move_backward(10).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 490);
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded
    ///
    /// let coordinate = Coordinate::new("seq0", 500, Strand::Negative);
    /// let result = coordinate.move_backward(10).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 510);
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded underflow
    ///
    /// let coordinate = Coordinate::new("seq0", 0, Strand::Positive);
    /// let result = coordinate.move_backward(10);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_backward(self, magnitude: usize) -> Option<Coordinate> {
        let position = match self.strand() {
            Strand::Positive => usize::checked_sub(self.position(), magnitude)?,
            Strand::Negative => usize::checked_add(self.position(), magnitude)?,
        };

        Some(Coordinate::new(
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
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::new("seq0", 0, Strand::Positive);
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(10, &interval).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 10);
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::new("seq0", 1000, Strand::Negative);
    /// let interval = "seq0:1000-0".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(10, &interval).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 990);
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded that _does not_ fall within interval
    ///
    /// let coordinate = Coordinate::new("seq0", 0, Strand::Positive);
    /// let interval = "seq0:0-10".parse::<Interval>()?;
    /// let result = coordinate.move_forward_checked_bounds(10, &interval);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_forward_checked_bounds(
        self,
        magnitude: usize,
        interval: &Interval,
    ) -> Option<Coordinate> {
        let position = match self.strand() {
            Strand::Positive => usize::checked_add(self.position(), magnitude)?,
            Strand::Negative => usize::checked_sub(self.position(), magnitude)?,
        };

        let result = Coordinate::new(self.contig(), position, self.strand().clone());

        match interval.contains(&result) {
            true => Some(result),
            false => None,
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
    /// use chain::core::Strand;
    ///
    /// // Positive-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::new("seq0", 500, Strand::Positive);
    /// let interval = "seq0:0-1000".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(10, &interval).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 490);
    /// assert_eq!(result.strand(), &Strand::Positive);
    ///
    /// // Negative-stranded that falls within interval
    ///
    /// let coordinate = Coordinate::new("seq0", 500, Strand::Negative);
    /// let interval = "seq0:1000-0".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(10, &interval).unwrap();
    ///
    /// assert_eq!(result.contig(), &String::from("seq0"));
    /// assert_eq!(result.position(), 510);
    /// assert_eq!(result.strand(), &Strand::Negative);
    ///
    /// // Positive-stranded that _does not_ fall within interval
    ///
    /// let coordinate = Coordinate::new("seq0", 0, Strand::Positive);
    /// let interval = "seq0:0-10".parse::<Interval>()?;
    /// let result = coordinate.move_backward_checked_bounds(10, &interval);
    ///
    /// assert_eq!(result, None);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn move_backward_checked_bounds(
        self,
        magnitude: usize,
        interval: &Interval,
    ) -> Option<Coordinate> {
        let position = match self.strand() {
            Strand::Positive => usize::checked_sub(self.position(), magnitude)?,
            Strand::Negative => usize::checked_add(self.position(), magnitude)?,
        };

        let result = Coordinate::new(self.contig(), position, self.strand().clone());

        match interval.contains(&result) {
            true => Some(result),
            false => None,
        }
    }
}
