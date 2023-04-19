//! Pairs of intervals that match contiguously.

use crate::core::interval;
use crate::core::Coordinate;
use crate::core::Interval;

/// An error related to constructing a contiguous interval pair.
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// The two intervals don't have the same size, so they can't linearly map
    /// to one another.
    DistancesDontMatch(usize, usize),
    /// An invalid interval was detected.
    InvalidInterval(interval::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::DistancesDontMatch(reference, query) => write!(
                f,
                "reference interval distance ({}) doesn't match query interval distance ({})",
                reference, query
            ),
            Error::InvalidInterval(err) => write!(f, "invalid interval: {}", err),
        }
    }
}

impl std::error::Error for Error {}

/// A utility struct which contains a linearly mapped segment of both the
/// reference and the query sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContiguousIntervalPair(Interval, Interval);

impl ContiguousIntervalPair {
    /// Attempts to create a new [`ContiguousIntervalPair`] from a reference
    /// interval and a query interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Interval;
    /// use chain::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:0-1000".parse::<Interval>()?;
    /// ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(reference: Interval, query: Interval) -> Result<Self, Error> {
        if reference.distance() != query.distance() {
            return Err(Error::DistancesDontMatch(
                reference.distance(),
                query.distance(),
            ));
        }

        Ok(Self(reference, query))
    }

    /// Gets the reference interval by reference for the interval pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Interval;
    /// use chain::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:0-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference.clone(), query)?;
    ///
    /// assert_eq!(pair.reference(), &reference);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference(&self) -> &Interval {
        &self.0
    }

    /// Gets the query interval by reference for the interval pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Interval;
    /// use chain::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:0-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query.clone())?;
    ///
    /// assert_eq!(pair.query(), &query);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query(&self) -> &Interval {
        &self.1
    }

    /// Lifts over a coordinate within the reference coordinate system to the
    /// query coordinate system.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chain::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    ///
    /// // Positive-stranded to positive-stranded
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:1000-2000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let old = Coordinate::new("seq0", 50, Strand::Positive);
    /// let new = Coordinate::new("seq1", 1050, Strand::Positive);
    /// let lifted = pair.liftover(&old).unwrap();
    ///
    /// assert_eq!(new, lifted);
    ///
    /// // Positive-stranded to negative-stranded
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:2000-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let old = Coordinate::new("seq0", 50, Strand::Positive);
    /// let new = Coordinate::new("seq1", 1950, Strand::Negative);
    /// let lifted = pair.liftover(&old).unwrap();
    ///
    /// assert_eq!(new, lifted);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn liftover(&self, coordinate: &Coordinate) -> Option<Coordinate> {
        let offset = self.reference().offset_from_start(coordinate)?;
        self.query().translate_offset_from_start(offset)
    }

    /// Consumes self to clamp an interval pair to the specified `interval` for
    /// the reference interval. The query interval is similarly lifted over and clamped.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile as chain;
    /// use chain::core::Coordinate;
    /// use chain::core::Interval;
    /// use chain::core::Strand;
    /// use chain::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    ///
    /// // Positive-stranded to positive-stranded
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:1000-2000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let interval = "seq0:50-51".parse::<Interval>()?;
    /// let result = pair.clamp(&interval)?;
    ///
    /// assert_eq!(result.reference().start(), &Coordinate::new("seq0", 50, Strand::Positive));
    /// assert_eq!(result.reference().end(), &Coordinate::new("seq0", 51, Strand::Positive));
    /// assert_eq!(result.query().start(), &Coordinate::new("seq1", 1050, Strand::Positive));
    /// assert_eq!(result.query().end(), &Coordinate::new("seq1", 1051, Strand::Positive));
    ///
    /// // Positive-stranded to negative-stranded
    ///
    /// let reference = "seq0:0-1000".parse::<Interval>()?;
    /// let query = "seq1:2000-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let interval = "seq0:50-51".parse::<Interval>()?;
    /// let result = pair.clamp(&interval)?;
    ///
    /// assert_eq!(result.reference().start(), &Coordinate::new("seq0", 50, Strand::Positive));
    /// assert_eq!(result.reference().end(), &Coordinate::new("seq0", 51, Strand::Positive));
    /// assert_eq!(result.query().start(), &Coordinate::new("seq1", 1950, Strand::Negative));
    /// assert_eq!(result.query().end(), &Coordinate::new("seq1", 1949, Strand::Negative));
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn clamp(self, interval: &Interval) -> Result<ContiguousIntervalPair, Error> {
        let reference = self
            .reference()
            .clone()
            .clamp(interval)
            .map_err(Error::InvalidInterval)?;

        let query_start = self.liftover(reference.start()).unwrap();

        // Note that the position is moved backward with a bounds check because
        // we always expect the _end_ of an interval minus one to fall within
        // the interval again. I don't feel that the bounds check is strictly
        // required since this assumption is trivially known, but we do it here
        // nonetheless.
        //
        // Adding back the one, however, is a completely different matter. Since
        // the ending coordinate is not included in the interval, we expect that
        // the end plus one will fall outside of the interval any time the end
        // coordinate pointed to the end of a [`ContiguousIntervalPair`]. Thus,
        // we _must_ do the unchecked `move_forward` method below after
        // liftover.
        let query_end = self
            .liftover(
                &reference
                    .end()
                    .clone()
                    .move_backward_checked_bounds(1, self.reference())
                    .unwrap(),
            )
            .unwrap()
            .move_forward(1)
            .unwrap();

        let query = Interval::try_new(query_start, query_end).unwrap();
        ContiguousIntervalPair::try_new(reference, query)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_it_constructs_a_valid_interval_pair() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "seq0:0-1000".parse::<Interval>()?;
        let query = "seq1:1000-2000".parse::<Interval>()?;

        ContiguousIntervalPair::try_new(reference, query)?;
        Ok(())
    }

    #[test]
    fn test_it_fails_when_the_interval_sizes_dont_match() -> Result<(), Box<dyn std::error::Error>>
    {
        let reference = "seq0:0-1000".parse::<Interval>()?;
        let query = "seq1:0-20000".parse::<Interval>()?;

        let err = ContiguousIntervalPair::try_new(reference, query).unwrap_err();
        assert!(matches!(err, Error::DistancesDontMatch(_, _)));
        assert_eq!(
            err.to_string(),
            "reference interval distance (1000) doesn't match query interval distance (20000)"
        );

        Ok(())
    }
}
