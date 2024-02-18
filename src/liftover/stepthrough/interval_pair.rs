//! Pairs of intervals that match contiguously.

use omics::coordinate::interval;
use omics::coordinate::interval::zero::Interval;
use omics::coordinate::zero::Coordinate;

/// An error related to constructing a contiguous interval pair.
#[derive(Debug)]
pub enum Error {
    /// The two intervals don't have the same size. As such, they can't
    /// contiguously map to one another.
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
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::zero::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
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
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::zero::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
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
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::zero::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
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
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::Strand;
    /// use omics::coordinate::interval::zero::Interval;
    /// use omics::coordinate::zero::Coordinate;
    ///
    /// // Positive-stranded to positive-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:1000-2000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let old = Coordinate::try_new("seq0", Strand::Positive, 50)?;
    /// let new = Coordinate::try_new("seq1", Strand::Positive, 1050)?;
    /// let lifted = pair.liftover(&old).unwrap().unwrap();
    ///
    /// assert_eq!(new, lifted);
    ///
    /// // Positive-stranded to negative-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:-:2000-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let old = Coordinate::try_new("seq0", Strand::Positive, 50)?;
    /// let new = Coordinate::try_new("seq1", Strand::Negative, 1950)?;
    /// let lifted = pair.liftover(&old).unwrap().unwrap();
    ///
    /// assert_eq!(new, lifted);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn liftover(&self, coordinate: &Coordinate) -> Option<Result<Coordinate, Error>> {
        let offset = self.reference().offset(coordinate)?;

        self.query()
            .translate(offset)
            .map_err(Error::InvalidInterval)
            .transpose()
    }

    /// Consumes self to clamp an interval pair to the specified `interval` for
    /// the reference interval. The query interval is similarly lifted over and
    /// clamped.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::Strand;
    /// use omics::coordinate::interval::zero::Interval;
    /// use omics::coordinate::zero::Coordinate;
    ///
    /// // Positive-stranded to positive-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:1000-2000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let interval = "seq0:+:50-51".parse::<Interval>()?;
    /// let result = pair.clamp(&interval)?;
    ///
    /// assert_eq!(
    ///     result.reference().start(),
    ///     &Coordinate::try_new("seq0", Strand::Positive, 50)?
    /// );
    /// assert_eq!(
    ///     result.reference().end(),
    ///     &Coordinate::try_new("seq0", Strand::Positive, 51)?
    /// );
    /// assert_eq!(
    ///     result.query().start(),
    ///     &Coordinate::try_new("seq1", Strand::Positive, 1050)?
    /// );
    /// assert_eq!(
    ///     result.query().end(),
    ///     &Coordinate::try_new("seq1", Strand::Positive, 1051)?
    /// );
    ///
    /// // Positive-stranded to negative-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:-:2000-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let interval = "seq0:+:50-51".parse::<Interval>()?;
    /// let result = pair.clamp(&interval)?;
    ///
    /// assert_eq!(
    ///     result.reference().start(),
    ///     &Coordinate::try_new("seq0", Strand::Positive, 50)?
    /// );
    /// assert_eq!(
    ///     result.reference().end(),
    ///     &Coordinate::try_new("seq0", Strand::Positive, 51)?
    /// );
    /// assert_eq!(
    ///     result.query().start(),
    ///     &Coordinate::try_new("seq1", Strand::Negative, 1950)?
    /// );
    /// assert_eq!(
    ///     result.query().end(),
    ///     &Coordinate::try_new("seq1", Strand::Negative, 1949)?
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn clamp(self, interval: &Interval) -> Result<ContiguousIntervalPair, Error> {
        let reference = self
            .reference()
            .clone()
            .clamp(interval)
            .map_err(Error::InvalidInterval)?;

        let query_start = self.liftover(reference.start()).unwrap().unwrap();

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
                    .unwrap()
                    .unwrap(),
            )
            .unwrap()
            .unwrap()
            .move_forward(1)
            .unwrap()
            .unwrap();

        let query = Interval::try_new(query_start, query_end).unwrap();
        ContiguousIntervalPair::try_new(reference, query)
    }
}

impl std::fmt::Display for ContiguousIntervalPair {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} -> {}", self.reference(), self.query())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_it_constructs_a_valid_interval_pair() -> Result<(), Box<dyn std::error::Error>> {
        let reference = "seq0:+:0-1000".parse::<Interval>()?;
        let query = "seq1:+:1000-2000".parse::<Interval>()?;

        ContiguousIntervalPair::try_new(reference, query)?;

        Ok(())
    }

    #[test]
    fn test_it_fails_when_the_interval_sizes_dont_match() -> Result<(), Box<dyn std::error::Error>>
    {
        let reference = "seq0:+:0-1000".parse::<Interval>()?;
        let query = "seq1:+:0-20000".parse::<Interval>()?;

        let err = ContiguousIntervalPair::try_new(reference, query).unwrap_err();
        assert!(matches!(err, Error::DistancesDontMatch(_, _)));
        assert_eq!(
            err.to_string(),
            "reference interval distance (1000) doesn't match query interval distance (20000)"
        );

        Ok(())
    }
}
