//! Pairs of intervals that match contiguously.

use omics::coordinate;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;

/// An error related to constructing a contiguous interval pair.
#[derive(Debug)]
pub enum Error {
    /// The two intervals don't have the same size. As such, they can't
    /// contiguously map to one another.
    EntityCountsDontMatch(u64, u64),

    /// An interval error.
    Interval(coordinate::interval::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::EntityCountsDontMatch(reference, query) => write!(
                f,
                "reference interval entity count ({reference}) doesn't match query interval entity count \
                 ({query})"
            ),
            Error::Interval(err) => write!(f, "interval error: {err}"),
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
    /// use omics::coordinate::interval::interbase::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
    /// ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn try_new(reference: Interval, query: Interval) -> Result<Self, Error> {
        if reference.count_entities() != query.count_entities() {
            return Err(Error::EntityCountsDontMatch(
                reference.count_entities(),
                query.count_entities(),
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
    /// use omics::coordinate::interval::interbase::Interval;
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

    /// Consumes `self` and gets the reference interval for the interval pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::interbase::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference.clone(), query)?;
    ///
    /// assert_eq!(pair.into_reference(), reference);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_reference(self) -> Interval {
        self.0
    }

    /// Gets the query interval by reference for the interval pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::interbase::Interval;
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

    /// Consumes `self` and returns the query interval for the interval pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::interbase::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query.clone())?;
    ///
    /// assert_eq!(pair.into_query(), query);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_query(self) -> Interval {
        self.1
    }

    /// Consumes `self` and returns the consituent contiguous intervals that
    /// make up the [`ContiguousIntervalPair`].
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::interval::interbase::Interval;
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:0-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference.clone(), query.clone())?;
    ///
    /// let (a, b) = pair.into_parts();
    /// assert_eq!(a, reference);
    /// assert_eq!(b, query);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn into_parts(self) -> (Interval, Interval) {
        (self.0, self.1)
    }

    /// Lifts over a coordinate within the reference coordinate system to the
    /// query coordinate system.
    ///
    /// # Examples
    ///
    /// ```
    /// use chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair;
    /// use omics::coordinate::Strand;
    /// use omics::coordinate::interbase::Coordinate;
    /// use omics::coordinate::interval::interbase::Interval;
    ///
    /// // Positive-stranded to positive-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:1000-2000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let old = Coordinate::new("seq0", Strand::Positive, 50u64);
    /// let new = Coordinate::new("seq1", Strand::Positive, 1050u64);
    /// let lifted = pair.liftover(&old).unwrap();
    ///
    /// assert_eq!(new, lifted);
    ///
    /// // Positive-stranded to negative-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:-:2000-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let old = Coordinate::new("seq0", Strand::Positive, 50u64);
    /// let new = Coordinate::new("seq1", Strand::Negative, 1950u64);
    /// let lifted = pair.liftover(&old).unwrap();
    ///
    /// assert_eq!(new, lifted);
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn liftover(&self, coordinate: &Coordinate) -> Option<Coordinate> {
        let offset = self.reference().coordinate_offset(coordinate)?;
        self.query().coordinate_at_offset(offset)
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
    /// use omics::coordinate::interbase::Coordinate;
    /// use omics::coordinate::interval::interbase::Interval;
    ///
    /// // Positive-stranded to positive-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:+:1000-2000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let interval = "seq0:+:50-51".parse::<Interval>()?;
    /// let result = pair.clamp(interval)?;
    ///
    /// assert_eq!(
    ///     result.reference().start(),
    ///     &Coordinate::new("seq0", Strand::Positive, 50u64)
    /// );
    /// assert_eq!(
    ///     result.reference().end(),
    ///     &Coordinate::new("seq0", Strand::Positive, 51u64)
    /// );
    /// assert_eq!(
    ///     result.query().start(),
    ///     &Coordinate::new("seq1", Strand::Positive, 1050u64)
    /// );
    /// assert_eq!(
    ///     result.query().end(),
    ///     &Coordinate::new("seq1", Strand::Positive, 1051u64)
    /// );
    ///
    /// // Positive-stranded to negative-stranded
    ///
    /// let reference = "seq0:+:0-1000".parse::<Interval>()?;
    /// let query = "seq1:-:2000-1000".parse::<Interval>()?;
    /// let pair = ContiguousIntervalPair::try_new(reference, query)?;
    ///
    /// let interval = "seq0:+:50-51".parse::<Interval>()?;
    /// let result = pair.clamp(interval)?;
    ///
    /// assert_eq!(
    ///     result.reference().start(),
    ///     &Coordinate::new("seq0", Strand::Positive, 50u64)
    /// );
    /// assert_eq!(
    ///     result.reference().end(),
    ///     &Coordinate::new("seq0", Strand::Positive, 51u64)
    /// );
    /// assert_eq!(
    ///     result.query().start(),
    ///     &Coordinate::new("seq1", Strand::Negative, 1950u64)
    /// );
    /// assert_eq!(
    ///     result.query().end(),
    ///     &Coordinate::new("seq1", Strand::Negative, 1949u64)
    /// );
    ///
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn clamp(self, interval: Interval) -> Result<ContiguousIntervalPair, Error> {
        let reference = self
            .reference()
            .clone()
            .clamp(interval)
            .map_err(Error::Interval)?;

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
                    .move_backward(1)
                    .filter(|coord| self.reference().contains_coordinate(coord))
                    .unwrap(),
            )
            .unwrap()
            .move_forward(1)
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
    fn valid_interval_pair() {
        let reference = "seq0:+:0-1000".parse::<Interval>().unwrap();
        let query = "seq1:+:1000-2000".parse::<Interval>().unwrap();
        ContiguousIntervalPair::try_new(reference, query).unwrap();
    }

    #[test]
    fn interval_sizes_dont_match() {
        let reference = "seq0:+:0-1000".parse::<Interval>().unwrap();
        let query = "seq1:+:0-20000".parse::<Interval>().unwrap();

        let err = ContiguousIntervalPair::try_new(reference, query).unwrap_err();
        assert!(matches!(err, Error::EntityCountsDontMatch(_, _)));
        assert_eq!(
            err.to_string(),
            "reference interval entity count (1000) doesn't match query interval entity count \
             (20000)"
        );
    }
}
