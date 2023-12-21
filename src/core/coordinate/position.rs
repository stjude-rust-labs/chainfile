//! Facilities for representing a zero-based position or the non-inclusive
//! negative bound position.

use std::ops::Add;
use std::ops::Sub;

/// An error related to a [`Position`].
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    /// Attempted to convert a [`Position::NegativeBound`] to a `usize`.
    CannotConvertNegativeBoundToUsize,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::CannotConvertNegativeBoundToUsize => {
                write!(f, "negative bound cannot be converted to a usize")
            }
        }
    }
}

impl std::error::Error for Error {}

/// The exact, 0-based position upon a contig which a coordinate is located.
///
/// In rare cases, you will need a value to represent the non-inclusive position
/// -1 (to indicate an interval on the negative strand that starts at an
/// arbitrary position and runs up until and includes zero). For those cases,
/// you can use [`Position::NegativeBound`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Position {
    /// A zero-based position.
    ZeroBased(usize),
    /// A special value representing the non-inclusive, -1 position.
    NegativeBound,
}

impl std::ops::Add for Position {
    type Output = usize;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Position::ZeroBased(a), Position::ZeroBased(b)) => a.add(b),
            (Position::ZeroBased(a), Position::NegativeBound) => {
                // Adding the negative bound from a zero-based position will end
                // up as `position + -1`, which is `position - 1`. However, if
                // position is zero, then the result would be -1, which is not
                // allowed for a usize. As such, we panic if that case is
                // encountered.
                match a {
                    0 => panic!("cannot add zero and negative bound"),
                    a => a.sub(1),
                }
            }
            (Position::NegativeBound, Position::ZeroBased(b)) => {
                // Adding a zero-based position to a negative bound will result
                // in `-1 + position` or `position - 1`. However, if position is
                // zero, then the result would be -1, which is not allowed for a
                // usize. As such, we panic if that case is encountered.
                match b {
                    0 => panic!("cannot add negative bound and zero"),
                    a => a.sub(1),
                }
            }
            (Position::NegativeBound, Position::NegativeBound) => {
                panic!("cannot add negative bound to negative bound")
            }
        }
    }
}

impl std::ops::Add for &Position {
    type Output = usize;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Position::ZeroBased(a), Position::ZeroBased(b)) => a.add(b),
            (Position::ZeroBased(a), Position::NegativeBound) => {
                // Adding the negative bound from a zero-based position will end
                // up as `position + -1`, which is `position - 1`. However, if
                // position is zero, then the result would be -1, which is not
                // allowed for a usize. As such, we panic if that case is
                // encountered.
                match a {
                    0 => panic!("cannot add zero and negative bound"),
                    a => a.sub(1),
                }
            }
            (Position::NegativeBound, Position::ZeroBased(b)) => {
                // Adding a zero-based position to a negative bound will result
                // in `-1 + position` or `position - 1`. However, if position is
                // zero, then the result would be -1, which is not allowed for a
                // usize. As such, we panic if that case is encountered.
                match b {
                    0 => panic!("cannot add negative bound and zero"),
                    a => a.sub(1),
                }
            }
            (Position::NegativeBound, Position::NegativeBound) => {
                panic!("cannot add negative bound to negative bound")
            }
        }
    }
}

impl std::ops::Sub for Position {
    type Output = usize;

    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Position::ZeroBased(a), Position::ZeroBased(b)) => a.sub(b),
            (Position::ZeroBased(a), Position::NegativeBound) => {
                // Subtracting the negative bound from a zero-based position
                // will end up as `position - -1`, which is `position + 1`.
                a.add(1)
            }
            (Position::NegativeBound, Position::ZeroBased(b)) => {
                panic!("cannot substract {} from negative bound", b)
            }
            (Position::NegativeBound, Position::NegativeBound) => 0,
        }
    }
}

impl std::ops::Sub for &Position {
    type Output = usize;

    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Position::ZeroBased(a), Position::ZeroBased(b)) => a.sub(b),
            (Position::ZeroBased(a), Position::NegativeBound) => {
                // Subtracting the negative bound from a zero-based position
                // will end up as `position - -1`, which is `position + 1`.
                a.add(1)
            }
            (Position::NegativeBound, Position::ZeroBased(b)) => {
                panic!("cannot substract {} from negative bound", b)
            }
            (Position::NegativeBound, Position::NegativeBound) => {
                panic!("cannot subtract negative bound from negative bound")
            }
        }
    }
}

impl PartialOrd for Position {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Position {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (Position::ZeroBased(a), Position::ZeroBased(b)) => a.cmp(b),
            (Position::ZeroBased(_), Position::NegativeBound) => std::cmp::Ordering::Greater,
            (Position::NegativeBound, Position::ZeroBased(_)) => std::cmp::Ordering::Less,
            (Position::NegativeBound, Position::NegativeBound) => std::cmp::Ordering::Equal,
        }
    }
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Position::ZeroBased(a) => write!(f, "{}", a),
            Position::NegativeBound => write!(f, "$"),
        }
    }
}

impl Position {
    /// Creates a new zero-based position with the given location.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Position;
    /// use chainfile as chain;
    ///
    /// let position = Position::new(0);
    /// assert!(matches!(position, Position::ZeroBased(0)));
    /// ```
    pub fn new(position: usize) -> Position {
        Position::ZeroBased(position)
    }

    /// Creates a special, negative bound position.
    ///
    /// # Examples
    ///
    /// ```
    /// use chain::core::Position;
    /// use chainfile as chain;
    ///
    /// let position = Position::negative_bound();
    /// assert!(matches!(position, Position::NegativeBound));
    /// ```
    pub fn negative_bound() -> Position {
        Position::NegativeBound
    }
}

impl From<usize> for Position {
    fn from(value: usize) -> Self {
        Position::ZeroBased(value)
    }
}

impl TryFrom<Position> for usize {
    type Error = Error;

    fn try_from(value: Position) -> Result<Self, Self::Error> {
        match value {
            Position::ZeroBased(v) => Ok(v),
            Position::NegativeBound => Err(Error::CannotConvertNegativeBoundToUsize),
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_zero_based_position_creation() {
        let position = Position::new(10);
        assert_eq!(
            <Position as TryInto<usize>>::try_into(position),
            Ok(10usize)
        );
    }

    #[test]
    fn test_negative_position_creation() {
        let position = Position::negative_bound();
        assert_eq!(
            <Position as TryInto<usize>>::try_into(position),
            Err(Error::CannotConvertNegativeBoundToUsize)
        );
    }

    #[test]
    fn test_adding_positions_together() {
        // Adding two zero-based positions
        let a = Position::new(0);
        let b = Position::new(10);
        assert_eq!(a + b, 10);

        // Adding a zero-based position and a negative bound
        let a = Position::new(0);
        let b = Position::negative_bound();
        let result = std::panic::catch_unwind(|| a + b);
        assert!(result.is_err());

        // Adding a negative bound to a zero-based position
        let a = Position::negative_bound();
        let b = Position::new(0);
        let result = std::panic::catch_unwind(|| a + b);
        assert!(result.is_err());

        // Adding a negative bound to negative bound
        let a = Position::negative_bound();
        let b = Position::negative_bound();
        let result = std::panic::catch_unwind(|| a + b);
        assert!(result.is_err());
    }

    #[test]
    fn test_subtracting_positions_from_one_another() {
        // Subtracting two zero-based positions
        let a = Position::new(10);
        let b = Position::new(5);
        assert_eq!(a - b, 5);

        // Subtracting with a zero-based position and a negative bound
        let a = Position::new(10);
        let b = Position::negative_bound();
        assert_eq!(a - b, 11);

        // Subtracting with a negative bound and a zero-based position
        let a = Position::negative_bound();
        let b = Position::new(10);
        let result = std::panic::catch_unwind(|| a - b);
        assert!(result.is_err());

        // Substracting a negative bound from a negative bound
        let a = Position::negative_bound();
        let b = Position::negative_bound();
        assert_eq!(a - b, 0);
    }

    #[test]
    fn test_ordering_of_positions() {
        // Ordering of two zero-based positions
        let a = Position::new(10);
        let b = Position::new(5);
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Greater);

        // Ordering of a zero-based position and a negative bound
        let a = Position::new(0);
        let b = Position::negative_bound();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Greater);

        // Ordering of a negative bound and a zero-based position
        let a = Position::negative_bound();
        let b = Position::new(0);
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Less);

        // Ordering of two negative bounds
        let a = Position::negative_bound();
        let b = Position::negative_bound();
        assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Equal));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Equal);
    }

    #[test]
    fn test_usize_into_position() {
        let position: Position = 30usize.into();
        assert_eq!(position, Position::ZeroBased(30));
    }
}
