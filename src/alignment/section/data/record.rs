//! Accessories to alignment data records.

/// A kind of alignment data record.
///
/// In other words, whether the alignment data record is terminating or not.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Kind {
    /// Every line except the last line in an alignment section.
    NonTerminating,

    /// The last line in an alignment section.
    Terminating,
}
