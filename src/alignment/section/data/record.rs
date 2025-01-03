//! Alignment data records.

/// A kind of alignment data record.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A non-terminating alignment data line.
    ///
    /// In other words, all alignment data records that precede the last line
    /// within an alignment section.
    NonTerminating,

    /// The last line in an alignment section.
    Terminating,
}
