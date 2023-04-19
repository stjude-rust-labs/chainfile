//! Core functionality used across the crate.

pub mod coordinate;
pub mod interval;
pub mod strand;

pub use coordinate::Contig;
pub use coordinate::Coordinate;
pub use interval::Interval;
pub use strand::Strand;
