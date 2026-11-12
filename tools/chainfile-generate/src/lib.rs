//! Random test data for chain files.
//!
//! This crate generates:
//!
//! * [fixtures](fixture): a synthetic reference genome, query genome, and chain
//!   file between them, and
//! * random [BED](bed) and [VCF](vcf) records for the reference genome of a
//!   chain file, each generated to exercise a specific [`Situation`] (such as
//!   spanning a gap in a chain or falling in an unaligned region).
//!
//! Every generated record is named after its situation (e.g., `spans_gap_12`)
//! and is verified with a [`chainfile::liftover::Machine`] to have the
//! intended liftover outcome. Generation is deterministic for a given random
//! number generator state.

#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(rust_2021_compatibility)]
#![warn(missing_debug_implementations)]

pub mod bed;
pub mod fixture;
pub mod vcf;

mod fasta;
mod index;
mod io;
mod sequence;
mod situation;

#[cfg(test)]
mod tests;

pub use situation::SharedSituations;
pub use situation::Situation;
pub use situation::VcfSituations;
pub use situation::percentage;
