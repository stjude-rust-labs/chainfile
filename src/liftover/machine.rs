//! A machine for lifting over coordinates within a reference genome.

use std::collections::HashMap;

use rust_lapper as lapper;

use crate::core::coordinate::Position;
use crate::core::Contig;
use crate::core::Interval;
use crate::core::Strand;
use crate::liftover::stepthrough::interval_pair::ContiguousIntervalPair;

pub mod builder;

pub use builder::Builder;

/// A machine for lifting over coordinates from a reference genome to a query
/// genome.
///
/// Generally, you will want to use a [`builder::Builder`] to construct one of
/// these.
pub struct Machine {
    inner: HashMap<Contig, lapper::Lapper<usize, ContiguousIntervalPair>>,
}

impl Machine {
    /// Performs a liftover from the specified `interval` to the query genome.
    pub fn liftover(&self, interval: &Interval) -> Option<Vec<ContiguousIntervalPair>> {
        let entry = self.inner.get(interval.contig())?;

        let (start, stop) = match interval.strand() {
            Strand::Positive => {
                let start = match interval.start().position() {
                    Position::ZeroBased(a) => *a,
                    Position::NegativeBound => {
                        unreachable!("negative bound not allowed on positive strand")
                    }
                };

                let stop = match interval.end().position() {
                    Position::ZeroBased(b) => *b,
                    Position::NegativeBound => {
                        unreachable!("negative bound not allowed on positive strand")
                    }
                };

                (start, stop)
            }
            Strand::Negative => {
                let start = match interval.end().position() {
                    Position::ZeroBased(a) => a + 1,
                    Position::NegativeBound => 0,
                };

                let stop = match interval.start().position() {
                    Position::ZeroBased(b) => b + 1,
                    Position::NegativeBound => unreachable!(
                        "negative bound will never be the start of a negative-stranded interval"
                    ),
                };

                (start, stop)
            }
        };

        let results = entry
            .find(start, stop)
            .map(|e| e.val.clone())
            .map(|pair| pair.clamp(interval).unwrap())
            .collect::<Vec<_>>();

        match results.is_empty() {
            true => None,
            false => Some(results),
        }
    }
}
