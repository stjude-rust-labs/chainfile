//! A machine for lifting over coordinates within a reference genome.

use std::collections::HashMap;

use rust_lapper as lapper;

use crate::core::Contig;
use crate::core::Interval;
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
    pub fn liftover(&self, interval: Interval) -> Option<Vec<ContiguousIntervalPair>> {
        let entry = self.inner.get(interval.contig())?;

        let results = entry
            .find(interval.start().position(), interval.end().position())
            .map(|e| e.val.clone())
            .map(|pair| pair.clamp(&interval).unwrap())
            .collect::<Vec<_>>();

        match results.is_empty() {
            true => None,
            false => Some(results),
        }
    }
}
