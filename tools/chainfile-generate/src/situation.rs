//! The situations that generated records exercise.

use std::fmt;

use anyhow::Result;
use anyhow::bail;

/// The tolerance used when checking that percentages sum to at most 100.
const PERCENTAGE_TOLERANCE: f64 = 1e-9;

/// Parses a percentage between 0 and 100 (inclusive).
pub fn percentage(value: &str) -> std::result::Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("`{value}` is not a number"))?;

    if !(0.0..=100.0).contains(&parsed) {
        return Err(format!("`{value}` is not a percentage between 0 and 100"));
    }

    Ok(parsed)
}

/// The situation a generated record was generated to exercise.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Situation {
    /// The record is contained within a single aligned block of a chain on the
    /// positive query strand.
    Contained,

    /// The record spans a gap between two aligned blocks of a chain.
    SpansGap,

    /// The record is in a reference region that no chain aligns.
    Unaligned,

    /// The record is aligned by two or more chains (such as a record that
    /// overlaps a duplicated region or spans the junction between two chains).
    MultipleChains,

    /// The record is contained within a single aligned block of a chain on the
    /// negative query strand.
    OppositeStrand,

    /// The record is an insertion or deletion.
    Indel,

    /// The record's alternate allele matches the query genome, and its
    /// reference allele does not.
    SwappedAlleles,

    /// Neither of the record's alleles matches the query genome.
    ReferenceMismatch,

    /// The record has two alternate alleles.
    Multiallelic,
}

impl Situation {
    /// All situations in declaration order.
    pub const ALL: [Situation; 9] = [
        Situation::Contained,
        Situation::SpansGap,
        Situation::Unaligned,
        Situation::MultipleChains,
        Situation::OppositeStrand,
        Situation::Indel,
        Situation::SwappedAlleles,
        Situation::ReferenceMismatch,
        Situation::Multiallelic,
    ];

    /// Gets the snake case name of the situation.
    pub fn name(&self) -> &'static str {
        match self {
            Situation::Contained => "contained",
            Situation::SpansGap => "spans_gap",
            Situation::Unaligned => "unaligned",
            Situation::MultipleChains => "multiple_chains",
            Situation::OppositeStrand => "opposite_strand",
            Situation::Indel => "indel",
            Situation::SwappedAlleles => "swapped_alleles",
            Situation::ReferenceMismatch => "reference_mismatch",
            Situation::Multiallelic => "multiallelic",
        }
    }

    /// Parses the situation from a generated record name such as
    /// `spans_gap_12`.
    ///
    /// Returns [`None`] if the name is not a generated record name.
    pub fn from_record_name(name: &str) -> Option<Situation> {
        let (situation, index) = name.rsplit_once('_')?;

        if index.is_empty() || !index.bytes().all(|byte| byte.is_ascii_digit()) {
            return None;
        }

        Situation::ALL
            .into_iter()
            .find(|candidate| candidate.name() == situation)
    }

    /// Gets the name of the `index`th record generated for the situation.
    pub(crate) fn record_name(&self, index: usize) -> String {
        format!("{}_{index}", self.name())
    }
}

impl fmt::Display for Situation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// The percentages of records generated for situations that apply to every
/// output format.
///
/// The remaining records are contained within a single aligned block of a
/// chain on the positive query strand.
#[derive(Clone, Debug, Default, clap::Args)]
pub struct SharedSituations {
    /// The percentage of records that span a gap between two aligned blocks of
    /// a chain.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub spans_gap: f64,

    /// The percentage of records in reference regions that no chain aligns.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub unaligned: f64,

    /// The percentage of records that two or more chains align (such as
    /// records that overlap duplicated regions or span the junction between
    /// two chains).
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub multiple_chains: f64,

    /// The percentage of records within a single aligned block of a chain on
    /// the negative query strand.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub opposite_strand: f64,
}

impl SharedSituations {
    /// Gets the percentage of records for each situation.
    pub(crate) fn percentages(&self) -> [(Situation, f64); 4] {
        [
            (Situation::SpansGap, self.spans_gap),
            (Situation::Unaligned, self.unaligned),
            (Situation::MultipleChains, self.multiple_chains),
            (Situation::OppositeStrand, self.opposite_strand),
        ]
    }
}

/// The percentages of records generated for situations that only apply to
/// VCF output.
#[derive(Clone, Debug, Default, clap::Args)]
pub struct VcfSituations {
    /// The percentage of records that are insertions or deletions within a
    /// single aligned block of a chain.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub indels: f64,

    /// The percentage of records whose alternate allele (but not reference
    /// allele) matches the query genome.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub swapped_alleles: f64,

    /// The percentage of records where neither allele matches the query
    /// genome.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub reference_mismatch: f64,

    /// The percentage of records with two alternate alleles.
    #[arg(long, default_value_t = 0.0, value_parser = percentage)]
    pub multiallelic: f64,
}

impl VcfSituations {
    /// Gets the percentage of records for each situation.
    pub(crate) fn percentages(&self) -> [(Situation, f64); 4] {
        [
            (Situation::Indel, self.indels),
            (Situation::SwappedAlleles, self.swapped_alleles),
            (Situation::ReferenceMismatch, self.reference_mismatch),
            (Situation::Multiallelic, self.multiallelic),
        ]
    }
}

/// Splits `total` records between the situations.
///
/// Each situation receives `floor(total * percentage / 100)` records, and
/// the remaining records are [`Situation::Contained`]. The counts are
/// returned in the order of [`Situation::ALL`].
pub(crate) fn counts(
    total: usize,
    percentages: &[(Situation, f64)],
) -> Result<Vec<(Situation, usize)>> {
    let mut sum = 0.0;

    for &(situation, percentage) in percentages {
        if !(0.0..=100.0).contains(&percentage) {
            bail!("the percentage for `{situation}` ({percentage}) is not between 0 and 100");
        }

        sum += percentage;
    }

    if sum > 100.0 + PERCENTAGE_TOLERANCE {
        bail!("the situation percentages sum to {sum}, which is more than 100");
    }

    let mut counts = percentages
        .iter()
        .map(|&(situation, percentage)| {
            (
                situation,
                (total as f64 * percentage / 100.0).floor() as usize,
            )
        })
        .collect::<Vec<_>>();

    let assigned = counts.iter().map(|(_, count)| count).sum::<usize>();
    counts.push((Situation::Contained, total.saturating_sub(assigned)));
    counts.sort();

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip() {
        for situation in Situation::ALL {
            assert_eq!(situation.to_string(), situation.name());
            assert_eq!(
                Situation::from_record_name(&situation.record_name(12)),
                Some(situation)
            );
        }

        assert_eq!(Situation::from_record_name("spans_gap"), None);
        assert_eq!(Situation::from_record_name("spans_gap_"), None);
        assert_eq!(Situation::from_record_name("spans_gap_x"), None);
        assert_eq!(Situation::from_record_name("unknown_1"), None);
        assert_eq!(Situation::from_record_name("rs123"), None);
    }

    #[test]
    fn percentage_rejects_out_of_range_values() {
        assert_eq!(percentage("0"), Ok(0.0));
        assert_eq!(percentage("12.5"), Ok(12.5));
        assert_eq!(percentage("100"), Ok(100.0));
        assert!(percentage("-1").is_err());
        assert!(percentage("101").is_err());
        assert!(percentage("NaN").is_err());
        assert!(percentage("ten").is_err());
    }

    #[test]
    fn counts_split_records() {
        let counts = counts(
            1_000,
            &[(Situation::SpansGap, 8.0), (Situation::Unaligned, 0.15)],
        )
        .unwrap();

        assert_eq!(
            counts,
            [
                (Situation::Contained, 919),
                (Situation::SpansGap, 80),
                (Situation::Unaligned, 1),
            ]
        );
    }

    #[test]
    fn counts_reject_percentages_over_100() {
        assert!(
            counts(
                10,
                &[(Situation::SpansGap, 60.0), (Situation::Unaligned, 41.0)]
            )
            .is_err()
        );
        assert!(counts(10, &[(Situation::SpansGap, 101.0)]).is_err());
        assert!(
            counts(
                10,
                &[(Situation::SpansGap, 50.0), (Situation::Unaligned, 50.0)]
            )
            .is_ok()
        );
    }
}
