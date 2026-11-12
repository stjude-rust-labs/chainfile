//! Sequence and coordinate helpers.

use omics::coordinate::Strand;
use omics::coordinate::interval::interbase::Interval;
use rand::Rng;
use rand::RngExt as _;

/// The nucleotides used for generated sequences.
pub const BASES: [u8; 4] = *b"ACGT";

/// Returns whether a base is one of `A`, `C`, `G`, or `T`.
pub fn is_nucleotide(base: u8) -> bool {
    BASES.contains(&base)
}

/// Gets the complement of an uppercase nucleotide.
///
/// Unrecognized symbols are returned unchanged.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Reverse complements a sequence.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().copied().map(complement).collect()
}

/// Orients forward-strand bases to a strand by reverse complementing them
/// when the strand is negative.
pub fn oriented(bases: Vec<u8>, strand: Strand) -> Vec<u8> {
    match strand {
        Strand::Positive => bases,
        Strand::Negative => reverse_complement(&bases),
    }
}

/// Gets the bounds of an interval in forward-strand numbering as a 0-based,
/// half-open `(start, end)` pair.
///
/// Intervals on the negative strand are stored with a start that is greater
/// than the end, so the bounds are swapped.
pub fn forward_bounds(interval: &Interval) -> (u64, u64) {
    let start = interval.start().position().get();
    let end = interval.end().position().get();

    match interval.strand() {
        Strand::Positive => (start, end),
        Strand::Negative => (end, start),
    }
}

/// Generates a random nucleotide.
pub fn random_base(rng: &mut impl Rng) -> u8 {
    BASES[rng.random_range(0..BASES.len())]
}

/// Generates a random nucleotide that differs from every base in `excluded`.
///
/// # Panics
///
/// Panics if `excluded` contains every nucleotide.
pub fn random_base_except(rng: &mut impl Rng, excluded: &[u8]) -> u8 {
    let choices = BASES
        .iter()
        .copied()
        .filter(|base| !excluded.contains(base))
        .collect::<Vec<_>>();

    choices[rng.random_range(0..choices.len())]
}

/// Generates a random sequence of nucleotides.
pub fn random_sequence(rng: &mut impl Rng, length: u64) -> Vec<u8> {
    (0..length).map(|_| random_base(rng)).collect()
}

#[cfg(test)]
mod tests {
    use omics::coordinate::Contig;
    use omics::coordinate::interbase::Coordinate;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    use super::*;

    #[test]
    fn reverse_complement_reverses_and_complements() {
        assert_eq!(reverse_complement(b"AACGTN"), b"NACGTT");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn forward_bounds_swaps_negative_strand_intervals() {
        let interval = |strand, start: u64, end: u64| {
            Interval::try_new(
                Coordinate::new(Contig::new_unchecked("seq0"), strand, start),
                Coordinate::new(Contig::new_unchecked("seq0"), strand, end),
            )
            .unwrap()
        };

        assert_eq!(forward_bounds(&interval(Strand::Positive, 2, 12)), (2, 12));
        assert_eq!(forward_bounds(&interval(Strand::Negative, 12, 2)), (2, 12));
    }

    #[test]
    fn random_base_except_avoids_excluded_bases() {
        let mut rng = StdRng::seed_from_u64(0);

        for _ in 0..100 {
            assert_eq!(random_base_except(&mut rng, b"ACG"), b'T');
            assert!(!b"AT".contains(&random_base_except(&mut rng, b"AT")));
        }
    }
}
