//! Sequence and coordinate helpers.

use omics::coordinate::Strand;
use omics::coordinate::interval::interbase::Interval;

/// Gets the complement of a nucleotide (including IUPAC ambiguity codes).
///
/// Case is preserved, and unrecognized symbols are returned unchanged.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' | b'u' => b'a',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        other => other,
    }
}

/// Reverse complements a sequence.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().copied().map(complement).collect()
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

/// Orients forward-strand bases to a strand by reverse complementing them
/// when the strand is negative.
pub fn oriented(bases: Vec<u8>, strand: Strand) -> Vec<u8> {
    match strand {
        Strand::Positive => bases,
        Strand::Negative => reverse_complement(&bases),
    }
}

#[cfg(test)]
mod tests {
    use omics::coordinate::Contig;
    use omics::coordinate::interbase::Coordinate;

    use super::*;

    #[test]
    fn reverse_complement_handles_iupac_and_case() {
        assert_eq!(reverse_complement(b"ACGTN"), b"NACGT");
        assert_eq!(reverse_complement(b"acgtRYKM*"), b"*KMRYacgt");
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
}
