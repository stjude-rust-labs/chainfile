//! An index of the aligned and unaligned regions of a chain file.

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use chainfile::liftover::LiftoverResult;
use chainfile::liftover::Machine;
use chainfile::liftover::machine::builder::Builder;
use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;
use rand::Rng;
use rand::RngExt as _;
use rand::distr::Distribution as _;
use rand::distr::weighted::WeightedIndex;

use crate::io::read_to_end;
use crate::sequence::forward_bounds;
use crate::situation::Situation;

/// The number of attempts allowed per requested record before giving up.
const ATTEMPTS_PER_RECORD: usize = 1_000;

/// An aligned block of a chain.
#[derive(Clone, Debug)]
pub struct Block {
    /// The ID of the chain.
    pub chain: usize,

    /// The 0-based start of the block on the reference contig.
    pub start: u64,

    /// The 0-based, exclusive end of the block on the reference contig.
    pub end: u64,

    /// The query contig.
    pub query_contig: String,

    /// The query strand.
    pub strand: Strand,

    /// The 0-based start of the block on the forward strand of the query
    /// contig.
    pub query_start: u64,

    /// The 0-based, exclusive end of the block on the forward strand of the
    /// query contig.
    pub query_end: u64,
}

impl Block {
    /// Gets the forward-strand bounds on the query contig of the reference
    /// interval `[start, end)`, which must be within the block.
    ///
    /// On the negative strand, the reference offset `k` from the start of the
    /// block maps to the query position `query_end - 1 - k`.
    pub fn query_bounds(&self, start: u64, end: u64) -> (u64, u64) {
        debug_assert!(self.start <= start && start <= end && end <= self.end);

        match self.strand {
            Strand::Positive => (
                self.query_start + (start - self.start),
                self.query_start + (end - self.start),
            ),
            Strand::Negative => (
                self.query_end - (end - self.start),
                self.query_end - (start - self.start),
            ),
        }
    }
}

/// The chains that align a segment of a reference contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Coverage {
    /// No chain aligns the segment.
    Unaligned,

    /// Exactly one chain aligns the segment, with the block at the index.
    Clean(usize),

    /// Two or more chains align the segment.
    Overlap,
}

/// A maximal segment of a reference contig with the same coverage.
#[derive(Clone, Copy, Debug)]
struct Segment {
    /// The 0-based start of the segment.
    start: u64,

    /// The 0-based, exclusive end of the segment.
    end: u64,

    /// The chains that align the segment.
    coverage: Coverage,
}

/// A region of a reference contig.
#[derive(Clone, Copy, Debug)]
pub struct Region {
    /// The index of the contig.
    pub contig: usize,

    /// The 0-based start of the region.
    pub start: u64,

    /// The 0-based, exclusive end of the region.
    pub end: u64,

    /// The index of the block that aligns the region, if it is a clean piece of
    /// a single block.
    pub block: Option<usize>,
}

impl Region {
    /// Gets the length of the region.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
}

/// A pair of nearby blocks along with the window of reference intervals that
/// touch both of them.
///
/// Bridges describe both the gaps between consecutive blocks of a chain and
/// the junctions between blocks of different chains.
#[derive(Clone, Copy, Debug)]
pub struct Bridge {
    /// The index of the contig.
    pub contig: usize,

    /// The 0-based start of the window, which is within the block before the
    /// bridge.
    pub start: u64,

    /// The 0-based, exclusive end of the window, which is within the block
    /// after the bridge.
    pub end: u64,

    /// The 0-based start of the smallest interval that touches both blocks
    /// (the last base of the block before the bridge).
    pub core_start: u64,

    /// The 0-based, exclusive end of the smallest interval that touches both
    /// blocks (one past the first base of the block after the bridge).
    pub core_end: u64,
}

impl Bridge {
    /// Gets the minimum length of an interval that touches both blocks.
    pub fn min_length(&self) -> u64 {
        self.core_end - self.core_start
    }

    /// Gets the maximum length of an interval that touches both blocks and
    /// stays within the window.
    pub fn max_length(&self) -> u64 {
        self.end - self.start
    }

    /// Gets the lengths between `min` and `max` (inclusive) of intervals that
    /// touch both blocks and stay within the window, if there are any.
    pub fn lengths(&self, min: u64, max: u64) -> Option<(u64, u64)> {
        let min = min.max(self.min_length());
        let max = max.min(self.max_length());
        (min <= max).then_some((min, max))
    }

    /// Gets the range of starts of intervals of `length` bases that touch both
    /// blocks and stay within the window.
    ///
    /// The range is empty if no such interval exists.
    fn starts(&self, length: u64) -> (u64, u64) {
        // NOTE: the subtractions saturate because an interval near the start
        // of a contig may be longer than the position of the block after the
        // bridge.
        let lowest = self.start.max(self.core_end.saturating_sub(length));
        let highest = self.core_start.min(self.end.saturating_sub(length));
        (lowest, highest)
    }

    /// Gets the number of positions where an interval of `length` bases that
    /// touches both blocks and stays within the window can start.
    pub fn start_count(&self, length: u64) -> u64 {
        let (lowest, highest) = self.starts(length);
        (highest + 1).saturating_sub(lowest)
    }

    /// Samples an interval with a length between `min` and `max` (inclusive)
    /// that touches both blocks and stays within the window.
    ///
    /// Returns [`None`] if there is no such interval.
    pub fn sample(&self, rng: &mut impl Rng, min: u64, max: u64) -> Option<(u64, u64)> {
        let (min, max) = self.lengths(min, max)?;
        let length = rng.random_range(min..=max);
        let (lowest, highest) = self.starts(length);

        if lowest > highest {
            return None;
        }

        let start = rng.random_range(lowest..=highest);
        Some((start, start + length))
    }

    /// Gets the weight of the bridge when sampling intervals with a length
    /// between `min` and `max` (inclusive), which is the number of positions
    /// where the longest such interval can start (or zero if there are no such
    /// intervals).
    pub fn weight(&self, min: u64, max: u64) -> u64 {
        self.lengths(min, max)
            .map_or(0, |(_, max)| self.start_count(max))
    }
}

/// A candidate region for records that two or more chains align.
#[derive(Clone, Copy, Debug)]
pub enum MultipleChains {
    /// A region that two or more chains align.
    Overlap(Region),

    /// The junction between blocks of two different chains.
    Junction(Bridge),
}

/// The index of a single reference contig.
#[derive(Debug)]
pub struct ContigIndex {
    /// The name of the contig.
    pub name: String,

    /// The size of the contig.
    pub size: u64,

    /// The aligned blocks on the contig.
    pub blocks: Vec<Block>,

    /// The segments of the contig, ordered by position.
    segments: Vec<Segment>,

    /// The consecutive blocks within a chain that are separated by a gap.
    gaps: Vec<(usize, usize)>,
}

impl ContigIndex {
    /// Creates a contig index.
    fn new(name: String, size: u64, blocks: Vec<Block>, gaps: Vec<(usize, usize)>) -> Self {
        let segments = segments(size, &blocks);

        Self {
            name,
            size,
            blocks,
            segments,
            gaps,
        }
    }

    /// Returns whether a segment is aligned by no chain other than `chain`.
    fn is_free_for(&self, segment: &Segment, chain: usize) -> bool {
        match segment.coverage {
            Coverage::Unaligned => true,
            Coverage::Clean(block) => self.blocks[block].chain == chain,
            Coverage::Overlap => false,
        }
    }

    /// Gets the window around the gap between two consecutive blocks that no
    /// other chain aligns, if the gap itself is not aligned by another chain.
    fn gap_window(&self, contig: usize, before: usize, after: usize) -> Option<Bridge> {
        let before = &self.blocks[before];
        let after = &self.blocks[after];
        let chain = before.chain;

        let core_start = before.end - 1;
        let core_end = after.start + 1;

        let first = self
            .segments
            .partition_point(|segment| segment.end <= core_start);
        let mut last = first;

        loop {
            let segment = self.segments.get(last)?;

            if !self.is_free_for(segment, chain) {
                return None;
            }

            if segment.end >= core_end {
                break;
            }

            last += 1;
        }

        let mut lowest = first;

        while self.segments[lowest].start > before.start
            && lowest > 0
            && self.is_free_for(&self.segments[lowest - 1], chain)
        {
            lowest -= 1;
        }

        let mut highest = last;

        while self.segments[highest].end < after.end
            && highest + 1 < self.segments.len()
            && self.is_free_for(&self.segments[highest + 1], chain)
        {
            highest += 1;
        }

        Some(Bridge {
            contig,
            start: self.segments[lowest].start.max(before.start),
            end: self.segments[highest].end.min(after.end),
            core_start,
            core_end,
        })
    }
}

impl ContigIndex {
    /// Gets the junctions between blocks of different chains.
    ///
    /// A junction is a pair of blocks from different chains that are
    /// consecutive when the blocks are sorted by position, where the second
    /// block starts at or after the end of the first block. (Blocks that
    /// overlap are instead covered by the overlap regions.)
    fn junctions(&self, contig: usize) -> Vec<Bridge> {
        let mut order = (0..self.blocks.len()).collect::<Vec<_>>();
        order.sort_by_key(|&index| (self.blocks[index].start, self.blocks[index].end));

        order
            .windows(2)
            .filter_map(|pair| {
                let before = &self.blocks[pair[0]];
                let after = &self.blocks[pair[1]];

                (before.chain != after.chain && after.start >= before.end).then(|| Bridge {
                    contig,
                    start: before.start,
                    end: after.end,
                    core_start: before.end - 1,
                    core_end: after.start + 1,
                })
            })
            .collect()
    }
}

/// Splits a contig into maximal segments with the same coverage.
fn segments(size: u64, blocks: &[Block]) -> Vec<Segment> {
    let mut events = Vec::with_capacity(blocks.len() * 2);

    for (index, block) in blocks.iter().enumerate() {
        // NOTE: ends sort before starts at the same position (`false` before
        // `true`) so that abutting blocks never appear to overlap.
        events.push((block.start, true, index));
        events.push((block.end, false, index));
    }

    events.sort_unstable();

    let mut result = Vec::<Segment>::new();
    let mut active = BTreeSet::new();
    let mut cursor = 0;
    let mut events = events.into_iter().peekable();

    /// Appends a segment, merging it with the previous segment if they have
    /// the same coverage.
    fn push(result: &mut Vec<Segment>, start: u64, end: u64, active: &BTreeSet<usize>) {
        let coverage = match active.len() {
            0 => Coverage::Unaligned,
            1 => Coverage::Clean(*active.first().unwrap()),
            _ => Coverage::Overlap,
        };

        match result.last_mut() {
            Some(last) if last.coverage == coverage && last.end == start => last.end = end,
            _ => result.push(Segment {
                start,
                end,
                coverage,
            }),
        }
    }

    while let Some(&(position, ..)) = events.peek() {
        if position > cursor {
            push(&mut result, cursor, position.min(size), &active);
            cursor = position;
        }

        while let Some((_, is_start, index)) = events.next_if(|(p, ..)| *p == position) {
            if is_start {
                active.insert(index);
            } else {
                active.remove(&index);
            }
        }
    }

    if size > cursor {
        push(&mut result, cursor, size, &active);
    }

    result
}

/// An index of the aligned and unaligned regions of a chain file along with a
/// liftover machine to verify generated records.
#[derive(Debug)]
pub struct Index {
    /// The reference contigs, sorted by name.
    pub contigs: Vec<ContigIndex>,

    /// The liftover machine.
    machine: Machine,
}

impl Index {
    /// Loads the index from an (optionally gzip-compressed) chain file.
    pub fn load(path: &Path) -> Result<Self> {
        let data = read_to_end(path)?;
        Self::from_bytes(&data).with_context(|| format!("indexing `{}`", path.display()))
    }

    /// Builds the index from the contents of a chain file.
    fn from_bytes(data: &[u8]) -> Result<Self> {
        let machine = Builder
            .try_build_from(chainfile::Reader::new(data))
            .context("building the liftover machine")?;

        let mut contigs = BTreeMap::<String, (u64, Vec<Block>, Vec<(usize, usize)>)>::new();

        for result in chainfile::Reader::new(data).sections() {
            let section = result.context("reading a chain")?;
            let header = section.header();
            let reference = header.reference_sequence();
            let query = header.query_sequence();

            if reference.strand() != Strand::Positive {
                bail!(
                    "chain {} is on the negative strand of the reference genome, which is not \
                     supported",
                    header.id()
                );
            }

            let (_, blocks, gaps) = contigs
                .entry(reference.chromosome_name().to_string())
                .or_insert_with(|| (reference.chromosome_size(), Vec::new(), Vec::new()));

            let mut previous = None::<usize>;

            for result in section
                .stepthrough_with_data()
                .with_context(|| format!("stepping through chain {}", header.id()))?
            {
                let (pair, _) =
                    result.with_context(|| format!("stepping through chain {}", header.id()))?;

                let (start, end) = forward_bounds(pair.reference());
                let (query_start, query_end) = forward_bounds(pair.query());

                if start == end {
                    continue;
                }

                let block = Block {
                    chain: header.id(),
                    start,
                    end,
                    query_contig: query.chromosome_name().to_string(),
                    strand: pair.query().strand(),
                    query_start,
                    query_end,
                };

                if let Some(previous) = previous {
                    let before: &Block = &blocks[previous];
                    let contiguous = match block.strand {
                        Strand::Positive => before.query_end == block.query_start,
                        Strand::Negative => before.query_start == block.query_end,
                    };

                    if before.end != block.start || !contiguous {
                        gaps.push((previous, blocks.len()));
                    }
                }

                previous = Some(blocks.len());
                blocks.push(block);
            }
        }

        let contigs = contigs
            .into_iter()
            .map(|(name, (size, blocks, gaps))| ContigIndex::new(name, size, blocks, gaps))
            .collect();

        Ok(Self { contigs, machine })
    }

    /// Gets the clean pieces of blocks (those aligned by exactly one chain),
    /// optionally restricted to blocks on one query strand.
    pub fn clean_pieces(&self, strand: Option<Strand>) -> Vec<Region> {
        self.regions(|contig, coverage| match coverage {
            Coverage::Clean(block) => {
                let matches = strand.is_none_or(|strand| contig.blocks[block].strand == strand);
                matches.then_some(Some(block))
            }
            _ => None,
        })
    }

    /// Gets the regions that no chain aligns.
    pub fn unaligned(&self) -> Vec<Region> {
        self.regions(|_, coverage| (coverage == Coverage::Unaligned).then_some(None))
    }

    /// Gets the regions that two or more chains align.
    pub fn overlaps(&self) -> Vec<Region> {
        self.regions(|_, coverage| (coverage == Coverage::Overlap).then_some(None))
    }

    /// Gets the regions of segments selected by `select`, which returns the
    /// block for the region (if any) when the segment is selected.
    fn regions(
        &self,
        select: impl Fn(&ContigIndex, Coverage) -> Option<Option<usize>>,
    ) -> Vec<Region> {
        self.contigs
            .iter()
            .enumerate()
            .flat_map(|(index, contig)| {
                let select = &select;
                contig.segments.iter().filter_map(move |segment| {
                    select(contig, segment.coverage).map(|block| Region {
                        contig: index,
                        start: segment.start,
                        end: segment.end,
                        block,
                    })
                })
            })
            .collect()
    }

    /// Gets the gaps between consecutive blocks of a chain that no other chain
    /// aligns.
    pub fn gaps(&self) -> Vec<Bridge> {
        self.contigs
            .iter()
            .enumerate()
            .flat_map(|(index, contig)| {
                contig
                    .gaps
                    .iter()
                    .filter_map(move |&(before, after)| contig.gap_window(index, before, after))
            })
            .collect()
    }

    /// Gets the junctions between blocks of different chains.
    pub fn junctions(&self) -> Vec<Bridge> {
        self.contigs
            .iter()
            .enumerate()
            .flat_map(|(index, contig)| contig.junctions(index))
            .collect()
    }

    /// Gets the candidate regions for records that two or more chains align:
    /// the overlap regions and the junctions between blocks of different
    /// chains.
    pub fn multiple_chains(&self) -> Vec<MultipleChains> {
        self.overlaps()
            .into_iter()
            .map(MultipleChains::Overlap)
            .chain(self.junctions().into_iter().map(MultipleChains::Junction))
            .collect()
    }

    /// Gets a block.
    pub fn block(&self, contig: usize, block: usize) -> &Block {
        &self.contigs[contig].blocks[block]
    }

    /// Lifts over the interval `[start, end)` on the positive strand of a
    /// reference contig.
    pub fn liftover(&self, contig: usize, start: u64, end: u64) -> Option<Vec<LiftoverResult>> {
        let contig = Contig::new_unchecked(self.contigs[contig].name.as_str());

        let interval = Interval::try_new(
            Coordinate::new(contig.clone(), Strand::Positive, start),
            Coordinate::new(contig, Strand::Positive, end),
        )
        .expect("interval should be valid");

        self.machine.liftover(interval)
    }
}

/// Gets the single segment that the interval `[start, end)` maps to if the
/// interval maps to exactly one chain with exactly one segment that covers
/// the whole interval.
pub fn single_segment(
    results: &Option<Vec<LiftoverResult>>,
    start: u64,
    end: u64,
) -> Option<&chainfile::liftover::stepthrough::interval_pair::ContiguousIntervalPair> {
    let [result] = results.as_deref()? else {
        return None;
    };

    let [segment] = result.segments() else {
        return None;
    };

    (forward_bounds(segment.reference()) == (start, end)).then_some(segment)
}

/// Candidate regions sampled in proportion to their weights.
pub struct Candidates<T> {
    /// The candidates.
    items: Vec<T>,

    /// The distribution of candidate indices.
    distribution: Option<WeightedIndex<u64>>,
}

impl<T> Candidates<T> {
    /// Creates candidates, dropping any with a weight of zero.
    pub fn new(items: impl IntoIterator<Item = T>, weight: impl Fn(&T) -> u64) -> Self {
        let (items, weights): (Vec<_>, Vec<_>) = items
            .into_iter()
            .map(|item| {
                let weight = weight(&item);
                (item, weight)
            })
            .filter(|(_, weight)| *weight > 0)
            .unzip();

        let distribution = WeightedIndex::new(weights).ok();

        Self {
            items,
            distribution,
        }
    }

    /// Returns whether there are no candidates.
    pub fn is_empty(&self) -> bool {
        self.distribution.is_none()
    }

    /// Samples a candidate.
    ///
    /// # Panics
    ///
    /// Panics if there are no candidates.
    pub fn sample(&self, rng: &mut impl Rng) -> &T {
        let distribution = self
            .distribution
            .as_ref()
            .expect("there should be candidates");
        &self.items[distribution.sample(rng)]
    }
}

/// Generates `count` records for a situation.
///
/// Each call to `attempt` tries to generate one record from a sampled
/// candidate and returns [`None`] if the record does not have the intended
/// outcome.
pub fn generate_records<T, C, R: Rng>(
    situation: Situation,
    count: usize,
    candidates: &Candidates<C>,
    rng: &mut R,
    mut attempt: impl FnMut(&mut R, &C, String) -> Result<Option<T>>,
) -> Result<Vec<T>> {
    let mut records = Vec::with_capacity(count);

    if count == 0 {
        return Ok(records);
    }

    if candidates.is_empty() {
        bail!("the chain file has no {situation} candidates");
    }

    let attempts = count.saturating_mul(ATTEMPTS_PER_RECORD);

    for _ in 0..attempts {
        let candidate = candidates.sample(rng);

        if let Some(record) = attempt(rng, candidate, situation.record_name(records.len()))? {
            records.push(record);

            if records.len() == count {
                return Ok(records);
            }
        }
    }

    bail!(
        "could only generate {} of {count} {situation} records after {attempts} attempts",
        records.len()
    )
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng as _;

    use super::*;

    /// A chain file with a gap, an overlap, and a negative strand chain.
    fn index() -> Index {
        let data = [
            "chain 100 chr1 100 + 10 40 chrA 50 + 0 27 1",
            "10\t5\t2",
            "15",
            "",
            "chain 50 chr1 100 + 35 60 chrB 30 - 5 30 2",
            "25",
            "",
            "chain 10 chr1 100 + 80 90 chrB 30 + 0 10 3",
            "10",
            "",
        ]
        .join("\n");

        Index::from_bytes(data.as_bytes()).unwrap()
    }

    /// Gets the bounds of regions.
    fn bounds(regions: &[Region]) -> Vec<(u64, u64)> {
        regions.iter().map(|r| (r.start, r.end)).collect()
    }

    #[test]
    fn regions_are_classified() {
        let index = index();

        assert_eq!(
            bounds(&index.clean_pieces(Some(Strand::Positive))),
            [(10, 20), (25, 35), (80, 90)]
        );
        assert_eq!(
            bounds(&index.clean_pieces(Some(Strand::Negative))),
            [(40, 60)]
        );
        assert_eq!(bounds(&index.overlaps()), [(35, 40)]);
        assert_eq!(
            bounds(&index.unaligned()),
            [(0, 10), (20, 25), (60, 80), (90, 100)]
        );
    }

    #[test]
    fn blocks_map_to_the_query() {
        let index = index();
        let blocks = &index.contigs[0].blocks;

        assert_eq!(blocks[0].query_bounds(10, 20), (0, 10));
        assert_eq!(blocks[1].query_bounds(25, 35), (12, 22));
        assert_eq!(blocks[2].query_bounds(35, 36), (24, 25));
        assert_eq!(blocks[2].query_bounds(59, 60), (0, 1));
    }

    #[test]
    fn gaps_have_windows() {
        let index = index();
        let gaps = index.gaps();

        assert_eq!(gaps.len(), 1);
        assert_eq!((gaps[0].start, gaps[0].end), (10, 35));
        assert_eq!((gaps[0].core_start, gaps[0].core_end), (19, 26));
        assert_eq!(gaps[0].min_length(), 7);
    }

    #[test]
    fn junctions_join_different_chains() {
        let index = index();
        let junctions = index.junctions();

        // NOTE: the second block of chain 1 overlaps chain 2, so only the
        // junction between chains 2 and 3 is reported.
        assert_eq!(junctions.len(), 1);
        assert_eq!((junctions[0].start, junctions[0].end), (35, 90));
        assert_eq!((junctions[0].core_start, junctions[0].core_end), (59, 81));
        assert_eq!(index.multiple_chains().len(), 2);
    }

    #[test]
    fn bridges_sample_only_feasible_intervals() {
        // NOTE: the window around the gap starts at the start of the contig and
        // is much longer than the position of the block after the gap, which
        // used to underflow when computing the lowest start.
        let data = [
            "chain 100 chr1 1000 + 0 500 chrA 1000 + 0 499 1",
            "2\t1\t0",
            "497",
            "",
        ]
        .join("\n");
        let index = Index::from_bytes(data.as_bytes()).unwrap();
        let [gap] = index.gaps()[..] else {
            panic!("there should be one gap");
        };

        assert_eq!(
            (gap.start, gap.end, gap.core_start, gap.core_end),
            (0, 500, 1, 4)
        );
        assert_eq!(gap.lengths(1, 2), None);
        assert_eq!(gap.weight(1, 2), 0);
        assert_eq!(gap.lengths(1, 1_000), Some((3, 500)));
        assert_eq!(gap.start_count(3), 1);
        assert_eq!(gap.start_count(100), 2);
        assert_eq!(gap.start_count(500), 1);

        let mut rng = rand::rngs::StdRng::seed_from_u64(0);

        for _ in 0..1_000 {
            let (start, end) = gap.sample(&mut rng, 1, 1_000).unwrap();
            assert!(start <= 1 && (4..=500).contains(&end));
        }
    }
}
