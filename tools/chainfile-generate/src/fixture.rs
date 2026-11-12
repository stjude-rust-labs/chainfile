//! Synthetic reference and query genomes with a chain file between them.
//!
//! The reference genome is random, and the query genome is assembled from
//! copies of reference regions (with substitutions, insertions, and
//! deletions) separated by random filler. The chain file describes exactly
//! how the copies were made, so it exercises the following:
//!
//! * chains with multiple blocks separated by reference-only, query-only, and
//!   double-sided gaps,
//! * unaligned regions in both genomes,
//! * chains on both query strands and on different query contigs,
//! * reference regions aligned by more than one chain (duplications), and
//! * substitutions, so the reference and query bases sometimes differ.

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::ops::RangeInclusive;
use std::path::Path;
use std::path::PathBuf;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use omics::coordinate::Strand;
use rand::Rng;
use rand::RngExt as _;
use rand::seq::SliceRandom as _;

use crate::sequence::random_base;
use crate::sequence::random_base_except;
use crate::sequence::random_sequence;
use crate::sequence::reverse_complement;

/// The name of the reference FASTA file.
const REFERENCE_FASTA: &str = "reference.fa";

/// The name of the query FASTA file.
const QUERY_FASTA: &str = "query.fa";

/// The name of the chain file.
const CHAIN: &str = "reference-to-query.chain";

/// The lengths of unaligned stretches of the reference genome and of random
/// filler in the query genome.
const UNALIGNED_LENGTH: RangeInclusive<u64> = 100..=5_000;

/// The lengths of reference regions aligned by a chain.
const CHAIN_LENGTH: RangeInclusive<u64> = 5_000..=50_000;

/// The minimum length of the reference region aligned by a chain.
const MIN_CHAIN_LENGTH: u64 = 500;

/// The maximum length of the reference region aligned by a duplicate chain.
const MAX_DUPLICATE_LENGTH: u64 = 10_000;

/// The sizes of aligned blocks.
const BLOCK_SIZE: RangeInclusive<u64> = 200..=5_000;

/// The sizes of reference-only and query-only gaps.
const SINGLE_SIDED_GAP: RangeInclusive<u64> = 1..=20;

/// The sizes of each side of double-sided gaps.
const DOUBLE_SIDED_GAP: RangeInclusive<u64> = 1..=50;

/// The number of bases on each FASTA sequence line.
const LINE_BASES: usize = 60;

/// Options for generating a fixture.
#[derive(Clone, Debug, clap::Args)]
pub struct Options {
    /// The number of contigs in each genome.
    #[arg(long, default_value_t = 4)]
    pub contigs: usize,

    /// The length of each reference contig.
    #[arg(long, default_value_t = 1_000_000)]
    pub contig_length: u64,

    /// The percentage of aligned reference bases that are substituted in the
    /// query genome.
    #[arg(long, default_value_t = 1.0, value_parser = crate::percentage)]
    pub substitution_pct: f64,

    /// The percentage of chains whose reference region is also aligned by a
    /// second chain.
    #[arg(long, default_value_t = 10.0, value_parser = crate::percentage)]
    pub duplication_pct: f64,

    /// The percentage of chains that align to the negative strand of the query
    /// genome.
    #[arg(long, default_value_t = 30.0, value_parser = crate::percentage)]
    pub negative_strand_pct: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            contigs: 4,
            contig_length: 1_000_000,
            substitution_pct: 1.0,
            duplication_pct: 10.0,
            negative_strand_pct: 30.0,
        }
    }
}

/// Paths of a generated fixture.
#[derive(Clone, Debug)]
pub struct Fixture {
    /// The chain file from the reference genome to the query genome.
    pub chain: PathBuf,

    /// The reference genome FASTA file (indexed).
    pub reference_fasta: PathBuf,

    /// The query genome FASTA file (indexed).
    pub query_fasta: PathBuf,
}

/// An aligned block of a chain followed by the gap to the next block.
#[derive(Clone, Copy, Debug)]
struct Block {
    /// The 0-based start of the block on the reference contig.
    start: u64,

    /// The size of the block.
    size: u64,

    /// The number of reference bases in the gap after the block.
    dt: u64,

    /// The number of query bases in the gap after the block.
    dq: u64,
}

/// A chain under construction.
#[derive(Debug)]
struct Chain {
    /// The index of the reference contig.
    reference_contig: usize,

    /// The blocks.
    blocks: Vec<Block>,

    /// The index of the query contig.
    query_contig: usize,

    /// The query strand.
    strand: Strand,

    /// The aligned query sequence in the orientation of the reference genome.
    sequence: Vec<u8>,

    /// The score.
    score: u64,

    /// The 0-based offset of the sequence on the forward strand of the query
    /// contig.
    offset: u64,
}

/// Generates a fixture.
///
/// Writes `reference.fa`, `reference.fa.fai`, `query.fa`, `query.fa.fai`, and
/// `reference-to-query.chain` into `directory` (which is created if it does
/// not exist).
pub fn generate(options: &Options, rng: &mut impl Rng, directory: &Path) -> Result<Fixture> {
    validate(options)?;

    let reference = (0..options.contigs)
        .map(|_| random_sequence(rng, options.contig_length))
        .collect::<Vec<_>>();

    let mut chains = Vec::new();

    for (contig, bases) in reference.iter().enumerate() {
        let length = bases.len() as u64;
        let mut cursor = 0;

        loop {
            cursor += rng.random_range(UNALIGNED_LENGTH);

            if cursor >= length {
                break;
            }

            let end = (cursor + rng.random_range(CHAIN_LENGTH)).min(length);

            if end - cursor < MIN_CHAIN_LENGTH {
                break;
            }

            let blocks = split_blocks(rng, cursor, end);
            let score = blocks.iter().map(|block| block.size).sum();
            chains.push(chain(rng, options, bases, contig, blocks, score));

            if rng.random_bool(options.duplication_pct / 100.0) {
                let length =
                    rng.random_range(MIN_CHAIN_LENGTH..=(end - cursor).min(MAX_DUPLICATE_LENGTH));
                let start = rng.random_range(cursor..=end - length);
                let blocks = split_blocks(rng, start, start + length);

                // NOTE: duplicates score lower than the primary chains so that
                // the primary chains are preferred.
                let score = blocks.iter().map(|block| block.size).sum::<u64>() / 2;
                chains.push(chain(rng, options, bases, contig, blocks, score));
            }

            cursor = end;
        }
    }

    let query = assemble(rng, options.contigs, &mut chains);

    std::fs::create_dir_all(directory)
        .with_context(|| format!("creating `{}`", directory.display()))?;

    let fixture = Fixture {
        chain: directory.join(CHAIN),
        reference_fasta: directory.join(REFERENCE_FASTA),
        query_fasta: directory.join(QUERY_FASTA),
    };

    write_fasta(&fixture.reference_fasta, &reference)?;
    write_fasta(&fixture.query_fasta, &query)?;
    write_chains(&fixture.chain, &chains, &reference, &query)?;

    Ok(fixture)
}

/// Validates the options.
fn validate(options: &Options) -> Result<()> {
    if options.contigs == 0 {
        bail!("the number of contigs must be at least one");
    }

    if options.contig_length == 0 {
        bail!("the contig length must be at least one");
    }

    for (name, value) in [
        ("substitution", options.substitution_pct),
        ("duplication", options.duplication_pct),
        ("negative strand", options.negative_strand_pct),
    ] {
        if !(0.0..=100.0).contains(&value) {
            bail!("the {name} percentage ({value}) is not between 0 and 100");
        }
    }

    Ok(())
}

/// Gets the name of a contig.
fn contig_name(index: usize) -> String {
    format!("chr{}", index + 1)
}

/// Splits the reference region `[start, end)` into blocks separated by gaps.
///
/// The region must be at least [`MIN_CHAIN_LENGTH`] bases long.
fn split_blocks(rng: &mut impl Rng, start: u64, end: u64) -> Vec<Block> {
    let max_gap = *DOUBLE_SIDED_GAP.end().max(SINGLE_SIDED_GAP.end());
    let min_block = *BLOCK_SIZE.start();

    let mut blocks = Vec::new();
    let mut cursor = start;

    loop {
        let remaining = end - cursor;
        let size = rng.random_range(BLOCK_SIZE);

        // NOTE: the last block absorbs the rest of the region so that the
        // chain ends at `end`. Another block is only started when there is
        // room for the largest gap and a full-sized block after it.
        if size + max_gap + min_block > remaining {
            blocks.push(Block {
                start: cursor,
                size: remaining,
                dt: 0,
                dq: 0,
            });

            return blocks;
        }

        let (dt, dq) = match rng.random_range(0..3) {
            0 => (rng.random_range(SINGLE_SIDED_GAP), 0),
            1 => (0, rng.random_range(SINGLE_SIDED_GAP)),
            _ => (
                rng.random_range(DOUBLE_SIDED_GAP),
                rng.random_range(DOUBLE_SIDED_GAP),
            ),
        };

        blocks.push(Block {
            start: cursor,
            size,
            dt,
            dq,
        });

        cursor += size + dt;
    }
}

/// Creates a chain from blocks on a reference contig, choosing its query
/// contig and strand and generating its query sequence.
fn chain(
    rng: &mut impl Rng,
    options: &Options,
    reference: &[u8],
    reference_contig: usize,
    blocks: Vec<Block>,
    score: u64,
) -> Chain {
    let query_contig = rng.random_range(0..options.contigs);
    let strand = if rng.random_bool(options.negative_strand_pct / 100.0) {
        Strand::Negative
    } else {
        Strand::Positive
    };

    let substitution = options.substitution_pct / 100.0;
    let mut sequence = Vec::new();

    for block in &blocks {
        let bases = &reference[block.start as usize..(block.start + block.size) as usize];

        for &base in bases {
            if rng.random_bool(substitution) {
                sequence.push(random_base_except(rng, &[base]));
            } else {
                sequence.push(base);
            }
        }

        sequence.extend((0..block.dq).map(|_| random_base(rng)));
    }

    Chain {
        reference_contig,
        blocks,
        query_contig,
        strand,
        sequence,
        score,
        offset: 0,
    }
}

/// Assembles the query contigs from the chains' query sequences separated by
/// random filler, recording the offset of each chain.
fn assemble(rng: &mut impl Rng, contigs: usize, chains: &mut [Chain]) -> Vec<Vec<u8>> {
    let mut pieces = vec![Vec::new(); contigs];

    for (index, chain) in chains.iter().enumerate() {
        pieces[chain.query_contig].push(index);
    }

    pieces
        .into_iter()
        .map(|mut pieces| {
            pieces.shuffle(rng);

            let length = rng.random_range(UNALIGNED_LENGTH);
            let mut contig = random_sequence(rng, length);

            for index in pieces {
                let chain = &mut chains[index];
                chain.offset = contig.len() as u64;

                match chain.strand {
                    Strand::Positive => contig.extend_from_slice(&chain.sequence),
                    Strand::Negative => contig.extend(reverse_complement(&chain.sequence)),
                }

                let length = rng.random_range(UNALIGNED_LENGTH);
                contig.extend(random_sequence(rng, length));
            }

            contig
        })
        .collect()
}

/// Writes a FASTA file and its `.fai` index.
fn write_fasta(path: &Path, contigs: &[Vec<u8>]) -> Result<()> {
    let mut fasta = Vec::new();
    let mut index = String::new();

    for (contig, bases) in contigs.iter().enumerate() {
        let name = contig_name(contig);
        writeln!(fasta, ">{name}")?;

        let offset = fasta.len();

        for line in bases.chunks(LINE_BASES) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }

        index.push_str(&format!(
            "{name}\t{}\t{offset}\t{LINE_BASES}\t{}\n",
            bases.len(),
            LINE_BASES + 1
        ));
    }

    std::fs::write(path, fasta).with_context(|| format!("writing `{}`", path.display()))?;

    let mut index_path = path.as_os_str().to_owned();
    index_path.push(".fai");
    let index_path = PathBuf::from(index_path);

    std::fs::write(&index_path, index)
        .with_context(|| format!("writing `{}`", index_path.display()))
}

/// Writes the chain file.
fn write_chains(
    path: &Path,
    chains: &[Chain],
    reference: &[Vec<u8>],
    query: &[Vec<u8>],
) -> Result<()> {
    let file = File::create(path).with_context(|| format!("creating `{}`", path.display()))?;
    let mut writer = BufWriter::new(file);

    for (index, chain) in chains.iter().enumerate() {
        let first = chain.blocks.first().expect("chains should have blocks");
        let last = chain.blocks.last().expect("chains should have blocks");

        let query_size = query[chain.query_contig].len() as u64;
        let length = chain.sequence.len() as u64;

        let (query_start, query_end, strand) = match chain.strand {
            Strand::Positive => (chain.offset, chain.offset + length, '+'),
            Strand::Negative => (
                query_size - (chain.offset + length),
                query_size - chain.offset,
                '-',
            ),
        };

        writeln!(
            writer,
            "chain {} {} {} + {} {} {} {query_size} {strand} {query_start} {query_end} {}",
            chain.score,
            contig_name(chain.reference_contig),
            reference[chain.reference_contig].len(),
            first.start,
            last.start + last.size,
            contig_name(chain.query_contig),
            index + 1
        )?;

        for block in &chain.blocks[..chain.blocks.len() - 1] {
            writeln!(writer, "{}\t{}\t{}", block.size, block.dt, block.dq)?;
        }

        writeln!(writer, "{}", last.size)?;
        writeln!(writer)?;
    }

    writer
        .flush()
        .with_context(|| format!("writing `{}`", path.display()))
}
