//! Random VCF records.

use std::io::Write;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use omics::coordinate::Strand;
use rand::Rng;
use rand::RngExt as _;
use rand::seq::IndexedRandom as _;
use rand::seq::SliceRandom as _;

use crate::SharedSituations;
use crate::VcfSituations;
use crate::fasta::IndexedFasta;
use crate::index::Candidates;
use crate::index::Index;
use crate::index::MultipleChains;
use crate::index::Region;
use crate::index::generate_records;
use crate::index::single_segment;
use crate::sequence::BASES;
use crate::sequence::forward_bounds;
use crate::sequence::is_nucleotide;
use crate::sequence::oriented;
use crate::sequence::random_base_except;
use crate::sequence::random_sequence;
use crate::situation::Situation;
use crate::situation::counts;

/// The maximum number of bases scanned at once when searching for bases
/// that differ between the reference and query genomes.
const MISMATCH_WINDOW: u64 = 1_000;

/// The genotypes of biallelic records.
const BIALLELIC_GENOTYPES: [&str; 2] = ["0/1", "1/1"];

/// The genotype of multiallelic records.
const MULTIALLELIC_GENOTYPE: &str = "1/2";

/// The largest supported maximum indel length.
const MAX_INDEL_LENGTH: u64 = 1_000_000;

/// Options for generating VCF records.
#[derive(Clone, Debug, clap::Args)]
pub struct Options {
    /// The number of records.
    #[arg(short = 'n', long, default_value_t = 10_000)]
    pub count: usize,

    /// The maximum number of inserted or deleted bases in an indel.
    #[arg(long, default_value_t = 10)]
    pub max_indel_length: u64,

    /// The percentages of records for each situation.
    #[command(flatten)]
    pub situations: SharedSituations,

    /// The percentages of records for each VCF-specific situation.
    #[command(flatten)]
    pub vcf_situations: VcfSituations,
}

/// A generated VCF record.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Record {
    /// The reference contig.
    contig: String,

    /// The 0-based start of the reference allele.
    start: u64,

    /// The ID.
    name: String,

    /// The reference allele.
    reference: Vec<u8>,

    /// The alternate alleles.
    alternates: Vec<Vec<u8>>,

    /// The genotype of the sample.
    genotype: &'static str,
}

/// The bases where an interval maps to a single aligned block.
struct Mapped {
    /// The reference bases.
    reference: Vec<u8>,

    /// The query bases, oriented to the reference genome.
    query: Vec<u8>,

    /// The query strand.
    strand: Strand,

    /// The 0-based start of the mapped interval on the forward strand of the
    /// query contig.
    query_start: u64,
}

impl Mapped {
    /// Returns whether every reference base is a nucleotide that matches the
    /// query genome.
    fn is_matched(&self) -> bool {
        self.reference == self.query && self.reference.iter().all(|&base| is_nucleotide(base))
    }
}

/// The state shared while generating records.
struct Generator<'a> {
    /// The chain file index.
    index: &'a Index,

    /// The reference genome.
    reference: IndexedFasta,

    /// The query genome.
    query: IndexedFasta,

    /// The options.
    options: &'a Options,
}

impl Generator<'_> {
    /// Fetches reference bases.
    fn reference_bases(&mut self, contig: usize, start: u64, end: u64) -> Result<Vec<u8>> {
        self.reference
            .fetch(&self.index.contigs[contig].name, start, end)
    }

    /// Fetches the reference and oriented query bases of `[start, end)` if the
    /// interval maps to exactly one chain with exactly one segment that covers
    /// the whole interval.
    fn mapped(&mut self, contig: usize, start: u64, end: u64) -> Result<Option<Mapped>> {
        let results = self.index.liftover(contig, start, end);

        let Some(segment) = single_segment(&results, start, end) else {
            return Ok(None);
        };

        let strand = segment.query().strand();
        let (query_start, query_end) = forward_bounds(segment.query());
        let query = self
            .query
            .fetch(segment.query().contig().as_str(), query_start, query_end)?;

        Ok(Some(Mapped {
            reference: self.reference_bases(contig, start, end)?,
            query: oriented(query, strand),
            strand,
            query_start,
        }))
    }

    /// Creates a record.
    fn record(
        &self,
        contig: usize,
        start: u64,
        name: String,
        reference: Vec<u8>,
        alternates: Vec<Vec<u8>>,
        genotype: &'static str,
    ) -> Record {
        Record {
            contig: self.index.contigs[contig].name.clone(),
            start,
            name,
            reference,
            alternates,
            genotype,
        }
    }

    /// Generates the records for a situation.
    fn generate(
        &mut self,
        situation: Situation,
        count: usize,
        rng: &mut impl Rng,
    ) -> Result<Vec<Record>> {
        match situation {
            Situation::Contained | Situation::OppositeStrand => {
                let strand = if situation == Situation::Contained {
                    Strand::Positive
                } else {
                    Strand::Negative
                };

                let candidates =
                    Candidates::new(self.index.clean_pieces(Some(strand)), Region::len);

                generate_records(situation, count, &candidates, rng, |rng, region, name| {
                    let start = rng.random_range(region.start..region.end);

                    let Some(mapped) = self.mapped(region.contig, start, start + 1)? else {
                        return Ok(None);
                    };

                    if mapped.strand != strand || !mapped.is_matched() {
                        return Ok(None);
                    }

                    let alternate = random_base_except(rng, &mapped.reference);
                    let genotype = genotype(rng);

                    Ok(Some(self.record(
                        region.contig,
                        start,
                        name,
                        mapped.reference,
                        vec![vec![alternate]],
                        genotype,
                    )))
                })
            }
            Situation::Indel => {
                let max = self.options.max_indel_length;
                let candidates = Candidates::new(self.index.clean_pieces(None), Region::len);

                generate_records(situation, count, &candidates, rng, |rng, region, name| {
                    let deletion = rng.random_bool(0.5);

                    let length = if !deletion {
                        1
                    } else if region.len() >= 2 {
                        rng.random_range(2..=(max + 1).min(region.len()))
                    } else {
                        return Ok(None);
                    };

                    let start = rng.random_range(region.start..=region.end - length);

                    let Some(mapped) = self.mapped(region.contig, start, start + length)? else {
                        return Ok(None);
                    };

                    // NOTE: lifting an indel over to the negative strand moves
                    // its padding base to the base before the mapped interval,
                    // so that base must exist.
                    if !mapped.is_matched()
                        || (mapped.strand == Strand::Negative && mapped.query_start == 0)
                    {
                        return Ok(None);
                    }

                    let reference = mapped.reference;

                    let alternate = if deletion {
                        reference[..1].to_vec()
                    } else {
                        let inserted = rng.random_range(1..=max);
                        let mut alternate = reference.clone();
                        alternate.extend(random_sequence(rng, inserted));
                        alternate
                    };

                    let genotype = genotype(rng);

                    Ok(Some(self.record(
                        region.contig,
                        start,
                        name,
                        reference,
                        vec![alternate],
                        genotype,
                    )))
                })
            }
            Situation::SwappedAlleles => {
                let candidates = Candidates::new(self.index.clean_pieces(None), Region::len);

                generate_records(situation, count, &candidates, rng, |rng, region, name| {
                    let Some(start) = self.sample_mismatch(rng, region)? else {
                        return Ok(None);
                    };

                    let Some(mapped) = self.mapped(region.contig, start, start + 1)? else {
                        return Ok(None);
                    };

                    let (reference, query) = (mapped.reference[0], mapped.query[0]);

                    if reference == query || !is_nucleotide(reference) || !is_nucleotide(query) {
                        return Ok(None);
                    }

                    let genotype = genotype(rng);

                    Ok(Some(self.record(
                        region.contig,
                        start,
                        name,
                        vec![reference],
                        vec![vec![query]],
                        genotype,
                    )))
                })
            }
            Situation::ReferenceMismatch | Situation::Multiallelic => {
                let candidates = Candidates::new(self.index.clean_pieces(None), Region::len);

                generate_records(situation, count, &candidates, rng, |rng, region, name| {
                    let start = rng.random_range(region.start..region.end);

                    let Some(mapped) = self.mapped(region.contig, start, start + 1)? else {
                        return Ok(None);
                    };

                    let (reference, query) = (mapped.reference[0], mapped.query[0]);

                    if !is_nucleotide(reference) || !is_nucleotide(query) {
                        return Ok(None);
                    }

                    let (reference, alternates, genotype) =
                        if situation == Situation::ReferenceMismatch {
                            let mismatched = random_base_except(rng, &[reference, query]);
                            let alternate = random_base_except(rng, &[mismatched, query]);
                            (mismatched, vec![vec![alternate]], genotype(rng))
                        } else if reference == query {
                            let mut alternates = BASES
                                .into_iter()
                                .filter(|&base| base != reference)
                                .map(|base| vec![base])
                                .collect::<Vec<_>>();
                            alternates.shuffle(rng);
                            alternates.truncate(2);
                            (reference, alternates, MULTIALLELIC_GENOTYPE)
                        } else {
                            return Ok(None);
                        };

                    Ok(Some(self.record(
                        region.contig,
                        start,
                        name,
                        vec![reference],
                        alternates,
                        genotype,
                    )))
                })
            }
            Situation::SpansGap => {
                let max = self.options.max_indel_length + 1;
                let candidates = Candidates::new(self.index.gaps(), |gap| gap.weight(2, max));

                generate_records(situation, count, &candidates, rng, |rng, gap, name| {
                    let Some((start, end)) = gap.sample(rng, 2, max) else {
                        return Ok(None);
                    };

                    let results = self.index.liftover(gap.contig, start, end);

                    if !matches!(results.as_deref(), Some([result]) if result.segments().len() >= 2)
                    {
                        return Ok(None);
                    }

                    self.deletion(rng, gap.contig, start, end, name)
                })
            }
            Situation::Unaligned => {
                let candidates = Candidates::new(self.index.unaligned(), Region::len);

                generate_records(situation, count, &candidates, rng, |rng, region, name| {
                    let start = rng.random_range(region.start..region.end);

                    if self
                        .index
                        .liftover(region.contig, start, start + 1)
                        .is_some()
                    {
                        return Ok(None);
                    }

                    self.snv(rng, region.contig, start, name)
                })
            }
            Situation::MultipleChains => {
                let max = self.options.max_indel_length + 1;

                let candidates =
                    Candidates::new(self.index.multiple_chains(), |candidate| match candidate {
                        MultipleChains::Overlap(region) => region.len(),
                        MultipleChains::Junction(junction) => junction.weight(2, max),
                    });

                generate_records(
                    situation,
                    count,
                    &candidates,
                    rng,
                    |rng, candidate, name| {
                        // NOTE: an SNV cannot touch two chains at a junction,
                        // so junctions get a padded deletion that covers at
                        // least one base of each block instead.
                        let (contig, start, end) = match candidate {
                            MultipleChains::Overlap(region) => {
                                let start = rng.random_range(region.start..region.end);
                                (region.contig, start, start + 1)
                            }
                            MultipleChains::Junction(junction) => {
                                let Some((start, end)) = junction.sample(rng, 2, max) else {
                                    return Ok(None);
                                };

                                (junction.contig, start, end)
                            }
                        };

                        let accepted = self
                            .index
                            .liftover(contig, start, end)
                            .is_some_and(|results| results.len() >= 2);

                        if !accepted {
                            return Ok(None);
                        }

                        if end - start == 1 {
                            self.snv(rng, contig, start, name)
                        } else {
                            self.deletion(rng, contig, start, end, name)
                        }
                    },
                )
            }
        }
    }

    /// Creates a biallelic SNV at `start` with the reference base as the
    /// reference allele and a random other base as the alternate allele.
    ///
    /// Returns [`None`] if the reference base is not a nucleotide.
    fn snv(
        &mut self,
        rng: &mut impl Rng,
        contig: usize,
        start: u64,
        name: String,
    ) -> Result<Option<Record>> {
        let reference = self.reference_bases(contig, start, start + 1)?;

        if !is_nucleotide(reference[0]) {
            return Ok(None);
        }

        let alternate = random_base_except(rng, &reference);
        let genotype = genotype(rng);

        Ok(Some(self.record(
            contig,
            start,
            name,
            reference,
            vec![vec![alternate]],
            genotype,
        )))
    }

    /// Creates a padded deletion of the reference bases `[start, end)` (all
    /// but the first base are deleted).
    ///
    /// Returns [`None`] if any reference base is not a nucleotide.
    fn deletion(
        &mut self,
        rng: &mut impl Rng,
        contig: usize,
        start: u64,
        end: u64,
        name: String,
    ) -> Result<Option<Record>> {
        let reference = self.reference_bases(contig, start, end)?;

        if !reference.iter().all(|&base| is_nucleotide(base)) {
            return Ok(None);
        }

        let alternate = reference[..1].to_vec();
        let genotype = genotype(rng);

        Ok(Some(self.record(
            contig,
            start,
            name,
            reference,
            vec![alternate],
            genotype,
        )))
    }

    /// Samples a position within a clean piece of a block where the reference
    /// and query bases are different nucleotides.
    ///
    /// A window of the piece is scanned, and [`None`] is returned if the window
    /// has no such positions.
    fn sample_mismatch(&mut self, rng: &mut impl Rng, region: &Region) -> Result<Option<u64>> {
        let block = self.index.block(
            region.contig,
            region.block.expect("clean pieces have blocks"),
        );

        let length = region.len().min(MISMATCH_WINDOW);
        let start = rng.random_range(region.start..=region.end - length);
        let end = start + length;

        let reference = self.reference_bases(region.contig, start, end)?;
        let (query_start, query_end) = block.query_bounds(start, end);
        let query = self
            .query
            .fetch(&block.query_contig, query_start, query_end)?;
        let query = oriented(query, block.strand);

        let mismatches = reference
            .iter()
            .zip(&query)
            .enumerate()
            .filter(|&(_, (&r, &q))| r != q && is_nucleotide(r) && is_nucleotide(q))
            .map(|(offset, _)| start + offset as u64)
            .collect::<Vec<_>>();

        Ok(mismatches.choose(rng).copied())
    }
}

/// Chooses a random biallelic genotype.
fn genotype(rng: &mut impl Rng) -> &'static str {
    BIALLELIC_GENOTYPES[rng.random_range(0..BIALLELIC_GENOTYPES.len())]
}

/// Generates random VCF records for the reference genome of `chain`.
///
/// Each record's ID is the situation it exercises followed by its index within
/// that situation (e.g., `spans_gap_12`). Every record is verified to have the
/// intended liftover outcome against the chain file and both genomes. The
/// FASTA files must be indexed. The records are sorted by contig name and
/// position.
pub fn generate(
    chain: &Path,
    reference_fasta: &Path,
    query_fasta: &Path,
    options: &Options,
    rng: &mut impl Rng,
    writer: &mut impl Write,
) -> Result<()> {
    if !(1..=MAX_INDEL_LENGTH).contains(&options.max_indel_length) {
        bail!("the maximum indel length must be between 1 and {MAX_INDEL_LENGTH}");
    }

    let percentages = options
        .situations
        .percentages()
        .into_iter()
        .chain(options.vcf_situations.percentages())
        .collect::<Vec<_>>();

    let counts = counts(options.count, &percentages)?;
    let index = Index::load(chain)?;

    let mut generator = Generator {
        index: &index,
        reference: IndexedFasta::open(reference_fasta)?,
        query: IndexedFasta::open(query_fasta)?,
        options,
    };

    let mut records = Vec::with_capacity(options.count);

    for (situation, count) in counts {
        let generated = generator
            .generate(situation, count, rng)
            .with_context(|| format!("generating {situation} VCF records"))?;
        records.extend(generated);
    }

    records.sort();

    writeln!(writer, "##fileformat=VCFv4.3")?;

    for contig in &index.contigs {
        writeln!(
            writer,
            "##contig=<ID={},length={}>",
            contig.name, contig.size
        )?;
    }

    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        writer,
        "{}",
        [
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"
        ]
        .join("\t")
    )?;

    for record in records {
        let alternates = record
            .alternates
            .iter()
            .map(|allele| String::from_utf8_lossy(allele))
            .collect::<Vec<_>>()
            .join(",");

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{alternates}\t.\tPASS\t.\tGT\t{}",
            record.contig,
            record.start + 1,
            record.name,
            String::from_utf8_lossy(&record.reference),
            record.genotype
        )?;
    }

    writer.flush()?;

    Ok(())
}
