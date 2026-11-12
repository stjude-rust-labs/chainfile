//! Random BED records.

use std::io::Write;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use omics::coordinate::Strand;
use rand::Rng;
use rand::RngExt as _;

use crate::SharedSituations;
use crate::index::Candidates;
use crate::index::Index;
use crate::index::MultipleChains;
use crate::index::Region;
use crate::index::generate_records;
use crate::index::single_segment;
use crate::situation::Situation;
use crate::situation::counts;

/// Options for generating BED records.
#[derive(Clone, Debug, clap::Args)]
pub struct Options {
    /// The number of records.
    #[arg(short = 'n', long, default_value_t = 10_000)]
    pub count: usize,

    /// The minimum length of a record.
    #[arg(long, default_value_t = 1)]
    pub min_length: u64,

    /// The maximum length of a record.
    #[arg(long, default_value_t = 1_000)]
    pub max_length: u64,

    /// The percentages of records for each situation.
    #[command(flatten)]
    pub situations: SharedSituations,
}

/// A generated BED record.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Record {
    /// The reference contig.
    contig: String,

    /// The 0-based start.
    start: u64,

    /// The 0-based, exclusive end.
    end: u64,

    /// The name.
    name: String,
}

/// Generates random BED6 records for the reference genome of `chain`.
///
/// Each record's name is the situation it exercises followed by its index
/// within that situation (e.g., `spans_gap_12`). Every record is verified to
/// have the intended liftover outcome. The records are sorted by contig name
/// and start.
pub fn generate(
    chain: &Path,
    options: &Options,
    rng: &mut impl Rng,
    writer: &mut impl Write,
) -> Result<()> {
    if options.min_length == 0 {
        bail!("the minimum length must be at least one");
    }

    if options.min_length > options.max_length {
        bail!(
            "the minimum length ({}) is greater than the maximum length ({})",
            options.min_length,
            options.max_length
        );
    }

    let counts = counts(options.count, &options.situations.percentages())?;
    let index = Index::load(chain)?;

    let mut records = Vec::with_capacity(options.count);

    for (situation, count) in counts {
        let generated = generate_situation(&index, options, situation, count, rng)
            .with_context(|| format!("generating {situation} BED records"))?;
        records.extend(generated);
    }

    records.sort();

    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t0\t+",
            record.contig, record.start, record.end, record.name
        )?;
    }

    writer.flush()?;

    Ok(())
}

/// Generates the records for a situation.
fn generate_situation(
    index: &Index,
    options: &Options,
    situation: Situation,
    count: usize,
    rng: &mut impl Rng,
) -> Result<Vec<Record>> {
    let usable = |region: &Region| {
        if region.len() >= options.min_length {
            region.len()
        } else {
            0
        }
    };

    let record = |contig: usize, start: u64, end: u64, name: String| Record {
        contig: index.contigs[contig].name.clone(),
        start,
        end,
        name,
    };

    match situation {
        Situation::Contained | Situation::OppositeStrand => {
            let strand = if situation == Situation::Contained {
                Strand::Positive
            } else {
                Strand::Negative
            };

            let candidates = Candidates::new(index.clean_pieces(Some(strand)), usable);

            generate_records(situation, count, &candidates, rng, |rng, region, name| {
                let (start, end) = sample_within(rng, region, options);
                let results = index.liftover(region.contig, start, end);

                let accepted = single_segment(&results, start, end)
                    .is_some_and(|segment| segment.query().strand() == strand);

                Ok(accepted.then(|| record(region.contig, start, end, name)))
            })
        }
        Situation::Unaligned => {
            let candidates = Candidates::new(index.unaligned(), usable);

            generate_records(situation, count, &candidates, rng, |rng, region, name| {
                let (start, end) = sample_within(rng, region, options);
                let accepted = index.liftover(region.contig, start, end).is_none();
                Ok(accepted.then(|| record(region.contig, start, end, name)))
            })
        }
        Situation::MultipleChains => {
            let candidates =
                Candidates::new(index.multiple_chains(), |candidate| match candidate {
                    MultipleChains::Overlap(region) => usable(region),
                    MultipleChains::Junction(junction) => {
                        junction.weight(options.min_length, options.max_length)
                    }
                });

            generate_records(
                situation,
                count,
                &candidates,
                rng,
                |rng, candidate, name| {
                    let (contig, (start, end)) = match candidate {
                        MultipleChains::Overlap(region) => {
                            (region.contig, sample_within(rng, region, options))
                        }
                        MultipleChains::Junction(junction) => {
                            let Some(interval) =
                                junction.sample(rng, options.min_length, options.max_length)
                            else {
                                return Ok(None);
                            };

                            (junction.contig, interval)
                        }
                    };

                    let accepted = index
                        .liftover(contig, start, end)
                        .is_some_and(|results| results.len() >= 2);

                    Ok(accepted.then(|| record(contig, start, end, name)))
                },
            )
        }
        Situation::SpansGap => {
            let candidates = Candidates::new(index.gaps(), |gap| {
                gap.weight(options.min_length, options.max_length)
            });

            generate_records(situation, count, &candidates, rng, |rng, gap, name| {
                let Some((start, end)) = gap.sample(rng, options.min_length, options.max_length)
                else {
                    return Ok(None);
                };

                let results = index.liftover(gap.contig, start, end);

                let accepted = matches!(
                    results.as_deref(),
                    Some([result]) if result.segments().len() >= 2
                );

                Ok(accepted.then(|| record(gap.contig, start, end, name)))
            })
        }
        _ => bail!("{situation} records cannot be generated in BED format"),
    }
}

/// Samples an interval within a region with a length between the minimum and
/// maximum lengths.
///
/// The region must be at least the minimum length.
fn sample_within(rng: &mut impl Rng, region: &Region, options: &Options) -> (u64, u64) {
    let length = rng.random_range(options.min_length..=options.max_length.min(region.len()));
    let start = rng.random_range(region.start..=region.end - length);
    (start, start + length)
}
