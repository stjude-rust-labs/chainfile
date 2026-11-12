//! The `liftover bed` subcommand.

use std::collections::BTreeMap;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use chainfile::liftover::Machine;
use noodles::bed;
use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;
use tracing::debug;
use tracing::info;
use tracing::trace;
use tracing::warn;

use crate::io::load_machines;
use crate::io::open_input;
use crate::io::open_output;
use crate::parallel;
use crate::parallel::Chunk;
use crate::sequence::forward_bounds;

/// The index of the strand column (column 6) within the other fields (the
/// columns after column 3).
const STRAND_FIELD: usize = 2;

/// The maximum number of columns that are lifted over (BED6).
const MAX_LIFTED_COLUMNS: usize = 6;

/// Lift over BED records from the reference genome to the query genome.
///
/// A record that maps to multiple segments (for example, when it spans a gap
/// in the alignment or aligns to more than one chain) is written once for
/// each segment. Columns 1–3 are lifted over, and the strand in column 6 is
/// flipped when a segment maps to the opposite strand. Columns 4 and 5 are
/// copied unchanged.
///
/// Records that do not map are written to the unmapped file (if provided),
/// each preceded by a `#` line explaining why.
///
/// The output is not sorted.
#[derive(clap::Args, Debug)]
pub struct Args {
    /// The chain file (optionally gzipped).
    #[arg(short, long)]
    chain: PathBuf,

    /// The BED file to lift over (optionally gzipped); `-` reads from standard
    /// input.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// The output file; standard output is used when not provided.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// The file to write unmapped records to.
    #[arg(short, long)]
    unmapped: Option<PathBuf>,

    /// Copy columns 7 and beyond unchanged.
    ///
    /// Without this flag or `--ignore-bed-plus`, records with more than six
    /// columns are an error, as those columns may contain coordinates (for
    /// example, the thick and block columns of BED12) that would be incorrect
    /// after liftover. Only use this flag when those columns do not depend on
    /// the record's coordinates.
    #[arg(long, conflicts_with = "ignore_bed_plus")]
    include_bed_plus: bool,

    /// Drop columns 7 and beyond from lifted records.
    ///
    /// Unmapped records are still written with all of their columns.
    #[arg(long)]
    ignore_bed_plus: bool,

    /// The number of worker threads; defaults to the number of available
    /// cores.
    #[arg(short, long)]
    threads: Option<NonZeroUsize>,
}

/// The reason a record could not be lifted over.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum Unmapped {
    /// The record has no length.
    ZeroLength,

    /// The record's contig is not in the chain file.
    UnknownContig,

    /// No part of the record aligns to the query genome.
    NoAlignment,
}

impl Unmapped {
    /// Gets a description of the reason.
    fn description(&self) -> &'static str {
        match self {
            Unmapped::ZeroLength => "zero-length interval",
            Unmapped::UnknownContig => "contig not in chain file",
            Unmapped::NoAlignment => "no alignment in query genome",
        }
    }
}

/// A lifted segment of a record.
#[derive(Debug, Eq, PartialEq)]
struct Segment {
    /// The start on the positive strand of the query contig.
    start: Coordinate,

    /// The end on the positive strand of the query contig.
    end: Coordinate,

    /// Whether the segment maps to the opposite strand.
    flipped: bool,
}

/// Statistics about lifted over records.
#[derive(Debug, Default)]
struct Stats {
    /// The number of records read.
    records: u64,

    /// The number of records that were lifted over.
    mapped: u64,

    /// The number of records that could not be lifted over.
    unmapped: u64,

    /// The number of segments written.
    segments: u64,

    /// The number of segments that mapped to the opposite strand.
    flipped: u64,

    /// The number of records that were split into multiple segments.
    split: u64,

    /// The number of records with columns 7 and beyond.
    bed_plus: u64,

    /// The number and column count of the first record with columns 7 and
    /// beyond.
    first_bed_plus: Option<(u64, usize)>,

    /// The number of records that could not be lifted over for each reason.
    reasons: BTreeMap<Unmapped, u64>,
}

/// The output of lifting over a chunk of records.
#[derive(Debug, Default)]
struct Output {
    /// The lifted records.
    lifted: Vec<u8>,

    /// The unmapped records.
    unmapped: Vec<u8>,

    /// Statistics about the records.
    stats: Stats,
}

/// The state of a worker thread.
struct Worker<'a> {
    /// The liftover machine.
    machine: &'a Machine,

    /// The record buffer.
    record: bed::Record<3>,

    /// The line buffer for lifted records.
    line: Vec<u8>,

    /// The segment buffer for lifted records.
    segments: Vec<Segment>,

    /// The previous record's contig and whether it is in the chain file.
    ///
    /// Records are typically grouped by contig, so this is reused when the
    /// name matches.
    previous: Option<(Contig, bool)>,
}

/// Runs the `liftover bed` subcommand.
pub fn liftover(args: Args) -> Result<()> {
    let threads = parallel::threads(args.threads);

    info!(
        chain = %args.chain.display(),
        input = %args.input.display(),
        output = %args.output.as_deref().map(|p| p.display().to_string()).unwrap_or_else(|| String::from("-")),
        unmapped = %args
            .unmapped
            .as_deref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| String::from("none")),
        include_bed_plus = args.include_bed_plus,
        ignore_bed_plus = args.ignore_bed_plus,
        threads,
        "lifting over BED records"
    );

    let machines = load_machines(&args.chain, threads)?;

    let input = open_input(&args.input)?;
    let mut writer = open_output(args.output.as_deref())?;
    let mut unmapped_writer = args
        .unmapped
        .as_deref()
        .map(|path| open_output(Some(path)))
        .transpose()?;

    if unmapped_writer.is_none() {
        debug!("no unmapped file was provided, so unmapped records will be discarded");
    }

    let workers = machines
        .iter()
        .map(|machine| Worker {
            machine,
            record: bed::Record::default(),
            line: Vec::new(),
            segments: Vec::new(),
            previous: None,
        })
        .collect();

    let mut stats = Stats::default();

    parallel::process(
        input,
        "BED",
        Some(b'#'),
        workers,
        |worker, chunk, output| lift_chunk(worker, chunk, output, &args),
        |output: Output| {
            writer
                .write_all(&output.lifted)
                .context("writing lifted BED records")?;

            if let Some(writer) = unmapped_writer.as_mut() {
                writer
                    .write_all(&output.unmapped)
                    .context("writing unmapped BED records")?;
            }

            stats.merge(output.stats, args.ignore_bed_plus);
            Ok(())
        },
    )?;

    writer.flush()?;

    if let Some(mut writer) = unmapped_writer {
        writer.flush()?;
    }

    let Stats {
        records: n,
        mapped,
        unmapped,
        segments: segments_written,
        flipped,
        split,
        bed_plus,
        reasons,
        ..
    } = stats;

    info!(
        "lifted over {mapped} of {n} BED records into {segments_written} segments; {unmapped} \
         records were unmapped"
    );

    if split > 0 {
        info!("{split} records were split into multiple segments");
    }

    if flipped > 0 {
        info!("{flipped} segments mapped to the opposite strand");
    }

    if bed_plus > 0 {
        let action = if args.ignore_bed_plus {
            "dropped"
        } else {
            "copied unchanged"
        };
        info!("{bed_plus} records had columns 7 and beyond {action}");
    }

    for (reason, count) in reasons {
        info!("  {}: {count}", reason.description());
    }

    Ok(())
}

impl Stats {
    /// Merges the statistics of the next chunk of records.
    ///
    /// The first record with columns 7 and beyond is logged.
    fn merge(&mut self, other: Stats, ignore_bed_plus: bool) {
        if let (None, Some((record, columns))) = (self.first_bed_plus, other.first_bed_plus) {
            if ignore_bed_plus {
                info!(
                    record,
                    columns, "dropping columns 7 and beyond from lifted records"
                );
            } else {
                warn!(
                    record,
                    columns,
                    "copying columns 7 and beyond unchanged; ensure they do not depend on the \
                     record's coordinates"
                );
            }

            self.first_bed_plus = other.first_bed_plus;
        }

        self.records += other.records;
        self.mapped += other.mapped;
        self.unmapped += other.unmapped;
        self.segments += other.segments;
        self.flipped += other.flipped;
        self.split += other.split;
        self.bed_plus += other.bed_plus;

        for (reason, count) in other.reasons {
            *self.reasons.entry(reason).or_default() += count;
        }
    }
}

/// Lifts over a chunk of records.
fn lift_chunk(
    worker: &mut Worker<'_>,
    chunk: &Chunk,
    output: &mut Output,
    args: &Args,
) -> Result<()> {
    let Worker {
        machine,
        record,
        line,
        segments,
        previous,
    } = worker;

    let stats = &mut output.stats;
    let mut reader = bed::io::Reader::<3, _>::new(&chunk.data[..]);
    let mut unmapped_writer = args
        .unmapped
        .is_some()
        .then(|| bed::io::Writer::<3, _>::new(&mut output.unmapped));

    for n in chunk.first..chunk.first + chunk.lines as u64 {
        reader
            .read_record(record)
            .with_context(|| format!("reading BED record {n}"))?;

        stats.records += 1;

        let columns = 3 + record.other_fields().len();
        if columns > MAX_LIFTED_COLUMNS && !args.include_bed_plus && !args.ignore_bed_plus {
            bail!(
                "BED record {n} has {columns} columns, but only up to {MAX_LIFTED_COLUMNS} can be \
                 lifted over; columns 7 and beyond may contain coordinates (such as the thick and \
                 block columns of BED12) that would be incorrect after liftover. Rerun with \
                 `--ignore-bed-plus` to drop them or, if they do not depend on the record's \
                 coordinates, with `--include-bed-plus` to copy them unchanged"
            );
        }

        if columns > MAX_LIFTED_COLUMNS {
            if stats.bed_plus == 0 {
                stats.first_bed_plus = Some((n, columns));
            }

            stats.bed_plus += 1;
        }

        let name: &[u8] = record.reference_sequence_name();
        let (contig, known) = match previous.take() {
            Some(previous) if previous.0.as_str().as_bytes() == name => previous,
            _ => {
                let contig = Contig::try_new(record.reference_sequence_name().to_string())
                    .with_context(|| format!("parsing the contig of BED record {n}"))?;
                let known = machine
                    .reference_chromosomes()
                    .contains_key(contig.as_str());
                (contig, known)
            }
        };
        *previous = Some((contig.clone(), known));

        // NOTE: `noodles` converts BED's 0-based start into a 1-based
        // `Position` (and leaves the exclusive end as is), so one is subtracted
        // from the start to get back to a 0-based, half-open interbase
        // interval.
        let start = record
            .feature_start()
            .with_context(|| format!("parsing the start of BED record {n}"))?
            .get() as u64
            - 1;
        let end = record
            .feature_end()
            .transpose()
            .with_context(|| format!("parsing the end of BED record {n}"))?
            .map(|position| position.get() as u64)
            .unwrap_or_default();

        if start > end {
            bail!("BED record {n} has a start ({start}) greater than its end ({end})");
        }

        trace!(record = n, %contig, start, end, columns, "read BED record");

        let interval = Interval::try_new(
            Coordinate::new(contig.clone(), Strand::Positive, start),
            Coordinate::new(contig.clone(), Strand::Positive, end),
        )
        .with_context(|| format!("creating interval for BED record {n}"))?;

        match lift(machine, interval, known, segments) {
            Ok(()) => {
                stats.mapped += 1;

                if segments.len() > 1 {
                    stats.split += 1;
                    debug!(
                        record = n,
                        %contig,
                        start,
                        end,
                        segments = segments.len(),
                        "BED record lifted over into multiple segments"
                    );
                }

                for segment in segments.iter() {
                    trace!(
                        record = n,
                        contig = %segment.start.contig(),
                        start = segment.start.position().get(),
                        end = segment.end.position().get(),
                        flipped = segment.flipped,
                        "writing lifted segment"
                    );

                    if segment.flipped {
                        stats.flipped += 1;
                    }

                    write_segment(
                        &mut output.lifted,
                        line,
                        record,
                        segment,
                        args.ignore_bed_plus,
                    )
                    .with_context(|| format!("writing BED record {n}"))?;
                    stats.segments += 1;
                }
            }
            Err(reason) => {
                stats.unmapped += 1;
                *stats.reasons.entry(reason).or_default() += 1;

                debug!(
                    record = n,
                    %contig,
                    start,
                    end,
                    reason = reason.description(),
                    "BED record could not be lifted over"
                );

                if let Some(writer) = unmapped_writer.as_mut() {
                    writeln!(writer.get_mut(), "#{}", reason.description())?;
                    writer
                        .write_record(record)
                        .with_context(|| format!("writing unmapped BED record {n}"))?;
                }
            }
        }
    }

    Ok(())
}

/// Lifts over an interval on the positive strand of the reference genome.
///
/// `known` is whether the interval's contig is in the chain file. The lifted
/// segments are written to `segments`, which is reused between records.
fn lift(
    machine: &Machine,
    interval: Interval,
    known: bool,
    segments: &mut Vec<Segment>,
) -> std::result::Result<(), Unmapped> {
    segments.clear();

    if interval.start().position() == interval.end().position() {
        return Err(Unmapped::ZeroLength);
    }

    if !known {
        return Err(Unmapped::UnknownContig);
    }

    let results = machine
        .liftover(interval.clone())
        .ok_or(Unmapped::NoAlignment)?;

    if results.len() > 1 {
        debug!(
            %interval,
            chains = results.len(),
            "interval aligns to multiple chains"
        );
    }

    segments.extend(
        results
            .into_iter()
            .flat_map(|result| result.into_segments())
            .map(|pair| {
                let query = pair.query();
                let (start, end) = forward_bounds(query);

                Segment {
                    start: Coordinate::new(query.contig().clone(), Strand::Positive, start),
                    end: Coordinate::new(query.contig().clone(), Strand::Positive, end),
                    flipped: pair.query().strand() != pair.reference().strand(),
                }
            }),
    );

    Ok(())
}

/// Writes a lifted segment of a record.
///
/// The line is formatted into `line` (which is reused between records to
/// avoid allocations) and then written in one call.
fn write_segment<W: Write + ?Sized>(
    writer: &mut W,
    line: &mut Vec<u8>,
    record: &bed::Record<3>,
    segment: &Segment,
    ignore_bed_plus: bool,
) -> Result<()> {
    let limit = if ignore_bed_plus {
        // NOTE: the other fields exclude the three required columns (contig,
        // start, and end), so only the first three of them (columns 4–6) are
        // kept to stay within BED6.
        MAX_LIFTED_COLUMNS - 3
    } else {
        usize::MAX
    };

    line.clear();
    line.extend_from_slice(segment.start.contig().as_str().as_bytes());
    line.push(b'\t');
    push_number(line, segment.start.position().get());
    line.push(b'\t');
    push_number(line, segment.end.position().get());

    for (i, value) in record.other_fields().iter().take(limit).enumerate() {
        let value: &[u8] = value.as_ref();
        line.push(b'\t');

        if i == STRAND_FIELD && segment.flipped {
            // NOTE: strands other than `+` and `-` (such as `.`) are unknown,
            // so they are copied unchanged.
            let strand = match value {
                b"+" => Some(Strand::Positive),
                b"-" => Some(Strand::Negative),
                _ => None,
            };

            if let Some(strand) = strand {
                line.push(match strand.complement() {
                    Strand::Positive => b'+',
                    Strand::Negative => b'-',
                });
                continue;
            }
        }

        line.extend_from_slice(value);
    }

    line.push(b'\n');
    writer.write_all(line)?;
    Ok(())
}

/// Appends the decimal representation of a number to a buffer.
fn push_number(buffer: &mut Vec<u8>, mut value: u64) {
    let mut digits = [0u8; 20];
    let mut i = digits.len();

    loop {
        i -= 1;
        digits[i] = b'0' + (value % 10) as u8;
        value /= 10;

        if value == 0 {
            break;
        }
    }

    buffer.extend_from_slice(&digits[i..]);
}
