//! The `liftover bed` subcommand.

use std::collections::BTreeMap;
use std::io::Write;
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

use crate::io::load_machine;
use crate::io::open_input;
use crate::io::open_output;
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

/// Runs the `liftover bed` subcommand.
pub fn liftover(args: Args) -> Result<()> {
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
        "lifting over BED records"
    );

    let machine = load_machine(&args.chain)?;

    let mut reader = bed::io::Reader::<3, _>::new(open_input(&args.input)?);
    let mut writer = open_output(args.output.as_deref())?;
    let mut unmapped_writer = args
        .unmapped
        .as_deref()
        .map(|path| open_output(Some(path)).map(bed::io::Writer::<3, _>::new))
        .transpose()?;

    let mut record = bed::Record::<3>::default();
    let mut line = Vec::new();
    let mut segments = Vec::new();

    // NOTE: records are typically grouped by contig, so the previous record's
    // contig (and whether it is in the chain file) is reused when the name
    // matches.
    let mut previous: Option<(Contig, bool)> = None;
    let mut n = 0u64;
    let mut mapped = 0u64;
    let mut unmapped = 0u64;
    let mut segments_written = 0u64;
    let mut flipped = 0u64;
    let mut split = 0u64;
    let mut bed_plus = 0u64;
    let mut reasons = BTreeMap::<Unmapped, u64>::new();

    if unmapped_writer.is_none() {
        debug!("no unmapped file was provided, so unmapped records will be discarded");
    }

    loop {
        let read = reader
            .read_record(&mut record)
            .with_context(|| format!("reading BED record {}", n + 1))?;

        if read == 0 {
            break;
        }

        n += 1;

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
            if bed_plus == 0 {
                if args.ignore_bed_plus {
                    info!(
                        record = n,
                        columns, "dropping columns 7 and beyond from lifted records"
                    );
                } else {
                    warn!(
                        record = n,
                        columns,
                        "copying columns 7 and beyond unchanged; ensure they do not depend on the \
                         record's coordinates"
                    );
                }
            }

            bed_plus += 1;
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
        previous = Some((contig.clone(), known));

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

        match lift(&machine, interval, known, &mut segments) {
            Ok(()) => {
                mapped += 1;

                if segments.len() > 1 {
                    split += 1;
                    debug!(
                        record = n,
                        %contig,
                        start,
                        end,
                        segments = segments.len(),
                        "BED record lifted over into multiple segments"
                    );
                }

                for segment in &segments {
                    trace!(
                        record = n,
                        contig = %segment.start.contig(),
                        start = segment.start.position().get(),
                        end = segment.end.position().get(),
                        flipped = segment.flipped,
                        "writing lifted segment"
                    );

                    if segment.flipped {
                        flipped += 1;
                    }

                    write_segment(
                        &mut writer,
                        &mut line,
                        &record,
                        segment,
                        args.ignore_bed_plus,
                    )
                    .with_context(|| format!("writing BED record {n}"))?;
                    segments_written += 1;
                }
            }
            Err(reason) => {
                unmapped += 1;
                *reasons.entry(reason).or_default() += 1;

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
                        .write_record(&record)
                        .with_context(|| format!("writing unmapped BED record {n}"))?;
                }
            }
        }
    }

    writer.flush()?;

    if let Some(mut writer) = unmapped_writer {
        writer.get_mut().flush()?;
    }

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
