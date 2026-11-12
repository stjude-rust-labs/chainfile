//! The `liftover vcf` subcommand.

mod swap;

use std::collections::BTreeMap;
use std::io::Write as _;
use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::Context as _;
use anyhow::Result;
use chainfile::liftover::Machine;
use noodles::core::Position as NoodlesPosition;
use noodles::vcf;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::header::record::value::map::Contig as NoodlesContig;
use noodles::vcf::header::record::value::map::Info;
use noodles::vcf::header::record::value::map::info::Number;
use noodles::vcf::header::record::value::map::info::Type;
use noodles::vcf::variant::io::Write as _;
use noodles::vcf::variant::record_buf::AlternateBases;
use noodles::vcf::variant::record_buf::RecordBuf;
use noodles::vcf::variant::record_buf::info::field::Value;
use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;
use omics::coordinate::position::base::Position;
use tracing::debug;
use tracing::info;

use crate::fasta::IndexedFasta;
use crate::io::load_machines;
use crate::io::open_input;
use crate::io::open_output;
use crate::parallel;
use crate::parallel::Chunk;
use crate::sequence::forward_bounds;
use crate::sequence::reverse_complement;

/// The INFO field that records why a record was rejected.
const REJECT_REASON_KEY: &str = "LiftoverRejectReason";

/// The INFO field that flags records whose alleles were swapped.
const SWAPPED_ALLELES_KEY: &str = "LiftoverSwappedAlleles";

/// The INFO field for the end position of a record.
const END_KEY: &str = "END";

/// The allele that denotes an overlapping deletion.
const OVERLAPPING_DELETION: &str = "*";

/// Lift over VCF records from the reference genome to the query genome.
///
/// A record is only lifted over when its entire reference allele maps to a
/// single contiguous segment of a single chain. Alleles are reverse
/// complemented when the record maps to the opposite strand, and indels are
/// re-anchored to the preceding base in the query genome.
///
/// The reference allele is checked against the query genome FASTA, and records
/// with a differing reference allele are rejected. The exception is a
/// biallelic single-base substitution whose alternate allele matches the query
/// genome: its alleles are swapped, its allele-indexed fields (such as `GT`,
/// `AD`, `PL`, `AF`, and `AC`) are rewritten to match, fields that cannot be
/// rewritten are dropped, and the `LiftoverSwappedAlleles` INFO flag is set.
/// `QUAL` and `FILTER` are not changed.
///
/// The contig header lines of the output are replaced with the contigs in the
/// query genome FASTA index. The output is not sorted.
///
/// Rejected records are written unchanged to the rejected file (if provided)
/// with the reason in the `LiftoverRejectReason` INFO field.
#[derive(clap::Args, Debug)]
pub struct Args {
    /// The chain file (optionally gzipped).
    #[arg(short, long)]
    chain: PathBuf,

    /// The query genome FASTA (the genome lifted over to), which must have a
    /// `.fai` index (and a `.gzi` index when BGZF-compressed).
    #[arg(short = 'f', long)]
    query_fasta: PathBuf,

    /// The VCF file to lift over (optionally gzipped); `-` reads from standard
    /// input.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// The output file; standard output is used when not provided.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// The file to write rejected records to.
    #[arg(short, long)]
    rejected: Option<PathBuf>,

    /// The number of worker threads; defaults to the number of available
    /// cores.
    #[arg(short, long)]
    threads: Option<NonZeroUsize>,
}

/// The reason a record could not be lifted over.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum Rejection {
    /// The record has no position.
    MissingPosition,

    /// The record has a symbolic or breakend allele.
    SymbolicAllele,

    /// The record has an `END` INFO field.
    EndField,

    /// The record's contig is not in the chain file.
    UnknownContig,

    /// The record does not align to the query genome.
    NoAlignment,

    /// The record aligns to more than one chain.
    MultipleMappings,

    /// The record does not align as a single contiguous segment.
    PartialMapping,

    /// The mapped contig is not in the query genome FASTA.
    ContigNotInFasta,

    /// The reference allele differs from the query genome.
    ReferenceMismatch,

    /// The query genome's base matches an alternate allele, and the alleles
    /// could not be swapped (because the record is multiallelic or its
    /// allele-indexed fields are invalid).
    SwappedAlleles,

    /// An indel on the opposite strand maps to the start of a contig, so it
    /// cannot be re-anchored.
    MissingAnchor,
}

impl Rejection {
    /// Gets the code that describes the rejection.
    fn code(&self) -> &'static str {
        match self {
            Rejection::MissingPosition => "missing_position",
            Rejection::SymbolicAllele => "symbolic_allele",
            Rejection::EndField => "end_field",
            Rejection::UnknownContig => "unknown_contig",
            Rejection::NoAlignment => "no_alignment",
            Rejection::MultipleMappings => "multiple_mappings",
            Rejection::PartialMapping => "partial_mapping",
            Rejection::ContigNotInFasta => "contig_not_in_fasta",
            Rejection::ReferenceMismatch => "reference_mismatch",
            Rejection::SwappedAlleles => "swapped_alleles",
            Rejection::MissingAnchor => "missing_anchor",
        }
    }
}

/// A lifted variant.
#[derive(Debug, Eq, PartialEq)]
struct Lifted {
    /// The query contig.
    contig: Contig,

    /// The position on the forward strand of the query contig.
    position: Position,

    /// The reference allele.
    reference: String,

    /// The alternate alleles.
    alternates: Vec<String>,

    /// Whether the reference and alternate alleles were swapped.
    swapped: bool,
}

/// Transformed alleles.
#[derive(Debug, Eq, PartialEq)]
struct Alleles {
    /// The number of bases the record shifted towards the start of the contig.
    offset: u64,

    /// The reference allele.
    reference: String,

    /// The alternate alleles.
    alternates: Vec<String>,

    /// Whether the reference and alternate alleles were swapped.
    swapped: bool,
}

/// The output of lifting over a chunk of records.
#[derive(Debug, Default)]
struct Output {
    /// The lifted records.
    lifted: Vec<u8>,

    /// The rejected records.
    rejected: Vec<u8>,

    /// The number of records read.
    records: u64,

    /// The number of records that were lifted over.
    lifted_count: u64,

    /// The number of lifted records whose alleles were swapped.
    swapped_count: u64,

    /// The number of records rejected for each reason.
    rejections: BTreeMap<Rejection, u64>,
}

/// The state of a worker thread.
struct Worker<'a> {
    /// The liftover machine.
    machine: &'a Machine,

    /// The query genome FASTA.
    fasta: IndexedFasta,

    /// The record buffer.
    record: RecordBuf,
}

/// The headers of the input and outputs.
struct Headers {
    /// The input header.
    input: vcf::Header,

    /// The header of the lifted records.
    output: vcf::Header,

    /// The header of the rejected records.
    rejected: vcf::Header,
}

/// Runs the `liftover vcf` subcommand.
pub fn liftover(args: Args) -> Result<()> {
    let threads = parallel::threads(args.threads);
    debug!(threads, "lifting over VCF records");

    let machines = load_machines(&args.chain, threads)?;
    let workers = machines
        .iter()
        .map(|machine| {
            Ok(Worker {
                machine,
                fasta: IndexedFasta::open(&args.query_fasta)?,
                record: RecordBuf::default(),
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut reader = vcf::io::Reader::new(open_input(&args.input)?);
    let header = reader.read_header().context("reading VCF header")?;

    let mut output_header = header.clone();
    let contigs = output_header.contigs_mut();
    contigs.clear();
    for (name, length) in workers[0].fasta.contigs() {
        let mut contig = Map::<NoodlesContig>::new();
        *contig.length_mut() = Some(usize::try_from(*length)?);
        contigs.insert(name.clone(), contig);
    }

    output_header.infos_mut().insert(
        String::from(SWAPPED_ALLELES_KEY),
        Map::<Info>::new(
            Number::Count(0),
            Type::Flag,
            "The reference and alternate alleles were swapped during liftover",
        ),
    );

    let mut writer = vcf::io::Writer::new(open_output(args.output.as_deref())?);
    writer.write_header(&output_header)?;

    let mut rejected_header = header.clone();
    rejected_header.infos_mut().insert(
        String::from(REJECT_REASON_KEY),
        Map::<Info>::new(
            Number::Count(1),
            Type::String,
            "The reason the record could not be lifted over",
        ),
    );

    let mut rejected_writer = match args.rejected.as_deref() {
        Some(path) => {
            let mut writer = vcf::io::Writer::new(open_output(Some(path))?);
            writer.write_header(&rejected_header)?;
            Some(writer)
        }
        None => None,
    };

    let write_rejected = rejected_writer.is_some();
    let headers = Headers {
        input: header,
        output: output_header,
        rejected: rejected_header,
    };

    let mut n = 0u64;
    let mut lifted_count = 0u64;
    let mut swapped_count = 0u64;
    let mut rejections = BTreeMap::<Rejection, u64>::new();

    // NOTE: the inner reader is positioned at the first record after the
    // header is read.
    parallel::process(
        reader.into_inner(),
        "VCF",
        None,
        workers,
        |worker, chunk, output| lift_chunk(worker, chunk, output, &headers, write_rejected),
        |output: Output| {
            writer
                .get_mut()
                .write_all(&output.lifted)
                .context("writing lifted VCF records")?;

            if let Some(writer) = rejected_writer.as_mut() {
                writer
                    .get_mut()
                    .write_all(&output.rejected)
                    .context("writing rejected VCF records")?;
            }

            n += output.records;
            lifted_count += output.lifted_count;
            swapped_count += output.swapped_count;

            for (rejection, count) in output.rejections {
                *rejections.entry(rejection).or_default() += count;
            }

            Ok(())
        },
    )?;

    writer.get_mut().flush()?;

    if let Some(mut writer) = rejected_writer {
        writer.get_mut().flush()?;
    }

    let rejected = n - lifted_count;
    info!(
        "lifted over {lifted_count} of {n} VCF records ({swapped_count} with swapped alleles); \
         {rejected} records were rejected"
    );

    for (rejection, count) in rejections {
        info!("  {}: {count}", rejection.code());
    }

    Ok(())
}

/// Lifts over a chunk of records.
///
/// Rejected records are only written when `write_rejected` is set.
fn lift_chunk(
    worker: &mut Worker<'_>,
    chunk: &Chunk,
    output: &mut Output,
    headers: &Headers,
    write_rejected: bool,
) -> Result<()> {
    let Worker {
        machine,
        fasta,
        record,
    } = worker;

    let mut reader = vcf::io::Reader::new(&chunk.data[..]);
    let mut writer = vcf::io::Writer::new(&mut output.lifted);
    let mut rejected_writer = vcf::io::Writer::new(&mut output.rejected);

    for n in chunk.first..chunk.first + chunk.lines as u64 {
        // NOTE: records are read into a reused buffer because `record_bufs()`
        // clones every record it reads. Parsing resets every field of the
        // buffer.
        reader
            .read_record_buf(&headers.input, record)
            .with_context(|| format!("reading VCF record {n}"))?;

        output.records += 1;

        let lifted =
            lift(machine, fasta, record).with_context(|| format!("lifting over VCF record {n}"))?;

        // NOTE: the swap is applied to a copy so that a record whose fields
        // cannot be swapped is written to the rejected file unchanged.
        let lifted = match lifted {
            Ok(lifted) if lifted.swapped => {
                let mut swapped = record.clone();

                if swap::swap_alleles(&headers.input, &mut swapped) {
                    debug!("swapped the alleles of VCF record {n}");
                    swapped
                        .info_mut()
                        .insert(String::from(SWAPPED_ALLELES_KEY), Some(Value::Flag));
                    *record = swapped;
                    Ok(lifted)
                } else {
                    debug!("could not swap the allele-indexed fields of VCF record {n}");
                    Err(Rejection::SwappedAlleles)
                }
            }
            lifted => lifted,
        };

        match lifted {
            Ok(lifted) => {
                output.lifted_count += 1;
                output.swapped_count += u64::from(lifted.swapped);

                *record.reference_sequence_name_mut() = lifted.contig.to_string();
                *record.variant_start_mut() = Some(
                    NoodlesPosition::try_from(usize::try_from(lifted.position.get())?)
                        .context("invalid lifted position")?,
                );
                *record.reference_bases_mut() = lifted.reference;
                *record.alternate_bases_mut() = AlternateBases::from(lifted.alternates);

                writer
                    .write_variant_record(&headers.output, record)
                    .with_context(|| format!("writing VCF record {n}"))?;
            }
            Err(rejection) => {
                *output.rejections.entry(rejection).or_default() += 1;

                if write_rejected {
                    record.info_mut().insert(
                        String::from(REJECT_REASON_KEY),
                        Some(Value::String(rejection.code().into())),
                    );

                    rejected_writer
                        .write_variant_record(&headers.rejected, record)
                        .with_context(|| format!("writing rejected VCF record {n}"))?;
                }
            }
        }
    }

    Ok(())
}

/// Returns whether an alternate allele is symbolic (including breakends).
fn is_symbolic(allele: &str) -> bool {
    allele.starts_with('<')
        || allele.contains(['[', ']'])
        || (allele.len() > 1 && (allele.starts_with('.') || allele.ends_with('.')))
}

/// Lifts over a variant record.
///
/// The outer result reports I/O failures, and the inner result reports
/// whether the record could be lifted over.
fn lift(
    machine: &Machine,
    fasta: &mut IndexedFasta,
    record: &RecordBuf,
) -> Result<std::result::Result<Lifted, Rejection>> {
    let Some(position) = record.variant_start() else {
        return Ok(Err(Rejection::MissingPosition));
    };

    let alternates = record.alternate_bases().as_ref();

    if alternates.iter().any(|allele| is_symbolic(allele)) {
        return Ok(Err(Rejection::SymbolicAllele));
    }

    if record.info().get(END_KEY).is_some() {
        return Ok(Err(Rejection::EndField));
    }

    let contig = record.reference_sequence_name();

    if !machine.reference_chromosomes().contains_key(contig) {
        return Ok(Err(Rejection::UnknownContig));
    }

    let reference = record.reference_bases().to_ascii_uppercase();

    // NOTE: VCF positions are 1-based and name the first base of the reference
    // allele, so one is subtracted to get the 0-based interbase start of the
    // interval covering the reference allele (whose exclusive end is the start
    // plus the allele's length). `noodles` positions are never zero, so this
    // cannot underflow.
    let start = position.get() as u64 - 1;
    let end = start.checked_add(reference.len() as u64).with_context(|| {
        format!("VCF record at {contig}:{position} ends beyond the maximum position")
    })?;

    let interval = Interval::try_new(
        Coordinate::new(Contig::new_unchecked(contig), Strand::Positive, start),
        Coordinate::new(Contig::new_unchecked(contig), Strand::Positive, end),
    )
    .with_context(|| format!("creating interval for VCF record at {contig}:{position}"))?;

    let Some(results) = machine.liftover(interval) else {
        return Ok(Err(Rejection::NoAlignment));
    };

    let [result] = results.as_slice() else {
        return Ok(Err(Rejection::MultipleMappings));
    };

    let [segment] = result.segments() else {
        return Ok(Err(Rejection::PartialMapping));
    };

    if forward_bounds(segment.reference()) != (start, end) {
        return Ok(Err(Rejection::PartialMapping));
    }

    let query_contig = segment.query().contig();
    let (query_start, query_end) = forward_bounds(segment.query());
    let flipped = segment.query().strand() != segment.reference().strand();

    let Some(target) = fasta.fetch(query_contig, query_start, query_end)? else {
        return Ok(Err(Rejection::ContigNotInFasta));
    };

    let transformed = transform(
        &reference,
        alternates,
        &target,
        flipped,
        || match query_start.checked_sub(1) {
            Some(anchor_start) => fasta.fetch(query_contig, anchor_start, query_start),
            None => Ok(None),
        },
    )?;

    let Alleles {
        offset,
        reference,
        alternates,
        swapped,
    } = match transformed {
        Ok(alleles) => alleles,
        Err(rejection) => return Ok(Err(rejection)),
    };

    // NOTE: the interbase start is converted to a 1-based (in-base) position,
    // shifted back by any base used to re-anchor an indel.
    let position = Position::try_from(query_start + 1 - offset)
        .context("lifted position is not a valid 1-based position")?;

    Ok(Ok(Lifted {
        contig: query_contig.clone(),
        position,
        reference,
        alternates,
        swapped,
    }))
}

/// Transforms the alleles of a record that maps to `target` (the query
/// genome's forward-strand bases where the reference allele maps).
///
/// The alleles of a biallelic single-base substitution whose alternate allele
/// matches `target` are swapped.
fn transform(
    reference: &str,
    alternates: &[String],
    target: &[u8],
    flipped: bool,
    fetch_anchor: impl FnOnce() -> Result<Option<Vec<u8>>>,
) -> Result<std::result::Result<Alleles, Rejection>> {
    let orient = |allele: &str| -> String {
        if !flipped || allele == OVERLAPPING_DELETION {
            return allele.to_string();
        }

        String::from_utf8(reverse_complement(allele.as_bytes()))
            .expect("reverse complement of a string should be valid UTF-8")
    };

    let lifted_reference = orient(reference);
    let mut alternates = alternates.iter().map(|a| orient(a)).collect::<Vec<_>>();
    let target = String::from_utf8_lossy(target).into_owned();

    let is_substitution = reference.len() == 1
        && alternates
            .iter()
            .all(|a| a.len() == 1 || a == OVERLAPPING_DELETION);

    if !lifted_reference.eq_ignore_ascii_case(&target) {
        let is_swapped =
            is_substitution && alternates.iter().any(|a| a.eq_ignore_ascii_case(&target));

        if !is_swapped {
            return Ok(Err(Rejection::ReferenceMismatch));
        }

        // NOTE: only biallelic records are swapped, as the allele-indexed
        // fields of multiallelic records would need to be reordered rather
        // than simply exchanged.
        let [alternate] = alternates.as_mut_slice() else {
            return Ok(Err(Rejection::SwappedAlleles));
        };

        let reference = std::mem::replace(alternate, lifted_reference);

        return Ok(Ok(Alleles {
            offset: 0,
            reference,
            alternates,
            swapped: true,
        }));
    }

    let reference = lifted_reference;

    let is_indel = alternates
        .iter()
        .any(|a| a != OVERLAPPING_DELETION && a.len() != reference.len());

    if !(flipped && is_indel) {
        return Ok(Ok(Alleles {
            offset: 0,
            reference,
            alternates,
            swapped: false,
        }));
    }

    // NOTE: VCF indels are usually padded with the base before the event.
    // After reverse complementing, that padding base is at the end of each
    // allele, so every allele is anchored with the base before the event in
    // the query genome, and a trailing base shared by all alleles is then
    // removed. Indels without a shared padding base (e.g., `AC>G`) keep all
    // of their bases.
    let Some(anchor) = fetch_anchor()? else {
        return Ok(Err(Rejection::MissingAnchor));
    };

    let anchor = String::from_utf8_lossy(&anchor).into_owned();
    let padding = reference.as_bytes()[reference.len() - 1];
    let is_padded = alternates.iter().all(|a| {
        a == OVERLAPPING_DELETION
            || a.as_bytes()
                .last()
                .is_some_and(|b| b.eq_ignore_ascii_case(&padding))
    });

    let reanchor = |allele: &str| {
        let allele = if is_padded {
            &allele[..allele.len() - 1]
        } else {
            allele
        };

        format!("{anchor}{allele}")
    };

    let reference = reanchor(&reference);
    for allele in alternates.iter_mut() {
        if allele != OVERLAPPING_DELETION {
            *allele = reanchor(allele);
        }
    }

    Ok(Ok(Alleles {
        offset: 1,
        reference,
        alternates,
        swapped: false,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles(values: &[&str]) -> Vec<String> {
        values.iter().map(|v| v.to_string()).collect()
    }

    fn lifted(
        offset: u64,
        reference: &str,
        alternates: &[&str],
        swapped: bool,
    ) -> std::result::Result<Alleles, Rejection> {
        Ok(Alleles {
            offset,
            reference: reference.to_string(),
            alternates: alleles(alternates),
            swapped,
        })
    }

    fn run(
        reference: &str,
        alternates: &[&str],
        target: &str,
        flipped: bool,
        anchor: Option<&str>,
    ) -> std::result::Result<Alleles, Rejection> {
        transform(
            reference,
            &alleles(alternates),
            target.as_bytes(),
            flipped,
            || Ok(anchor.map(|a| a.as_bytes().to_vec())),
        )
        .unwrap()
    }

    #[test]
    fn substitutions() {
        assert_eq!(
            run("T", &["G"], "T", false, None),
            lifted(0, "T", &["G"], false)
        );
        assert_eq!(
            run("t", &["g"], "T", false, None),
            lifted(0, "t", &["g"], false)
        );
        assert_eq!(
            run("T", &["C"], "A", true, None),
            lifted(0, "A", &["G"], false)
        );
        assert_eq!(
            run("G", &["T", "*"], "A", false, None),
            Err(Rejection::ReferenceMismatch)
        );
        assert_eq!(
            run("AG", &["TT"], "CA", false, None),
            Err(Rejection::ReferenceMismatch)
        );
    }

    #[test]
    fn swapped_substitutions() {
        assert_eq!(
            run("A", &["C"], "C", false, None),
            lifted(0, "C", &["A"], true)
        );
        assert_eq!(
            run("T", &["G"], "C", true, None),
            lifted(0, "C", &["A"], true)
        );
        assert_eq!(
            run("A", &["C", "G"], "C", false, None),
            Err(Rejection::SwappedAlleles)
        );
    }

    #[test]
    fn indels() {
        assert_eq!(
            run("ACA", &["A"], "ACA", false, None),
            lifted(0, "ACA", &["A"], false)
        );
        assert_eq!(
            run("AC", &["A"], "GT", true, Some("C")),
            lifted(1, "CG", &["C"], false)
        );
        assert_eq!(
            run("T", &["TGG", "*"], "A", true, Some("G")),
            lifted(1, "G", &["GCC", "*"], false)
        );
        assert_eq!(
            run("AC", &["G"], "AC", false, None),
            lifted(0, "AC", &["G"], false)
        );
        assert_eq!(
            run("AC", &["G"], "GT", true, Some("C")),
            lifted(1, "CGT", &["CC"], false)
        );
        assert_eq!(
            run("A", &["TG", "*"], "T", true, Some("C")),
            lifted(1, "CT", &["CCA", "*"], false)
        );
        assert_eq!(
            run("AC", &["A"], "GT", true, None),
            Err(Rejection::MissingAnchor)
        );
        assert_eq!(
            run("AC", &["A"], "GG", true, Some("C")),
            Err(Rejection::ReferenceMismatch)
        );
    }
}
