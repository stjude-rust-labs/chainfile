//! A binary to comprehensively test that the `chainfile` crate is performing
//! liftover comparably to other liftover utilities.
//!
//! ```shell
//! cargo run --release --bin=compare-crossmap --features=binaries hg19ToHg38
//! ```
//!
//! It achieves this by carrying out the following:
//!
//! * Randomly generating `n` singular locations within a genome, performing
//!   liftover with the facilities provided by `chainfile` and [`CrossMap`], and
//!   ensuring all of the results match.
//!
//! [`CrossMap`]: https://crossmap.sourceforge.net/

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Write as _;
use std::ops::Deref;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::process::Stdio;
use std::sync::LazyLock;

use anyhow::Context;
use anyhow::Result;
use anyhow::bail;
use chainfile::alignment::section::header::sequence;
use chainfile::liftover;
use clap::Parser;
use clap_verbosity_flag::Verbosity;
use flate2::read::GzDecoder;
use omics::coordinate::Coordinate;
use omics::coordinate::Interval;
use omics::coordinate::Strand;
use omics::coordinate::position::Number;
use omics::coordinate::system::Interbase;
use rand::Rng;
use rand::rngs::ThreadRng;
use regex::Regex;
use tempdir::TempDir;
use tracing::error;
use tracing::info;
use tracing::warn;
use tracing_log::AsTrace as _;
use tracing_subscriber::EnvFilter;
use weighted_rand::builder::NewBuilder;
use weighted_rand::builder::WalkerTableBuilder;
use weighted_rand::table::WalkerTable;

const CHAINFILE_URL_PREFIX: &str = "https://hgdownload.soe.ucsc.edu/goldenPath";

////////////////////////////////////////////////////////////////////////////////////////
// Chain names
////////////////////////////////////////////////////////////////////////////////////////

pub struct ChainName {
    /// The genome we're converting from.
    ///
    /// This genome name will always be lowercase as per UCSC conventions.
    from: String,

    /// The genome we're converting To.
    ///
    /// This genome name will always be sentence case as per UCSC conventions.
    to: String,
}

static REGEX: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"([a-z0-9_]+)To([A-Z][a-z0-9_]*)").unwrap());

impl ChainName {
    /// Attempts to create a new chain name.
    ///
    /// [`None`] is returned if the chain name is not valid.
    pub fn try_new(value: impl AsRef<str>) -> Option<Self> {
        let value = value.as_ref();
        let groups = REGEX.captures(value)?;

        let from = groups.get(1).unwrap().as_str().to_string();
        let to = groups.get(2).unwrap().as_str().to_string();

        Some(Self { from, to })
    }

    /// Gets the "from" genome name.
    pub fn from(&self) -> &str {
        &self.from
    }

    /// Gets the "to" genome name.
    pub fn to(&self) -> &str {
        &self.to
    }

    /// Gets the full chain name.
    pub fn name(&self) -> String {
        format!("{}To{}", self.from, self.to)
    }

    /// Gets the full file name for the chain.
    pub fn file_name(&self) -> String {
        format!("{}.over.chain.gz", self.name())
    }

    /// Gets the download url for the file.
    pub fn download_url(&self) -> String {
        format!(
            "{}/{}/liftOver/{}",
            CHAINFILE_URL_PREFIX,
            self.from,
            self.file_name()
        )
    }
}

#[cfg(test)]
mod chain_name_tests {
    use super::ChainName;

    #[test]
    fn valid() {
        let name = ChainName::try_new("hg19ToHg38").unwrap();

        assert_eq!(name.from(), "hg19");
        assert_eq!(name.to(), "Hg38");
        assert_eq!(name.name(), "hg19ToHg38");
        assert_eq!(name.file_name(), "hg19ToHg38.over.chain.gz");
        assert_eq!(
            name.download_url(),
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
        );
    }

    #[test]
    fn invalid() {
        assert!(ChainName::try_new("hg19Tohg38").is_none());
        assert!(ChainName::try_new("hg19To").is_none());
        assert!(ChainName::try_new("hg19").is_none());
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Bed files
////////////////////////////////////////////////////////////////////////////////////////

/// A simple struct representing a single entry in a BED file.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct BedEntry {
    /// The contig.
    contig: String,

    /// The start position.
    start: Number,

    /// The end position.
    end: Number,
}

impl BedEntry {
    /// Creates a new, single-position BED entry.
    fn new_single_position(contig: String, position: Number) -> Self {
        // NOTE: BED wants a 0-based, inclusive start and a 1-based, exclusive
        // end. To make that happen here, we do this strange conversion below.
        Self {
            contig,
            start: position,
            end: position + 1,
        }
    }

    /// Turns the BED entry back into an interbase coordinate.
    fn into_coordinate(self, strand: Strand) -> Coordinate<Interbase> {
        omics::coordinate::Coordinate::<Interbase>::new(
            self.contig.as_str(),
            strand,
            self.start as Number,
        )
    }
}

fn write_bed_file(path: &Path, positions: impl IntoIterator<Item = BedEntry>) -> Result<()> {
    let mut file = File::create(path).context("creating BED file")?;

    for position in positions {
        writeln!(
            file,
            "{}\t{}\t{}",
            position.contig, position.start, position.end
        )
        .context("writing to BED file")?;
    }

    Ok(())
}

impl std::fmt::Display for BedEntry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} (BED format)",
            self.contig, self.start, self.end
        )
    }
}

/// The result of a liftover operation.
#[derive(Clone, Debug, Eq, PartialEq)]
enum LiftoverResult {
    /// The segment was lifted over successfully.
    Mapped(BedEntry),

    /// The segment could not be lifted over.
    Unmapped,
}

impl std::fmt::Display for LiftoverResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LiftoverResult::Mapped(entry) => write!(f, "{entry}"),
            LiftoverResult::Unmapped => write!(f, "<unmapped>"),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// CrossMap execution
////////////////////////////////////////////////////////////////////////////////////////

/// A struct for working with `CrossMap`.
struct CrossMap;

impl CrossMap {
    /// Ensures the `CrossMap` is installed.
    fn ensure_installed() -> Result<()> {
        Command::new("CrossMap")
            .arg("--version")
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .status()
            .context("retrieving `CrossMap`'s version")
            .map(|_| ())
    }

    /// Runs `crossmap bed` on a chain file and bed file.
    fn run_bed(chain_file: &Path, bed_file: &Path) -> Result<HashMap<BedEntry, LiftoverResult>> {
        let process = Command::new("CrossMap")
            .arg("bed")
            .arg(chain_file)
            .arg(bed_file)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()?;

        let mut results = HashMap::new();

        let output = process
            .wait_with_output()
            .context("running `CrossMap bed`")?;

        if !output.status.success() {
            bail!("the `CrossMap bed` command returned a non-zero exit code");
        }

        if !output.stderr.is_empty() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            for line in stderr.lines() {
                if !line.contains("[INFO]") {
                    panic!("`CrossMap` returned something on stderr: {stderr}");
                }
            }
        }

        let stdout =
            String::from_utf8(output.stdout).context("decoding stdout from byte stream")?;

        for line in stdout.lines() {
            let parts = line.split("\t").collect::<Vec<_>>();
            assert!(
                parts.len() > 3,
                "should always have at least four parts of a `CrossMap` result line; found \
                 `{line}`"
            );

            let from = BedEntry {
                contig: parts[0].to_owned(),
                start: parts[1]
                    .parse()
                    .context("parsing the from start position")?,
                end: parts[2].parse().context("parsing the from end position")?,
            };

            let to = match parts.len() {
                4 => LiftoverResult::Unmapped,
                7 => {
                    assert!(
                        parts[3] == "->",
                        "the fourth column should always be an arrow"
                    );

                    let to = BedEntry {
                        contig: parts[4].to_owned(),
                        start: parts[5].parse().context("parsing the to start position")?,
                        end: parts[6].parse().context("parsing the to end position")?,
                    };

                    LiftoverResult::Mapped(to)
                }
                _ => unreachable!("cannot parse `CrossMap` result: `{line}`"),
            };

            results.insert(from, to);
        }

        Ok(results)
    }
}
////////////////////////////////////////////////////////////////////////////////////////
// Choosing random locations
////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug)]
struct Chromosomes {
    /// The inner contigs.
    chromosomes: Box<[(String, u32)]>,

    /// Whether or not to only generate random positions from the canonical
    /// chromsomes.
    canonical_chromosomes_only: bool,

    /// The weighted distribution.
    weights: WalkerTable,

    /// The random number generator.
    rng: ThreadRng,
}

const CANONICAL_CHROMOSOMES: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY", "chrM", "chrMT", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M", "MT",
];

impl Chromosomes {
    /// Creates a new [`Chromosomes`].
    fn new(corpus: &HashMap<String, Number>, canonical_chromosomes_only: bool) -> Self {
        let chromosomes = corpus
            .iter()
            .map(|(chromosome, size)| {
                let size: u32 = (*size).try_into().expect("chromosome size to fit into u32");
                (chromosome.to_owned(), size)
            })
            .collect::<Vec<_>>()
            .into_boxed_slice();

        let weights = chromosomes
            .iter()
            .map(|(_, weight)| *weight)
            .collect::<Vec<_>>();

        Self {
            chromosomes,
            canonical_chromosomes_only,
            weights: WalkerTableBuilder::new(&weights).build(),
            rng: Default::default(),
        }
    }

    /// Picks a weighted random position from the contigs based on each contig's
    /// length.
    fn random_position(&mut self) -> BedEntry {
        loop {
            let index = self.weights.next();

            // SAFETY: we just picked a random contig from the map, so we know this
            // exists and will always unwrap.
            let (contig, size) = self.chromosomes[index].clone();

            if !self.canonical_chromosomes_only || CANONICAL_CHROMOSOMES.contains(&contig.as_str())
            {
                let position = self.rng.gen_range(0..size);
                return BedEntry::new_single_position(contig, position as Number);
            }
            // If the contig is not valid, continue to retry with a new index
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////
// Working directory
////////////////////////////////////////////////////////////////////////////////////////

struct WorkDirectory(PathBuf);

impl WorkDirectory {
    /// Creates a new [`WorkDirectory`] from the provided directory.
    fn new(directory: PathBuf) -> Result<Self> {
        let work_dir = Self(directory);

        std::fs::create_dir_all(work_dir.reference_dir())
            .context("creating reference directory")?;

        Ok(work_dir)
    }

    /// Creates a new [`WorkDirectory`] from a temporary directory generated
    /// only for this execution.
    fn new_temporary() -> Result<Self> {
        Self::new(
            TempDir::new("work")
                .context("creating (temporary) working directory")?
                .into_path(),
        )
    }

    /// The reference directory.
    fn reference_dir(&self) -> PathBuf {
        self.0.join("reference")
    }

    /// The input positions bed file to `CrossMap`.
    fn input_bed_file(&self) -> PathBuf {
        self.0.join("input.bed")
    }

    /// Attempts to download a chain file.
    fn download_chain_file(&self, chain_name: &ChainName) -> Result<PathBuf> {
        assert!(
            self.reference_dir().is_dir(),
            "reference directory must exist before downloading chain file"
        );

        let filename = chain_name.file_name();
        let filepath = self.reference_dir().join(&filename);

        if filepath.exists() {
            info!("file path already exists: {}!", filepath.display());
            info!("assuming this file is correct and eliding downloading");
            return Ok(filepath);
        }

        let url = chain_name.download_url();

        info!("chain file: dowloading {url}");

        let mut out = File::create(&filepath).context("creating new chain file file")?;
        let resp = reqwest::blocking::get(url).context("downloading chain file")?;

        if resp.status().is_success() {
            let bytes = resp.bytes().context("parsing chain file bytes")?;
            std::io::copy(&mut bytes.as_ref(), &mut out).context("writing chain file to disk")?;
            info!("chain file: download completed");
        } else {
            error!(
                "chain file: download failed with status `{}`",
                resp.status()
            );
        }

        Ok(filepath)
    }

    /// Write a set of positions to a BED file.
    fn write_input_bed_file(&self, positions: impl IntoIterator<Item = BedEntry>) -> Result<()> {
        write_bed_file(&self.input_bed_file(), positions)
    }
}

impl From<PathBuf> for WorkDirectory {
    fn from(value: PathBuf) -> Self {
        Self(value)
    }
}

impl Deref for WorkDirectory {
    type Target = PathBuf;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Comparsion between `chainfile` and`CrossMap`
////////////////////////////////////////////////////////////////////////////////////////

/// A comparison of liftover results between `chainfile` and `CrossMap`.
#[derive(Clone, Debug)]
struct Comparison {
    /// The original position.
    from: BedEntry,

    /// The lifted position from `chainfile`.
    chainfile: LiftoverResult,

    /// The lifted position from `CrossMap`.
    crossmap: LiftoverResult,
}

impl Comparison {
    /// Creates a new [`Comparison`].
    fn new(from: BedEntry, chainfile: LiftoverResult, crossmap: LiftoverResult) -> Self {
        Self {
            from,
            chainfile,
            crossmap,
        }
    }

    /// Returns whether or not the results _exactly_ matched.
    fn matched(&self) -> bool {
        self.chainfile == self.crossmap
    }

    fn from(&self) -> &BedEntry {
        &self.from
    }

    fn chainfile(&self) -> &LiftoverResult {
        &self.chainfile
    }

    fn crossmap(&self) -> &LiftoverResult {
        &self.crossmap
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////////////

/// Throws the proverbial kitchen sink at `chainfile`.
#[derive(Parser)]
struct Args {
    /// The name of the chain file to compare (e.g., `hg19ToHg38`—don't include
    /// the `.over.chain.gz`).
    chain: String,

    /// Whether or not to only look at the so-called "canonical" contigs
    /// (chr1-22, chrX, chrY, and chrM).
    #[arg(short, long, default_value_t = false)]
    canonical_chromosomes_only: bool,

    /// If desired, a permanent directory within which to do the work. This is
    /// generally for debugging only.
    #[arg(short, long)]
    directory: Option<PathBuf>,

    /// If desired, the number of mismatches to explore.
    #[arg(short, long)]
    explore_mismatches: Option<usize>,

    /// The number of positions to generate.
    #[arg(short, default_value_t = 1_000_000)]
    n: usize,

    #[command(flatten)]
    verbose: Verbosity,
}

fn throw(args: &Args) -> Result<()> {
    let chain_name = match ChainName::try_new(args.chain.clone()) {
        Some(name) => name,
        None => bail!("invalid chain name: {}", args.chain),
    };

    let work_dir = match &args.directory {
        Some(dir) => WorkDirectory::new(dir.clone()),
        None => WorkDirectory::new_temporary(),
    }
    .context("creating working directory")?;

    info!("working directory: {}", work_dir.display());

    info!("crossmap: ensuring `CrossMap` is installed");
    CrossMap::ensure_installed().context("ensuring `CrossMap` is installed")?;

    let chain_file_path = work_dir
        .download_chain_file(&chain_name)
        .with_context(|| format!("chain file: downloading {}", &chain_name.name()))?;

    let chain = File::open(&chain_file_path)
        .map(GzDecoder::new)
        .map(BufReader::new)
        .map(chainfile::Reader::new)?;

    let machine = liftover::machine::builder::Builder.try_build_from(chain)?;
    let query_contigs =
        Chromosomes::new(machine.query_chromosomes(), args.canonical_chromosomes_only);
    let mut reference_contigs = Chromosomes::new(
        machine.reference_chromosomes(),
        args.canonical_chromosomes_only,
    );

    let mut locations = Vec::new();

    for _ in 0..args.n {
        locations.push(reference_contigs.random_position());
    }

    info!("crossmap: writing BED file of positions");
    work_dir
        .write_input_bed_file(locations)
        .context("writing BED file")?;

    let mut comparisons = Vec::new();

    info!("crossmap: running liftover");
    for (from, crossmap_result) in CrossMap::run_bed(&chain_file_path, &work_dir.input_bed_file())
        .context("getting the `CrossMap` mappings")?
    {
        let start =
            Coordinate::<Interbase>::new(from.contig.as_str(), Strand::Positive, from.start);
        let end = Coordinate::<Interbase>::new(from.contig.as_str(), Strand::Positive, from.end);
        let from_interval =
            Interval::try_new(start, end).expect("interval to be able to be created");

        let chainfile_result = machine
            .liftover(from_interval)
            .map(|mut results| {
                assert!(
                    results.len() == 1,
                    "single valued inputs should only ever produce one output contiguous interval \
                     pair"
                );

                // SAFETY: we just asserted that the length is one, so this will
                // always unwrap.
                let result = results.pop().unwrap();
                assert!(
                    result.reference().strand() == Strand::Positive,
                    "strand should always be positive for reference"
                );

                let query = result.into_query();

                // NOTE: `CrossMap` always reports the results on the positive
                // strand. This alters how the nucleotide is represented as a
                // BED coordinate, as the start position on the positive strand
                // is the end position on the negative strand. Thus, we need to
                // reverse complement the results before comparison.
                let query = match query.strand() {
                    Strand::Positive => query,
                    Strand::Negative => query.reverse_complement(),
                };

                let (contig, _, position) = query.into_start().into_parts();

                LiftoverResult::Mapped(BedEntry::new_single_position(
                    contig.into_inner(),
                    // SAFETY: this should always unwrap, as the lower bound isn't
                    // going to be used in this kind of thing.
                    position.get() as Number,
                ))
            })
            .unwrap_or(LiftoverResult::Unmapped);

        comparisons.push(Comparison::new(from, chainfile_result, crossmap_result));
    }

    let total_comparisons = comparisons.len();
    let mismatches = comparisons
        .into_iter()
        .filter(|comparison| !comparison.matched())
        .collect::<Vec<_>>();
    let n_matches = total_comparisons - mismatches.len();
    let matched_percent = (n_matches as f64 * 100.0) / total_comparisons as f64;

    println!(
        "`chainfile` and `CrossMap` matched in {matched_percent:.2}% of cases (n = \
         {total_comparisons})"
    );

    if args.explore_mismatches.is_some() && !mismatches.is_empty() {
        let sections = File::open(&chain_file_path)
            .map(GzDecoder::new)
            .map(BufReader::new)
            .map(chainfile::Reader::new)?
            .sections()
            .collect::<Result<Vec<_>, _>>()
            .context("rereading the chainfile sections")?
            .into_iter()
            .map(|section| {
                section
                    .reference_sequence()
                    .interval()
                    .map(|interval| (section, interval))
            })
            .collect::<Result<Vec<_>, sequence::Error>>()
            .context("creating intervals from chainfile")?;

        // SAFETY: we just checked about that this is [`Some`] in the `if`
        // statement, so this will always unwrap.
        let take_n_mismatches = args.explore_mismatches.unwrap();

        for (i, comparison) in mismatches.into_iter().enumerate().take(take_n_mismatches) {
            warn!("== unmatched example #{} ==", i + 1);
            warn!("reference: {}", comparison.from());
            warn!("chainfile: {}", comparison.chainfile());
            warn!("crossmap:  {}", comparison.crossmap());
            warn!("  ↳ relevant alignment sections:");

            let coordinate = comparison
                .from()
                .clone()
                .into_coordinate(Strand::Positive)
                .nudge_forward()
                .unwrap();

            for (section, interval) in &sections {
                if interval.contains_entity(&coordinate) {
                    warn!("    ↳ {}", section.header());

                    for result in section.stepthrough().expect("step-through to be created") {
                        let pairs = result?;

                        if pairs.reference().contains_entity(&coordinate) {
                            warn!("      ↳ {}", pairs);

                            let query = pairs.into_query();
                            if query.strand() == Strand::Negative {
                                // (1) Find the size of the chromosome the query
                                // sits on.
                                let (_, size) = query_contigs
                                    .chromosomes
                                    .iter()
                                    .find(|(chromosome, _)| chromosome == query.contig().as_str())
                                    .expect("this should be found");

                                // (2) Grab the existing coordinate parts.
                                let (old_start, old_end) = query.into_coordinates();
                                let (old_start_contig, old_start_strand, old_start_pos) =
                                    old_start.into_parts();
                                let (old_end_contig, old_end_strand, old_end_pos) =
                                    old_end.into_parts();

                                // (3) Invert them using the current chromosome
                                // size.

                                // SAFETY: this should always unwrap, as we
                                // constructed this coordinate on the opposite
                                // strand successfully (and none of the
                                // operations here should cause it to not
                                // construct—unless, of course, there is a bug).
                                let new_start = Coordinate::<Interbase>::new(
                                    old_start_contig,
                                    old_start_strand.complement(),
                                    *size as Number - old_start_pos.get(),
                                );

                                // SAFETY: this should always unwrap, as we
                                // constructed this coordinate on the opposite
                                // strand successfully (and none of the
                                // operations here should cause it to not
                                // construct—unless, of course, there is a bug).
                                let new_end = Coordinate::<Interbase>::new(
                                    old_end_contig,
                                    old_end_strand.complement(),
                                    *size as Number - old_end_pos.get(),
                                );

                                // SAFETY: this should always unwrap, as we
                                // constructed this interval on the opposite
                                // strand successfully (and none of the
                                // operations here should cause it to not
                                // construct—unless, of course, there is a bug).
                                let interval = Interval::try_new(new_start, new_end).unwrap();
                                warn!("        ↳ in other words, {}", interval);
                            }
                        };
                    }
                }
            }

            if i < take_n_mismatches - 1 {
                warn!("");
            }
        }
    }

    if matched_percent != 100.0 {
        std::process::exit(1);
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    assert!(args.n > 0, "`n` must be greater than 0!");

    match std::env::var("RUST_LOG") {
        Ok(_) => tracing_subscriber::fmt()
            .with_env_filter(EnvFilter::from_default_env())
            .init(),
        Err(_) => tracing_subscriber::fmt()
            .with_max_level(args.verbose.log_level_filter().as_trace())
            .init(),
    };

    throw(&args)
}
