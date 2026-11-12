//! Shared state and helpers for comparisons.

use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write as _;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use rand::RngExt as _;
use rand::SeedableRng as _;
use rand::rngs::StdRng;
use tempfile::TempDir;

/// The directory in which files are written.
enum Directory {
    /// A directory that is kept afterward.
    Kept(PathBuf),

    /// A temporary directory that is removed afterward.
    Temporary(TempDir),
}

/// Shared state for a comparison.
pub struct Context {
    /// The random number generator.
    pub rng: StdRng,

    /// The `chainfile` executable.
    chainfile: PathBuf,

    /// The directory in which files are written.
    directory: Directory,

    /// The maximum number of disagreements to print.
    pub examples: usize,
}

impl Context {
    /// Creates a new context.
    pub fn new(
        seed: Option<u64>,
        chainfile: PathBuf,
        work_dir: Option<PathBuf>,
        examples: usize,
    ) -> Result<Self> {
        let seed = seed.unwrap_or_else(|| rand::rng().random());
        eprintln!("seed: {seed}");

        let directory = match work_dir {
            Some(path) => {
                fs::create_dir_all(&path)
                    .with_context(|| format!("creating `{}`", path.display()))?;
                Directory::Kept(path)
            }
            None => Directory::Temporary(TempDir::new().context("creating temporary directory")?),
        };

        Ok(Self {
            rng: StdRng::seed_from_u64(seed),
            chainfile,
            directory,
            examples,
        })
    }

    /// Gets the path of a file within the working directory.
    pub fn path(&self, name: &str) -> PathBuf {
        match &self.directory {
            Directory::Kept(path) => path.join(name),
            Directory::Temporary(directory) => directory.path().join(name),
        }
    }

    /// Creates a command that runs `chainfile`.
    pub fn chainfile(&self) -> Command {
        Command::new(&self.chainfile)
    }
}

/// Creates a file and writes to it.
pub fn create(path: &Path, write: impl FnOnce(&mut BufWriter<File>) -> Result<()>) -> Result<()> {
    let file = File::create(path).with_context(|| format!("creating `{}`", path.display()))?;
    let mut writer = BufWriter::new(file);
    write(&mut writer).with_context(|| format!("writing `{}`", path.display()))?;
    writer.flush()?;
    Ok(())
}

/// Runs a command and fails if it does not succeed.
pub fn execute(command: &mut Command) -> Result<()> {
    eprintln!("running {command:?}");

    let status = command
        .status()
        .with_context(|| format!("running {command:?}"))?;

    if !status.success() {
        bail!("{command:?} failed ({status})");
    }

    Ok(())
}

/// Displays a path for a report title.
pub fn display(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}
