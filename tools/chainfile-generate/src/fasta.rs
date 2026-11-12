//! Indexed FASTA access.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;
use noodles::core::Position;
use noodles::core::Region;
use noodles::fasta;

/// An indexed FASTA file.
///
/// The FASTA file must have a `.fai` index. BGZF-compressed FASTA files
/// additionally require a `.gzi` index.
pub struct IndexedFasta {
    /// The inner reader.
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,

    /// The contig lengths keyed by name.
    lengths: HashMap<String, u64>,
}

impl IndexedFasta {
    /// Opens an indexed FASTA file.
    pub fn open(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| {
                format!(
                    "opening indexed FASTA `{}` (a `.fai` index is required)",
                    path.display()
                )
            })?;

        let lengths = reader
            .index()
            .as_ref()
            .iter()
            .map(|record| (record.name().to_string(), record.length()))
            .collect();

        Ok(Self { reader, lengths })
    }

    /// Fetches the uppercased bases within the 0-based, half-open interval
    /// `[start, end)` on the forward strand of a contig.
    pub fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let Some(&length) = self.lengths.get(contig) else {
            bail!("contig `{contig}` is not present in the FASTA file");
        };

        if start > end || end > length {
            bail!(
                "interval [{start}, {end}) is out of bounds for contig `{contig}` with length \
                 {length}"
            );
        }

        if start == end {
            return Ok(Vec::new());
        }

        let region = Region::new(contig, to_position(start + 1)?..=to_position(end)?);
        let record = self
            .reader
            .query(&region)
            .with_context(|| format!("fetching {contig}:{}-{end}", start + 1))?;

        let mut bases = record.sequence().as_ref().to_vec();
        bases.make_ascii_uppercase();

        if bases.len() as u64 != end - start {
            bail!(
                "fetching {contig}:{}-{end} returned {} bases instead of {}; is the FASTA index \
                 up to date?",
                start + 1,
                bases.len(),
                end - start
            );
        }

        Ok(bases)
    }
}

/// Converts a 1-based position to a [`Position`].
fn to_position(value: u64) -> Result<Position> {
    usize::try_from(value)
        .ok()
        .and_then(Position::new)
        .with_context(|| format!("invalid position: {value}"))
}
