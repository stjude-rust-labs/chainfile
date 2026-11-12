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

/// The number of bases cached on either side of a fetched interval.
///
/// Records are typically sorted, so nearby records (in either direction, as
/// chains can map to the opposite strand) are fetched from the cache.
const WINDOW_FLANK: u64 = 32 * 1024;

/// A cached window of a contig's uppercased bases.
struct Window {
    /// The contig name.
    contig: String,

    /// The 0-based start of the window.
    start: u64,

    /// The uppercased bases.
    bases: Vec<u8>,
}

impl Window {
    /// Gets the bases within `[start, end)` if they are in the window.
    fn get(&self, contig: &str, start: u64, end: u64) -> Option<&[u8]> {
        if self.contig != contig || start < self.start {
            return None;
        }

        let from = usize::try_from(start - self.start).ok()?;
        let to = usize::try_from(end - self.start).ok()?;
        self.bases.get(from..to)
    }
}

/// An indexed FASTA file.
///
/// The FASTA file must have a `.fai` index. BGZF-compressed FASTA files
/// (`.gz` or `.bgz`) additionally require a `.gzi` index.
pub struct IndexedFasta {
    /// The inner reader.
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,

    /// The contig names and lengths in index order.
    contigs: Vec<(String, u64)>,

    /// The contig lengths keyed by name.
    lengths: HashMap<String, u64>,

    /// The most recently read window.
    window: Option<Window>,
}

impl IndexedFasta {
    /// Opens an indexed FASTA file.
    pub fn open(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| {
                format!(
                    "opening indexed FASTA `{}` (a `.fai` index is required, plus a `.gzi` index \
                     for BGZF-compressed files; see `samtools faidx`)",
                    path.display()
                )
            })?;

        let contigs = reader
            .index()
            .as_ref()
            .iter()
            .map(|record| (record.name().to_string(), record.length()))
            .collect::<Vec<_>>();

        let lengths = contigs.iter().cloned().collect();

        Ok(Self {
            reader,
            contigs,
            lengths,
            window: None,
        })
    }

    /// Gets the contig names and lengths in index order.
    pub fn contigs(&self) -> &[(String, u64)] {
        &self.contigs
    }

    /// Fetches the uppercased bases within the 0-based, half-open interval
    /// `[start, end)` on the forward strand of a contig.
    ///
    /// Returns [`None`] if the contig is not present in the FASTA file.
    pub fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<Option<Vec<u8>>> {
        let Some(&length) = self.lengths.get(contig) else {
            return Ok(None);
        };

        if start > end || end > length {
            bail!(
                "interval [{start}, {end}) is out of bounds for contig `{contig}` with length \
                 {length}"
            );
        }

        if start == end {
            return Ok(Some(Vec::new()));
        }

        if let Some(bases) = self.window.as_ref().and_then(|w| w.get(contig, start, end)) {
            return Ok(Some(bases.to_vec()));
        }

        let window_start = start.saturating_sub(WINDOW_FLANK);
        let window_end = end.saturating_add(WINDOW_FLANK).min(length);

        let region = Region::new(
            contig,
            to_position(window_start + 1)?..=to_position(window_end)?,
        );
        let record = self
            .reader
            .query(&region)
            .with_context(|| format!("fetching {contig}:{}-{window_end}", window_start + 1))?;

        let mut bases = record.sequence().as_ref().to_vec();
        bases.make_ascii_uppercase();

        if bases.len() as u64 != window_end - window_start {
            bail!(
                "fetching {contig}:{}-{window_end} returned {} bases instead of {}; is the FASTA \
                 index up to date?",
                window_start + 1,
                bases.len(),
                window_end - window_start
            );
        }

        let window = self.window.insert(Window {
            contig: contig.to_string(),
            start: window_start,
            bases,
        });

        // SAFETY: the window was just read to cover `[start, end)`.
        Ok(Some(window.get(contig, start, end).unwrap().to_vec()))
    }
}

/// Converts a 1-based position to a [`Position`].
fn to_position(value: u64) -> Result<Position> {
    usize::try_from(value)
        .ok()
        .and_then(Position::new)
        .with_context(|| format!("invalid position: {value}"))
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;

    use super::*;

    /// The number of bases per line in the test FASTA.
    const LINE_WIDTH: usize = 60;

    /// Writes an indexed FASTA with the given contigs to `dir`.
    fn write_fasta(dir: &Path, contigs: &[(&str, &[u8])]) -> std::path::PathBuf {
        let mut fasta = Vec::new();
        let mut index = String::new();

        for (name, bases) in contigs {
            fasta.extend_from_slice(format!(">{name}\n").as_bytes());
            let offset = fasta.len();

            for line in bases.chunks(LINE_WIDTH) {
                fasta.extend_from_slice(line);
                fasta.push(b'\n');
            }

            writeln!(
                index,
                "{name}\t{}\t{offset}\t{LINE_WIDTH}\t{}",
                bases.len(),
                LINE_WIDTH + 1
            )
            .unwrap();
        }

        let path = dir.join("test.fa");
        std::fs::write(&path, fasta).unwrap();
        std::fs::write(dir.join("test.fa.fai"), index).unwrap();
        path
    }

    #[test]
    fn fetches_match_the_sequence_across_windows_and_contigs() {
        let length = 4 * WINDOW_FLANK as usize;
        let first = (0..length).map(|i| b"ACGTacgtN"[i % 9]).collect::<Vec<_>>();
        let second = b"ggccaatt".to_vec();

        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &[("first", &first), ("second", &second)]);
        let mut fasta = IndexedFasta::open(&path).unwrap();

        let length = length as u64;
        let fetches = [
            ("first", 2 * WINDOW_FLANK, 2 * WINDOW_FLANK + 10),
            ("first", 2 * WINDOW_FLANK - 1, 2 * WINDOW_FLANK),
            ("first", 0, 5),
            ("second", 2, 6),
            ("first", length - 7, length),
            ("first", 3 * WINDOW_FLANK, length),
            ("first", 0, length),
        ];

        for (contig, start, end) in fetches {
            let bases = if contig == "first" { &first } else { &second };
            let expected = bases[start as usize..end as usize].to_ascii_uppercase();

            assert_eq!(
                fasta.fetch(contig, start, end).unwrap(),
                Some(expected),
                "{contig}:{start}-{end}"
            );
        }

        assert!(fasta.fetch("first", 0, length + 1).is_err());
        assert_eq!(fasta.fetch("missing", 0, 1).unwrap(), None);
    }
}
