//! Reads BED records and compares lifted BED segments.

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::BufRead as _;
use std::io::BufReader;
use std::path::Path;

use anyhow::Context as _;
use anyhow::Result;
use anyhow::bail;

use crate::report::Outcome;

/// A lifted BED segment.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Segment {
    /// The contig.
    pub contig: String,

    /// The 0-based start.
    pub start: u64,

    /// The exclusive end.
    pub end: u64,

    /// The strand (`+` or `-`).
    pub strand: String,
}

impl fmt::Display for Segment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}-{}({})",
            self.contig, self.start, self.end, self.strand
        )
    }
}

/// Formats a list of segments for a disagreement message.
pub fn describe(segments: &[Segment]) -> String {
    if segments.is_empty() {
        return String::from("unmapped");
    }

    segments
        .iter()
        .map(|segment| segment.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

/// Reads the segments of a BED6 file keyed by name.
///
/// Lines that are empty or that start with `#` are skipped.
pub fn read_segments(path: &Path) -> Result<HashMap<String, Vec<Segment>>> {
    let mut segments = HashMap::<String, Vec<Segment>>::new();

    for (fields, line_number) in read_records(path)? {
        let [contig, start, end, name, _, strand, ..] = fields.as_slice() else {
            bail!(
                "{}:{line_number}: expected at least six columns",
                path.display()
            );
        };

        let parse = |value: &str| {
            value
                .parse::<u64>()
                .with_context(|| format!("{}:{line_number}: invalid position", path.display()))
        };

        segments.entry(name.clone()).or_default().push(Segment {
            contig: contig.clone(),
            start: parse(start)?,
            end: parse(end)?,
            strand: strand.clone(),
        });
    }

    Ok(segments)
}

/// Reads the names of the records in a BED file in order.
pub fn read_names(path: &Path) -> Result<Vec<String>> {
    read_records(path)?
        .into_iter()
        .map(|(fields, line_number)| {
            fields.get(3).cloned().with_context(|| {
                format!("{}:{line_number}: expected a name column", path.display())
            })
        })
        .collect()
}

/// Reads the tab-separated fields of the records in a BED file along with
/// their 1-based line numbers.
fn read_records(path: &Path) -> Result<Vec<(Vec<String>, usize)>> {
    let file = File::open(path).with_context(|| format!("opening `{}`", path.display()))?;
    let mut records = Vec::new();

    for (i, line) in BufReader::new(file).lines().enumerate() {
        let line = line.with_context(|| format!("reading `{}`", path.display()))?;

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        records.push((line.split('\t').map(String::from).collect(), i + 1));
    }

    Ok(records)
}

/// Compares lifted segments that are expected to be identical (such as the
/// output of CrossMap, which also reports one segment per aligned block).
pub fn compare_exact(chainfile: &[Segment], other: &[Segment]) -> Outcome {
    match (chainfile.is_empty(), other.is_empty()) {
        (true, true) => Outcome::BothUnmapped,
        (false, true) => Outcome::OnlyChainfile,
        (true, false) => Outcome::OnlyOther,
        (false, false) => {
            let mut chainfile = chainfile.to_vec();
            let mut other = other.to_vec();
            chainfile.sort();
            other.sort();

            if chainfile == other {
                Outcome::Identical
            } else {
                Outcome::Differs
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segment(contig: &str, start: u64, end: u64, strand: &str) -> Segment {
        Segment {
            contig: contig.into(),
            start,
            end,
            strand: strand.into(),
        }
    }

    #[test]
    fn exact_comparison_ignores_order() {
        let a = [segment("chr1", 1, 2, "+"), segment("chr2", 3, 4, "-")];
        let b = [segment("chr2", 3, 4, "-"), segment("chr1", 1, 2, "+")];
        assert_eq!(compare_exact(&a, &b), Outcome::Identical);
        assert_eq!(compare_exact(&a, &b[..1]), Outcome::Differs);
        assert_eq!(compare_exact(&a, &[]), Outcome::OnlyChainfile);
        assert_eq!(compare_exact(&[], &b), Outcome::OnlyOther);
        assert_eq!(compare_exact(&[], &[]), Outcome::BothUnmapped);
    }
}
