//! Tallies and reports comparison outcomes.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use anyhow::Result;
use anyhow::bail;
use chainfile_generate::Situation;

/// The outcome of comparing how two tools lifted over a record.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Outcome {
    /// Both tools produced the same coordinates.
    Identical,

    /// Neither tool mapped the record.
    BothUnmapped,

    /// Only `chainfile` mapped the record.
    OnlyChainfile,

    /// Only the other tool mapped the record.
    OnlyOther,

    /// The tools mapped the record to different coordinates.
    Differs,
}

impl Outcome {
    /// Gets the label of the outcome.
    pub fn label(&self, tool: &str) -> String {
        match self {
            Outcome::Identical => String::from("identical"),
            Outcome::BothUnmapped => String::from("both_unmapped"),
            Outcome::OnlyChainfile => String::from("only_chainfile"),
            Outcome::OnlyOther => format!("only_{tool}"),
            Outcome::Differs => String::from("differs"),
        }
    }

    /// Returns whether the outcome is a disagreement.
    pub fn is_failure(&self) -> bool {
        matches!(
            self,
            Outcome::OnlyChainfile | Outcome::OnlyOther | Outcome::Differs
        )
    }
}

/// A tally of comparison outcomes by situation.
pub struct Report {
    /// The name of the tool `chainfile` is compared with.
    tool: &'static str,

    /// The number of records by situation and outcome.
    counts: BTreeMap<String, BTreeMap<Outcome, u64>>,

    /// Descriptions of the first disagreements.
    examples: Vec<String>,

    /// The maximum number of disagreements to describe.
    max_examples: usize,

    /// The number of disagreements.
    failures: u64,
}

impl Report {
    /// Creates a new report.
    pub fn new(tool: &'static str, max_examples: usize) -> Self {
        Self {
            tool,
            counts: BTreeMap::new(),
            examples: Vec::new(),
            max_examples,
            failures: 0,
        }
    }

    /// Records the outcome for a record.
    ///
    /// The `describe` closure is only called for disagreements that are shown
    /// as examples.
    pub fn record(&mut self, name: &str, outcome: Outcome, describe: impl FnOnce() -> String) {
        let situation = Situation::from_record_name(name)
            .map(|situation| situation.name().to_string())
            .unwrap_or_else(|| String::from("unknown"));

        if outcome.is_failure() {
            self.failures += 1;

            if self.examples.len() < self.max_examples {
                self.examples.push(format!(
                    "{name} ({}): {}",
                    outcome.label(self.tool),
                    describe()
                ));
            }
        }

        *self
            .counts
            .entry(situation)
            .or_default()
            .entry(outcome)
            .or_default() += 1;
    }

    /// Renders the report as an aligned table followed by any examples.
    pub fn render(&self) -> String {
        let mut outcomes = self
            .counts
            .values()
            .flat_map(|counts| counts.keys().cloned())
            .collect::<Vec<_>>();
        outcomes.sort();
        outcomes.dedup();

        let mut header = vec![String::from("situation"), String::from("records")];
        header.extend(outcomes.iter().map(|outcome| outcome.label(self.tool)));

        let row = |situation: &str, counts: &BTreeMap<Outcome, u64>| {
            let mut row = vec![
                situation.to_string(),
                counts.values().sum::<u64>().to_string(),
            ];
            row.extend(
                outcomes
                    .iter()
                    .map(|outcome| counts.get(outcome).copied().unwrap_or(0).to_string()),
            );
            row
        };

        let mut total = BTreeMap::<Outcome, u64>::new();

        for counts in self.counts.values() {
            for (outcome, count) in counts {
                *total.entry(outcome.clone()).or_default() += count;
            }
        }

        let mut rows = vec![header];
        rows.extend(
            self.counts
                .iter()
                .map(|(situation, counts)| row(situation, counts)),
        );
        rows.push(row("total", &total));

        let widths = (0..rows[0].len())
            .map(|i| rows.iter().map(|row| row[i].len()).max().unwrap_or(0))
            .collect::<Vec<_>>();

        let mut rendered = String::new();

        for row in rows {
            let line = row
                .iter()
                .zip(&widths)
                .enumerate()
                .map(|(i, (cell, width))| {
                    if i == 0 {
                        format!("{cell:<width$}")
                    } else {
                        format!("{cell:>width$}")
                    }
                })
                .collect::<Vec<_>>()
                .join("  ");
            writeln!(rendered, "{}", line.trim_end()).unwrap();
        }

        if !self.examples.is_empty() {
            writeln!(rendered).unwrap();
            writeln!(rendered, "disagreements:").unwrap();

            for example in &self.examples {
                writeln!(rendered, "  {example}").unwrap();
            }
        }

        rendered
    }

    /// Prints the report and fails if there were any disagreements.
    pub fn finish(self, title: &str) -> Result<()> {
        println!("{title}");
        println!();
        print!("{}", self.render());

        if self.failures > 0 {
            bail!(
                "{} records disagree between `chainfile` and {}",
                self.failures,
                self.tool
            );
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn report_tallies_by_situation_and_fails_on_disagreements() {
        let mut report = Report::new("crossmap", 1);
        report.record("contained_0", Outcome::Identical, || unreachable!());
        report.record("spans_gap_0", Outcome::BothUnmapped, || unreachable!());
        report.record("spans_gap_1", Outcome::Differs, || String::from("first"));
        report.record("other", Outcome::OnlyOther, || String::from("second"));

        let rendered = report.render();
        assert!(rendered.contains("situation"));
        assert!(rendered.contains("only_crossmap"));
        assert!(rendered.contains("spans_gap_1 (differs): first"));
        assert!(!rendered.contains("second"));
        assert!(rendered.lines().any(|line| line.starts_with("unknown")));
        assert!(rendered.lines().any(|line| line.starts_with("total ")));

        assert!(report.finish("title").is_err());
    }

    #[test]
    fn report_without_disagreements_succeeds() {
        let mut report = Report::new("crossmap", 5);
        report.record("contained_0", Outcome::Identical, || unreachable!());
        report.record("unaligned_0", Outcome::BothUnmapped, || unreachable!());
        assert!(report.finish("title").is_ok());
    }
}
