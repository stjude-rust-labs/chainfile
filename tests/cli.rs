//! Integration tests for the `chainfile` command line tool.

use std::io::Write as _;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::process::Output;
use std::process::Stdio;

use flate2::Compression;
use flate2::write::GzEncoder;

/// Gets the path to a test fixture.
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

/// Runs the `chainfile` binary with the given arguments and standard input.
fn run(args: &[&str], stdin: Option<&[u8]>) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chainfile"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("command should spawn");

    {
        let mut handle = child.stdin.take().expect("stdin should be piped");
        if let Some(stdin) = stdin {
            handle.write_all(stdin).expect("stdin should be written");
        }
    }

    child.wait_with_output().expect("command should complete")
}

/// Runs the `chainfile` binary and asserts that it succeeded.
fn run_ok(args: &[&str], stdin: Option<&[u8]>) -> String {
    let output = run(args, stdin);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("stdout should be UTF-8")
}

/// Gets the non-meta lines of a VCF.
fn vcf_body(vcf: &str) -> String {
    vcf.lines()
        .filter(|line| !line.starts_with("##"))
        .map(|line| format!("{line}\n"))
        .collect()
}

/// Builds tab-delimited output from rows of fields.
fn tsv(rows: &[&[&str]]) -> String {
    rows.iter()
        .map(|row| format!("{}\n", row.join("\t")))
        .collect()
}

fn chain() -> String {
    fixture("test.chain").display().to_string()
}

#[test]
fn liftover_bed() {
    let dir = tempfile::tempdir().unwrap();
    let unmapped = dir.path().join("unmapped.bed");
    let input = fixture("input.bed");

    let stdout = run_ok(
        &[
            "liftover",
            "bed",
            "--chain",
            &chain(),
            input.to_str().unwrap(),
            "--unmapped",
            unmapped.to_str().unwrap(),
        ],
        None,
    );

    assert_eq!(
        stdout,
        tsv(&[
            &["chrB", "5", "15", "span_gap", "0", "+"],
            &["chrB", "15", "23", "span_gap", "0", "+"],
            &["chrC", "5", "11", "flip", "0", "-"],
            &["chrB", "13", "14"],
        ])
    );

    assert_eq!(
        std::fs::read_to_string(unmapped).unwrap(),
        tsv(&[
            &["#contig not in chain file"],
            &["chrZ", "0", "5", "unknown"],
            &["#no alignment in query genome"],
            &["chrA", "10", "12", "gap"],
            &["#zero-length interval"],
            &["chrA", "5", "5", "zero"],
        ])
    );
}

#[test]
fn liftover_bed_from_stdin_with_gzipped_chain() {
    let dir = tempfile::tempdir().unwrap();
    let chain_gz = dir.path().join("test.chain.gz");

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(&std::fs::read(fixture("test.chain")).unwrap())
        .unwrap();
    std::fs::write(&chain_gz, encoder.finish().unwrap()).unwrap();

    let stdout = run_ok(
        &["liftover", "bed", "--chain", chain_gz.to_str().unwrap()],
        Some(b"chrA\t8\t9\n"),
    );

    assert_eq!(stdout, "chrB\t13\t14\n");
}

#[test]
fn liftover_bed_extra_columns() {
    let input = b"chrA\t0\t5\tname\t0\t+\t0\t5\n";
    let chain = chain();
    let args = ["liftover", "bed", "--chain", chain.as_str()];

    let output = run(&args, Some(input));
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--include-bed-plus"));
    assert!(stderr.contains("--ignore-bed-plus"));

    let stdout = run_ok(&[&args[..], &["--include-bed-plus"]].concat(), Some(input));
    assert_eq!(stdout, "chrB\t5\t10\tname\t0\t+\t0\t5\n");

    let stdout = run_ok(&[&args[..], &["--ignore-bed-plus"]].concat(), Some(input));
    assert_eq!(stdout, "chrB\t5\t10\tname\t0\t+\n");

    let output = run(
        &[&args[..], &["--include-bed-plus", "--ignore-bed-plus"]].concat(),
        Some(input),
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot be used with"));
}

#[test]
fn liftover_vcf() {
    let dir = tempfile::tempdir().unwrap();
    let rejected = dir.path().join("rejected.vcf");
    let input = fixture("input.vcf");
    let fasta = fixture("query.fa");

    let stdout = run_ok(
        &[
            "liftover",
            "vcf",
            "--chain",
            &chain(),
            "--query-fasta",
            fasta.to_str().unwrap(),
            input.to_str().unwrap(),
            "--rejected",
            rejected.to_str().unwrap(),
        ],
        None,
    );

    assert!(stdout.contains("##contig=<ID=chrB,length=30>\n"));
    assert!(stdout.contains("##contig=<ID=chrC,length=15>\n"));
    assert!(!stdout.contains("##contig=<ID=chrA"));
    assert!(stdout.contains("##INFO=<ID=LiftoverSwappedAlleles,Number=0,Type=Flag,"));

    assert_eq!(
        vcf_body(&stdout),
        tsv(&[
            &[
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"
            ],
            &[
                "chrB",
                "8",
                "snv_positive",
                "T",
                "G",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
            &[
                "chrB",
                "11",
                "snv_swapped",
                "C",
                "A",
                "50",
                "PASS",
                "AF=0.75;AC=2;AN=2;LiftoverSwappedAlleles",
                "GT:AD:PL",
                "1/1:1,9:300,30,0"
            ],
            &[
                "chrB",
                "17",
                "deletion_positive",
                "ACA",
                "A",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
            &[
                "chrB",
                "17",
                "delins_positive",
                "AC",
                "G",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
            &[
                "chrC",
                "12",
                "snv_negative",
                "A",
                "G",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
            &[
                "chrC",
                "7",
                "deletion_negative",
                "CG",
                "C",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
            &[
                "chrC",
                "7",
                "delins_negative",
                "CGT",
                "CC",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
            &[
                "chrC",
                "5",
                "insertion_negative",
                "G",
                "GCC",
                "50",
                "PASS",
                ".",
                "GT",
                "0/1"
            ],
        ])
    );

    let rejected = std::fs::read_to_string(rejected).unwrap();
    assert!(rejected.contains("##INFO=<ID=LiftoverRejectReason,"));
    assert_eq!(
        vcf_body(&rejected),
        tsv(&[
            &[
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"
            ],
            &[
                "chrA",
                "7",
                "snv_mismatch",
                "G",
                "T",
                "50",
                "PASS",
                "LiftoverRejectReason=reference_mismatch",
                "GT",
                "1/1"
            ],
            &[
                "chrA",
                "10",
                "spans_gap",
                "AATT",
                "A",
                "50",
                "PASS",
                "LiftoverRejectReason=partial_mapping",
                "GT",
                "0/1"
            ],
            &[
                "chrA",
                "6",
                "mnp_mismatch",
                "AG",
                "TT",
                "50",
                "PASS",
                "LiftoverRejectReason=reference_mismatch",
                "GT",
                "0/1"
            ],
            &[
                "chrZ",
                "2",
                "unknown_contig",
                "C",
                "A",
                "50",
                "PASS",
                "LiftoverRejectReason=unknown_contig",
                "GT",
                "0/1"
            ],
            &[
                "chrA",
                "2",
                "symbolic",
                "C",
                "<DEL>",
                "50",
                "PASS",
                "END=5;LiftoverRejectReason=symbolic_allele",
                "GT",
                "0/1"
            ],
        ])
    );
}

#[test]
fn view() {
    let stdout = run_ok(&["view", &chain()], None);
    assert_eq!(
        stdout,
        tsv(&[
            &[
                "chain_id",
                "reference_contig",
                "reference_strand",
                "reference_start",
                "reference_end",
                "query_contig",
                "query_strand",
                "query_start",
                "query_end",
                "size",
                "dt",
                "dq"
            ],
            &[
                "1", "chrA", "+", "0", "10", "chrB", "+", "5", "15", "10", "2", "0"
            ],
            &[
                "1", "chrA", "+", "12", "20", "chrB", "+", "15", "23", "8", ".", "."
            ],
            &[
                "2", "chrA", "+", "20", "30", "chrC", "-", "3", "13", "10", ".", "."
            ],
        ])
    );

    let stdout = run_ok(&["view", "--headers", "--id", "2", &chain()], None);
    assert_eq!(stdout, "chain 500 chrA 30 + 20 30 chrC 15 - 2 12 2\n");
}

#[test]
fn stats() {
    let stdout = run_ok(&["stats", &chain()], None);
    assert_eq!(
        stdout,
        tsv(&[
            &["name", "value"],
            &["chains", "2"],
            &["positive_strand_chains", "1"],
            &["negative_strand_chains", "1"],
            &["blocks", "3"],
            &["aligned_bases", "28"],
            &["reference_gap_bases", "2"],
            &["query_gap_bases", "0"],
            &["reference_contigs", "1"],
            &["query_contigs", "2"],
        ])
    );
}

#[test]
fn check() {
    let reference = fixture("reference.fa");
    let query = fixture("query.fa");
    let chain = chain();
    let args = [
        "check",
        "--chain",
        &chain,
        "--reference-fasta",
        reference.to_str().unwrap(),
        "--query-fasta",
        query.to_str().unwrap(),
    ];

    let stdout = run_ok(&args, None);
    assert_eq!(
        stdout,
        tsv(&[
            &[
                "query_strand",
                "aligned_bases",
                "compared_bases",
                "mismatches",
                "identity"
            ],
            &["all", "28", "28", "2", "0.928571"],
            &["+", "18", "18", "2", "0.888889"],
            &["-", "10", "10", "0", "1.000000"],
        ])
    );

    let output = run(&[&args[..], &["--min-identity", "0.9"]].concat(), None);
    assert!(output.status.success());

    let output = run(&[&args[..], &["--min-identity", "0.99"]].concat(), None);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("below the minimum"));
}
