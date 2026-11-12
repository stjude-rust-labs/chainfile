//! End-to-end tests of the generated data.

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use chainfile::liftover::LiftoverResult;
use chainfile::liftover::Machine;
use chainfile::liftover::machine::builder::Builder;
use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;
use rand::SeedableRng as _;
use rand::rngs::StdRng;
use tempfile::TempDir;

use crate::Situation;
use crate::bed;
use crate::fasta::IndexedFasta;
use crate::fixture;
use crate::fixture::Fixture;
use crate::index::Index;
use crate::sequence::forward_bounds;
use crate::sequence::oriented;
use crate::vcf;

/// The percentage of records generated for each situation.
const PERCENTAGE: f64 = 8.0;

/// The number of records generated.
const COUNT: usize = 1_000;

/// The maximum indel length.
const MAX_INDEL_LENGTH: u64 = 10;

/// Gets the options for a small fixture.
fn small_options() -> fixture::Options {
    fixture::Options {
        contigs: 2,
        contig_length: 200_000,
        ..Default::default()
    }
}

/// Generates a small fixture.
fn small_fixture(seed: u64) -> (TempDir, Fixture) {
    small_fixture_with(seed, &small_options())
}

/// Generates a small fixture with the given options.
fn small_fixture_with(seed: u64, options: &fixture::Options) -> (TempDir, Fixture) {
    let directory = tempfile::tempdir().unwrap();

    let fixture = fixture::generate(
        options,
        &mut StdRng::seed_from_u64(seed),
        &directory.path().join("fixture"),
    )
    .unwrap();

    (directory, fixture)
}

/// Gets the BED options with every situation at [`PERCENTAGE`].
fn bed_options() -> bed::Options {
    bed::Options {
        count: COUNT,
        min_length: 1,
        max_length: 1_000,
        situations: crate::SharedSituations {
            spans_gap: PERCENTAGE,
            unaligned: PERCENTAGE,
            multiple_chains: PERCENTAGE,
            opposite_strand: PERCENTAGE,
        },
    }
}

/// Gets the VCF options with every situation at [`PERCENTAGE`].
fn vcf_options() -> vcf::Options {
    vcf::Options {
        count: COUNT,
        max_indel_length: MAX_INDEL_LENGTH,
        situations: bed_options().situations,
        vcf_situations: crate::VcfSituations {
            indels: PERCENTAGE,
            swapped_alleles: PERCENTAGE,
            reference_mismatch: PERCENTAGE,
            multiallelic: PERCENTAGE,
        },
    }
}

/// Generates BED records.
fn generate_bed(chain: &Path, options: &bed::Options, seed: u64) -> anyhow::Result<String> {
    let mut output = Vec::new();
    bed::generate(
        chain,
        options,
        &mut StdRng::seed_from_u64(seed),
        &mut output,
    )?;
    Ok(String::from_utf8(output).unwrap())
}

/// Generates VCF records.
fn generate_vcf(fixture: &Fixture, options: &vcf::Options, seed: u64) -> anyhow::Result<String> {
    let mut output = Vec::new();
    vcf::generate(
        &fixture.chain,
        &fixture.reference_fasta,
        &fixture.query_fasta,
        options,
        &mut StdRng::seed_from_u64(seed),
        &mut output,
    )?;
    Ok(String::from_utf8(output).unwrap())
}

/// Builds a liftover machine from a chain file.
fn machine(chain: &Path) -> Machine {
    let reader = chainfile::Reader::new(BufReader::new(File::open(chain).unwrap()));
    Builder.try_build_from(reader).unwrap()
}

/// Lifts over the interval `[start, end)` on the positive strand of a
/// reference contig.
fn liftover(machine: &Machine, contig: &str, start: u64, end: u64) -> Option<Vec<LiftoverResult>> {
    let interval = Interval::try_new(
        Coordinate::new(Contig::new_unchecked(contig), Strand::Positive, start),
        Coordinate::new(Contig::new_unchecked(contig), Strand::Positive, end),
    )
    .unwrap();

    machine.liftover(interval)
}

/// Gets the query strand and forward-strand query bounds when an interval
/// maps to exactly one segment of one chain covering the whole interval.
fn single(
    results: &Option<Vec<LiftoverResult>>,
    start: u64,
    end: u64,
) -> Option<(String, Strand, u64, u64)> {
    let [result] = results.as_deref()? else {
        return None;
    };

    let [segment] = result.segments() else {
        return None;
    };

    if forward_bounds(segment.reference()) != (start, end) {
        return None;
    }

    let (query_start, query_end) = forward_bounds(segment.query());

    Some((
        segment.query().contig().as_str().to_string(),
        segment.query().strand(),
        query_start,
        query_end,
    ))
}

/// Counts the records for each situation.
fn situations<'a>(names: impl Iterator<Item = &'a str>) -> BTreeMap<Situation, usize> {
    let mut counts = BTreeMap::new();

    for name in names {
        let situation = Situation::from_record_name(name)
            .unwrap_or_else(|| panic!("`{name}` is not a generated record name"));
        *counts.entry(situation).or_default() += 1;
    }

    counts
}

#[test]
fn fixture_is_consistent() {
    let options = fixture::Options {
        duplication_pct: 50.0,
        ..small_options()
    };
    let (_directory, fixture) = small_fixture_with(1, &options);

    let machine = machine(&fixture.chain);
    let mut reference = IndexedFasta::open(&fixture.reference_fasta).unwrap();
    let mut query = IndexedFasta::open(&fixture.query_fasta).unwrap();

    let data = std::fs::read(&fixture.chain).unwrap();
    let mut reader = chainfile::Reader::new(&data[..]);

    let mut chains = 0;
    let mut strands = BTreeSet::new();
    let mut gaps = BTreeSet::new();
    let mut other_contig = false;
    let mut compared = 0u64;
    let mut matches = 0u64;

    for result in reader.sections() {
        let section = result.unwrap();
        let header = section.header();
        let reference_sequence = header.reference_sequence();
        let query_sequence = header.query_sequence();

        chains += 1;
        assert_eq!(header.id(), chains);
        strands.insert(query_sequence.strand());
        other_contig |= reference_sequence.chromosome_name() != query_sequence.chromosome_name();

        let mut reference_length = 0;
        let mut query_length = 0;

        for result in section.stepthrough_with_data().unwrap() {
            let (pair, data) = result.unwrap();

            if let (Some(dt), Some(dq)) = (data.dt(), data.dq()) {
                assert!(dt > 0 || dq > 0, "chain {} has an empty gap", header.id());
                gaps.insert((dt > 0, dq > 0));
                reference_length += dt;
                query_length += dq;
            }

            reference_length += data.size();
            query_length += data.size();

            let (start, end) = forward_bounds(pair.reference());
            let (query_start, query_end) = forward_bounds(pair.query());
            assert_eq!(end - start, data.size());
            assert_eq!(query_end - query_start, data.size());

            let reference_bases = reference
                .fetch(reference_sequence.chromosome_name(), start, end)
                .unwrap();
            let query_bases = oriented(
                query
                    .fetch(query_sequence.chromosome_name(), query_start, query_end)
                    .unwrap(),
                pair.query().strand(),
            );

            let block_matches = reference_bases
                .iter()
                .zip(&query_bases)
                .filter(|(r, q)| r == q)
                .count() as u64;

            assert!(
                block_matches as f64 / data.size() as f64 >= 0.95,
                "block {start}-{end} of chain {} has low identity",
                header.id()
            );

            compared += data.size();
            matches += block_matches;

            let results = liftover(&machine, reference_sequence.chromosome_name(), start, end);
            assert!(
                results
                    .unwrap()
                    .iter()
                    .any(|result| result.chain().id() == header.id())
            );
        }

        assert_eq!(
            reference_length,
            reference_sequence.alignment_end() - reference_sequence.alignment_start()
        );
        assert_eq!(
            query_length,
            query_sequence.alignment_end() - query_sequence.alignment_start()
        );
    }

    assert!(chains > 10);
    assert_eq!(strands.len(), 2);
    assert!(other_contig);
    assert_eq!(
        gaps,
        BTreeSet::from([(true, false), (false, true), (true, true)])
    );

    let identity = matches as f64 / compared as f64;
    assert!((0.98..1.0).contains(&identity), "identity is {identity}");

    let index = Index::load(&fixture.chain).unwrap();
    assert!(!index.overlaps().is_empty());
    assert!(!index.unaligned().is_empty());
    assert!(!index.gaps().is_empty());
}

#[test]
fn bed_records_have_their_situations() {
    let (_directory, fixture) = small_fixture(2);
    let output = generate_bed(&fixture.chain, &bed_options(), 3).unwrap();
    let machine = machine(&fixture.chain);

    let records = output
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 6);
            assert_eq!(&fields[4..], ["0", "+"]);
            (
                fields[0].to_string(),
                fields[1].parse::<u64>().unwrap(),
                fields[2].parse::<u64>().unwrap(),
                fields[3].to_string(),
            )
        })
        .collect::<Vec<_>>();

    let mut sorted = records.clone();
    sorted.sort();
    assert_eq!(records, sorted);

    let counts = situations(records.iter().map(|(.., name)| name.as_str()));
    assert_eq!(
        counts,
        BTreeMap::from([
            (Situation::Contained, 680),
            (Situation::SpansGap, 80),
            (Situation::Unaligned, 80),
            (Situation::MultipleChains, 80),
            (Situation::OppositeStrand, 80),
        ])
    );

    for (contig, start, end, name) in &records {
        assert!((1..=1_000).contains(&(end - start)));

        let results = liftover(&machine, contig, *start, *end);

        let ok = match Situation::from_record_name(name).unwrap() {
            Situation::Contained => single(&results, *start, *end)
                .is_some_and(|(_, strand, ..)| strand == Strand::Positive),
            Situation::OppositeStrand => single(&results, *start, *end)
                .is_some_and(|(_, strand, ..)| strand == Strand::Negative),
            Situation::SpansGap => {
                matches!(results.as_deref(), Some([result]) if result.segments().len() >= 2)
            }
            Situation::Unaligned => results.is_none(),
            Situation::MultipleChains => results.is_some_and(|results| results.len() >= 2),
            situation => panic!("unexpected situation: {situation}"),
        };

        assert!(ok, "{name} at {contig}:{start}-{end} has the wrong outcome");
    }
}

#[test]
fn vcf_records_have_their_situations() {
    let (_directory, fixture) = small_fixture(4);
    let output = generate_vcf(&fixture, &vcf_options(), 5).unwrap();
    let machine = machine(&fixture.chain);
    let mut reference = IndexedFasta::open(&fixture.reference_fasta).unwrap();
    let mut query = IndexedFasta::open(&fixture.query_fasta).unwrap();

    let (header, body): (Vec<_>, Vec<_>) = output.lines().partition(|line| line.starts_with('#'));

    assert_eq!(
        header,
        [
            "##fileformat=VCFv4.3",
            "##contig=<ID=chr1,length=200000>",
            "##contig=<ID=chr2,length=200000>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            &[
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample",
            ]
            .join("\t"),
        ]
    );

    let counts = situations(body.iter().map(|line| line.split('\t').nth(2).unwrap()));
    let mut expected = BTreeMap::from(Situation::ALL.map(|situation| (situation, 80)));
    expected.insert(Situation::Contained, 360);
    assert_eq!(counts, expected);

    let mut previous = None;

    for line in body {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), 10);
        assert_eq!(&fields[5..9], [".", "PASS", ".", "GT"]);

        let contig = fields[0];
        let start = fields[1].parse::<u64>().unwrap() - 1;
        let name = fields[2];
        let reference_allele = fields[3].as_bytes();
        let alternates = fields[4]
            .split(',')
            .map(|allele| allele.as_bytes())
            .collect::<Vec<_>>();
        let genotype = fields[9];
        let end = start + reference_allele.len() as u64;

        let key = (contig.to_string(), start, name.to_string());
        assert!(previous.as_ref().is_none_or(|previous| *previous <= key));
        previous = Some(key);

        assert!(
            fields[3..5]
                .iter()
                .flat_map(|allele| allele.bytes())
                .all(|base| b"ACGT,".contains(&base))
        );

        let reference_bases = reference.fetch(contig, start, end).unwrap();
        let results = liftover(&machine, contig, start, end);

        let mut oriented_query = || {
            let (query_contig, strand, query_start, query_end) =
                single(&results, start, end).expect("record should map to a single segment");
            let bases = query.fetch(&query_contig, query_start, query_end).unwrap();
            (oriented(bases, strand), strand, query_start)
        };

        let situation = Situation::from_record_name(name).unwrap();

        if situation == Situation::Multiallelic {
            assert_eq!(genotype, "1/2");
        } else {
            assert!(["0/1", "1/1"].contains(&genotype));
            assert_eq!(alternates.len(), 1, "{name} should be biallelic");
        }

        match situation {
            Situation::Contained | Situation::OppositeStrand => {
                let (query_bases, strand, _) = oriented_query();
                let expected = if situation == Situation::Contained {
                    Strand::Positive
                } else {
                    Strand::Negative
                };
                assert_eq!(strand, expected);
                assert_eq!(reference_allele.len(), 1);
                assert_eq!(reference_allele, reference_bases);
                assert_eq!(reference_allele, query_bases);
                assert_eq!(alternates[0].len(), 1);
                assert_ne!(alternates[0], reference_allele);
            }
            Situation::Indel => {
                let (query_bases, strand, query_start) = oriented_query();
                assert_eq!(reference_allele, reference_bases);
                assert_eq!(reference_allele, query_bases);
                assert_eq!(alternates[0][0], reference_allele[0]);

                let (short, long) = if reference_allele.len() == 1 {
                    (reference_allele, alternates[0])
                } else {
                    (alternates[0], reference_allele)
                };
                assert_eq!(short.len(), 1);
                assert!((2..=MAX_INDEL_LENGTH as usize + 1).contains(&long.len()));

                if strand == Strand::Negative {
                    assert!(query_start >= 1);
                }
            }
            Situation::SwappedAlleles => {
                let (query_bases, ..) = oriented_query();
                assert_eq!(reference_allele, reference_bases);
                assert_eq!(alternates[0], query_bases);
                assert_ne!(reference_allele, query_bases);
            }
            Situation::ReferenceMismatch => {
                let (query_bases, ..) = oriented_query();
                assert_eq!(reference_allele.len(), 1);
                assert_ne!(reference_allele, reference_bases);
                assert_ne!(reference_allele, query_bases);
                assert_ne!(alternates[0], query_bases);
                assert_ne!(alternates[0], reference_allele);
            }
            Situation::Multiallelic => {
                let (query_bases, ..) = oriented_query();
                assert_eq!(reference_allele, reference_bases);
                assert_eq!(reference_allele, query_bases);
                assert_eq!(alternates.len(), 2);
                assert_ne!(alternates[0], alternates[1]);
                assert!(
                    alternates
                        .iter()
                        .all(|alternate| alternate.len() == 1 && *alternate != reference_allele)
                );
            }
            Situation::SpansGap => {
                assert!(
                    matches!(results.as_deref(), Some([result]) if result.segments().len() >= 2)
                );
                assert_eq!(reference_allele, reference_bases);
                assert!(reference_allele.len() <= MAX_INDEL_LENGTH as usize + 1);
                assert_eq!(alternates[0], &reference_allele[..1]);
            }
            Situation::Unaligned => {
                assert!(results.is_none());
                assert_eq!(reference_allele, reference_bases);
            }
            Situation::MultipleChains => {
                assert!(results.is_some_and(|results| results.len() >= 2));
                assert_eq!(reference_allele, reference_bases);

                if reference_allele.len() > 1 {
                    assert!(reference_allele.len() <= MAX_INDEL_LENGTH as usize + 1);
                    assert_eq!(alternates[0], &reference_allele[..1]);
                } else {
                    assert_ne!(alternates[0], reference_allele);
                }
            }
        }
    }
}

#[test]
fn generation_is_deterministic() {
    let read = |fixture: &Fixture| -> Vec<Vec<u8>> {
        let mut paths = vec![
            fixture.chain.clone(),
            fixture.reference_fasta.clone(),
            fixture.query_fasta.clone(),
        ];

        for path in [&fixture.reference_fasta, &fixture.query_fasta] {
            let mut index = path.as_os_str().to_owned();
            index.push(".fai");
            paths.push(PathBuf::from(index));
        }

        paths
            .iter()
            .map(|path| std::fs::read(path).unwrap())
            .collect()
    };

    let (_first_directory, first) = small_fixture(6);
    let (_second_directory, second) = small_fixture(6);
    let (_third_directory, third) = small_fixture(7);

    assert_eq!(read(&first), read(&second));
    assert_ne!(read(&first), read(&third));

    let bed_options = bed::Options {
        count: 200,
        ..bed_options()
    };

    let bed = generate_bed(&first.chain, &bed_options, 8).unwrap();
    assert_eq!(bed, generate_bed(&second.chain, &bed_options, 8).unwrap());
    assert_ne!(bed, generate_bed(&first.chain, &bed_options, 9).unwrap());

    let vcf_options = vcf::Options {
        count: 200,
        ..vcf_options()
    };

    let vcf = generate_vcf(&first, &vcf_options, 8).unwrap();
    assert_eq!(vcf, generate_vcf(&second, &vcf_options, 8).unwrap());
    assert_ne!(vcf, generate_vcf(&first, &vcf_options, 9).unwrap());
}

#[test]
fn invalid_options_are_rejected() {
    let (_directory, fixture) = small_fixture(10);

    let mut options = bed_options();
    options.situations.spans_gap = 60.0;
    options.situations.unaligned = 50.0;
    let error = generate_bed(&fixture.chain, &options, 0).unwrap_err();
    assert!(error.to_string().contains("more than 100"), "{error}");

    let mut options = vcf_options();
    options.situations.spans_gap = 50.0;
    options.vcf_situations.indels = 50.0;
    let error = generate_vcf(&fixture, &options, 0).unwrap_err();
    assert!(error.to_string().contains("more than 100"), "{error}");

    let options = bed::Options {
        min_length: 0,
        ..bed_options()
    };
    assert!(generate_bed(&fixture.chain, &options, 0).is_err());

    let options = bed::Options {
        min_length: 10,
        max_length: 9,
        ..bed_options()
    };
    assert!(generate_bed(&fixture.chain, &options, 0).is_err());

    let options = vcf::Options {
        max_indel_length: 0,
        ..vcf_options()
    };
    assert!(generate_vcf(&fixture, &options, 0).is_err());

    let options = vcf::Options {
        max_indel_length: u64::MAX,
        ..vcf_options()
    };
    assert!(generate_vcf(&fixture, &options, 0).is_err());

    let options = fixture::Options {
        contigs: 0,
        ..Default::default()
    };
    let directory = tempfile::tempdir().unwrap();
    assert!(fixture::generate(&options, &mut StdRng::seed_from_u64(0), directory.path()).is_err());
}

#[test]
fn bed_records_are_generated_for_the_test_fixture() {
    let chain = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/test.chain");

    let options = bed::Options {
        count: 100,
        min_length: 1,
        max_length: 5,
        situations: crate::SharedSituations {
            spans_gap: 20.0,
            unaligned: 20.0,
            multiple_chains: 20.0,
            opposite_strand: 20.0,
        },
    };

    let output = generate_bed(&chain, &options, 11).unwrap();
    let mut counts = BTreeMap::<Situation, usize>::new();

    for line in output.lines() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[0], "chrA");

        let start = fields[1].parse::<u64>().unwrap();
        let end = fields[2].parse::<u64>().unwrap();
        let situation = Situation::from_record_name(fields[3]).unwrap();
        *counts.entry(situation).or_default() += 1;

        // NOTE: the chain file aligns `chrA:0-10` and `chrA:12-20` to `chrB`
        // on the positive strand (with a gap at `chrA:10-12`) and `chrA:20-30`
        // to `chrC` on the negative strand (with a junction between the chains
        // at `chrA:20`).
        let ok = match situation {
            Situation::Contained => (end <= 10) || (12 <= start && end <= 20),
            Situation::SpansGap => start <= 9 && (13..=20).contains(&end),
            Situation::Unaligned => 10 <= start && end <= 12,
            Situation::OppositeStrand => 20 <= start && end <= 30,
            Situation::MultipleChains => (12..=19).contains(&start) && 21 <= end,
            situation => panic!("unexpected situation: {situation}"),
        };

        assert!(ok, "{line}");
    }

    assert_eq!(
        counts,
        BTreeMap::from([
            (Situation::Contained, 20),
            (Situation::SpansGap, 20),
            (Situation::Unaligned, 20),
            (Situation::MultipleChains, 20),
            (Situation::OppositeStrand, 20),
        ])
    );

    let directory = tempfile::tempdir().unwrap();
    let compressed = directory.path().join("test.chain.gz");
    let mut encoder =
        flate2::write::GzEncoder::new(File::create(&compressed).unwrap(), Default::default());
    encoder.write_all(&std::fs::read(&chain).unwrap()).unwrap();
    encoder.finish().unwrap();
    assert_eq!(generate_bed(&compressed, &options, 11).unwrap(), output);

    let options = bed::Options {
        situations: crate::SharedSituations {
            multiple_chains: 10.0,
            ..Default::default()
        },
        ..options
    };
    let single = directory.path().join("single.chain");
    std::fs::write(
        &single,
        ["chain 10 chrA 100 + 10 20 chrB 100 + 0 10 1", "10", ""].join("\n"),
    )
    .unwrap();
    let error = generate_bed(&single, &options, 11).unwrap_err();
    assert!(
        format!("{error:#}").contains("the chain file has no multiple_chains candidates"),
        "{error:#}"
    );
}

#[test]
fn vcf_records_span_the_junction_between_chains_in_the_test_fixture() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures");
    let chain = fixtures.join("test.chain");
    let reference_fasta = fixtures.join("reference.fa");
    let query_fasta = fixtures.join("query.fa");

    let options = vcf::Options {
        count: 20,
        max_indel_length: 3,
        situations: crate::SharedSituations {
            multiple_chains: 100.0,
            ..Default::default()
        },
        vcf_situations: Default::default(),
    };

    let mut output = Vec::new();
    vcf::generate(
        &chain,
        &reference_fasta,
        &query_fasta,
        &options,
        &mut StdRng::seed_from_u64(0),
        &mut output,
    )
    .unwrap();

    let output = String::from_utf8(output).unwrap();
    let mut reference = IndexedFasta::open(&reference_fasta).unwrap();
    let records = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(records.len(), 20);

    for line in records {
        let fields = line.split('\t').collect::<Vec<_>>();
        let start = fields[1].parse::<u64>().unwrap() - 1;
        let end = start + fields[3].len() as u64;

        // NOTE: chain 1 ends at `chrA:20` and chain 2 starts there, so each
        // record is a padded deletion that covers bases of both chains.
        assert_eq!(fields[0], "chrA");
        assert!(
            (12..=19).contains(&start) && (21..=start + 4).contains(&end),
            "{line}"
        );
        assert_eq!(
            fields[3].as_bytes(),
            reference.fetch("chrA", start, end).unwrap()
        );
        assert_eq!(fields[4], &fields[3][..1]);
    }
}

#[test]
fn bed_records_span_gaps_near_the_start_of_a_contig() {
    // NOTE: this is a regression test for sampling intervals that span a gap
    // whose window is much longer than the position of the block after it.
    let directory = tempfile::tempdir().unwrap();
    let chain = directory.path().join("start.chain");
    std::fs::write(
        &chain,
        [
            "chain 100 chrA 1000 + 0 500 chrB 1000 + 0 499 1",
            "2\t1\t0",
            "497",
            "",
        ]
        .join("\n"),
    )
    .unwrap();

    let options = bed::Options {
        count: 200,
        min_length: 1,
        max_length: 1_000,
        situations: crate::SharedSituations {
            spans_gap: 100.0,
            ..Default::default()
        },
    };

    let output = generate_bed(&chain, &options, 0).unwrap();
    assert_eq!(output.lines().count(), 200);

    for line in output.lines() {
        let fields = line.split('\t').collect::<Vec<_>>();
        let start = fields[1].parse::<u64>().unwrap();
        let end = fields[2].parse::<u64>().unwrap();
        assert!(start <= 1 && (4..=500).contains(&end), "{line}");
    }

    let options = bed::Options {
        max_length: 2,
        ..options
    };
    let error = generate_bed(&chain, &options, 0).unwrap_err();
    assert!(
        format!("{error:#}").contains("the chain file has no spans_gap candidates"),
        "{error:#}"
    );
}
