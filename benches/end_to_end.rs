mod support;

use std::hint::black_box;

use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use omics::coordinate::Contig;
use omics::coordinate::Strand;
use omics::coordinate::interbase::Coordinate;
use omics::coordinate::interval::interbase::Interval;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use support::generate_chain_data;

const NUM_CONTIGS: usize = 25;
const BLOCKS_PER_CONTIG: usize = 1_000;
const BLOCK_SIZE: u64 = 1_000;
const GAP_SIZE: u64 = 100;
const SEED: u64 = 42;

fn total_contig_size() -> u64 {
    BLOCKS_PER_CONTIG as u64 * BLOCK_SIZE + (BLOCKS_PER_CONTIG as u64 - 1) * GAP_SIZE
}

fn generate_query_intervals(count: usize) -> Vec<Interval> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let total = total_contig_size();

    (0..count)
        .map(|_| {
            let contig_idx = rng.gen_range(1..=NUM_CONTIGS);
            let contig = Contig::new_unchecked(format!("chr{contig_idx}"));
            let pos = rng.gen_range(0..total - 1);

            let from = Coordinate::new(contig.clone(), Strand::Positive, pos);
            let to = Coordinate::new(contig, Strand::Positive, pos + 1);
            Interval::try_new(from, to).unwrap()
        })
        .collect()
}

fn bench_query_only(c: &mut Criterion) {
    let data = generate_chain_data(NUM_CONTIGS, BLOCKS_PER_CONTIG, BLOCK_SIZE, GAP_SIZE);
    let reader = chainfile::Reader::new(&data[..]);
    let machine = chainfile::liftover::machine::Builder
        .try_build_from(reader)
        .unwrap();

    let mut group = c.benchmark_group("query_only");

    for num_queries in [1_000, 10_000, 100_000] {
        let intervals = generate_query_intervals(num_queries);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{num_queries}_queries")),
            &intervals,
            |b, intervals| {
                b.iter(|| {
                    for interval in intervals {
                        black_box(machine.liftover(interval.clone()));
                    }
                });
            },
        );
    }

    group.finish();
}

fn bench_end_to_end(c: &mut Criterion) {
    let data = generate_chain_data(NUM_CONTIGS, BLOCKS_PER_CONTIG, BLOCK_SIZE, GAP_SIZE);

    let mut group = c.benchmark_group("end_to_end");

    for num_queries in [1_000, 10_000, 100_000] {
        let intervals = generate_query_intervals(num_queries);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{num_queries}_queries")),
            &intervals,
            |b, intervals| {
                b.iter(|| {
                    let reader = chainfile::Reader::new(&data[..]);
                    let machine = chainfile::liftover::machine::Builder
                        .try_build_from(reader)
                        .unwrap();

                    for interval in intervals {
                        black_box(machine.liftover(interval.clone()));
                    }
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_query_only, bench_end_to_end);
criterion_main!(benches);
