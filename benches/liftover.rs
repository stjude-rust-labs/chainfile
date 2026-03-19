mod support;

use std::hint::black_box;

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
const NUM_QUERIES: usize = 1_000;

fn build_machine() -> chainfile::liftover::Machine {
    let data = generate_chain_data(NUM_CONTIGS, BLOCKS_PER_CONTIG, BLOCK_SIZE, GAP_SIZE);
    let reader = chainfile::Reader::new(&data[..]);
    chainfile::liftover::machine::Builder
        .try_build_from(reader)
        .unwrap()
}

fn total_contig_size() -> u64 {
    BLOCKS_PER_CONTIG as u64 * BLOCK_SIZE + (BLOCKS_PER_CONTIG as u64 - 1) * GAP_SIZE
}

fn generate_single_position_intervals(count: usize) -> Vec<Interval> {
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

fn generate_small_intervals(count: usize) -> Vec<Interval> {
    let mut rng = StdRng::seed_from_u64(SEED);

    (0..count)
        .map(|_| {
            let contig_idx = rng.gen_range(1..=NUM_CONTIGS);
            let contig = Contig::new_unchecked(format!("chr{contig_idx}"));
            let pos = rng.gen_range(0..BLOCK_SIZE - 100);

            let from = Coordinate::new(contig.clone(), Strand::Positive, pos);
            let to = Coordinate::new(contig, Strand::Positive, pos + 100);
            Interval::try_new(from, to).unwrap()
        })
        .collect()
}

fn bench_single_position_hit(c: &mut Criterion) {
    let machine = build_machine();
    let intervals = generate_single_position_intervals(NUM_QUERIES);

    c.bench_function("liftover/single_position_hit", |b| {
        b.iter(|| {
            for interval in &intervals {
                black_box(machine.liftover(interval.clone()));
            }
        });
    });
}

fn bench_single_position_miss(c: &mut Criterion) {
    let machine = build_machine();
    let contig = Contig::new_unchecked("nonexistent");
    let from = Coordinate::new(contig.clone(), Strand::Positive, 0u64);
    let to = Coordinate::new(contig, Strand::Positive, 1u64);
    let interval = Interval::try_new(from, to).unwrap();

    c.bench_function("liftover/single_position_miss", |b| {
        b.iter(|| black_box(machine.liftover(interval.clone())));
    });
}

fn bench_small_interval(c: &mut Criterion) {
    let machine = build_machine();
    let intervals = generate_small_intervals(NUM_QUERIES);

    c.bench_function("liftover/small_interval", |b| {
        b.iter(|| {
            for interval in &intervals {
                black_box(machine.liftover(interval.clone()));
            }
        });
    });
}

fn bench_large_interval_straddling_gaps(c: &mut Criterion) {
    let machine = build_machine();
    let total = total_contig_size();
    let contig = Contig::new_unchecked("chr1");
    let from = Coordinate::new(contig.clone(), Strand::Positive, 0u64);
    let to = Coordinate::new(contig, Strand::Positive, total);
    let interval = Interval::try_new(from, to).unwrap();

    c.bench_function("liftover/large_interval_straddling_gaps", |b| {
        b.iter(|| black_box(machine.liftover(interval.clone())));
    });
}

criterion_group!(
    benches,
    bench_single_position_hit,
    bench_single_position_miss,
    bench_small_interval,
    bench_large_interval_straddling_gaps
);
criterion_main!(benches);
