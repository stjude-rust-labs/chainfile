mod support;

use std::hint::black_box;

use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use support::generate_chain_data;

fn bench_machine_building(c: &mut Criterion) {
    let mut group = c.benchmark_group("machine_building");

    let configs = [
        ("small", 5, 100, 500, 50),
        ("medium", 25, 1_000, 1_000, 100),
        ("large", 25, 10_000, 1_000, 100),
    ];

    for (label, num_contigs, blocks, block_size, gap_size) in configs {
        let data = generate_chain_data(num_contigs, blocks, block_size, gap_size);

        group.bench_with_input(BenchmarkId::from_parameter(label), &data, |b, data| {
            b.iter(|| {
                let reader = chainfile::Reader::new(&data[..]);
                black_box(
                    chainfile::liftover::machine::Builder
                        .try_build_from(reader)
                        .unwrap(),
                )
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_machine_building);
criterion_main!(benches);
