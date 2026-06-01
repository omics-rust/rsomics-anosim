use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_anosim::{DistanceMatrix, anosim};

fn synth(n: usize, ngroups: usize) -> (DistanceMatrix, Vec<String>) {
    let mut data = vec![0.0_f64; n * n];
    let g: Vec<String> = (0..n).map(|i| format!("G{}", i % ngroups)).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            let same = i % ngroups == j % ngroups;
            let base = if same { 1.0 } else { 5.0 };
            let d = base + ((i * 31 + j * 17) % 7) as f64 * 0.1;
            data[i * n + j] = d;
            data[j * n + i] = d;
        }
    }
    let ids = (0..n).map(|i| format!("s{i}")).collect();
    (DistanceMatrix { ids, data }, g)
}

fn bench(c: &mut Criterion) {
    let (dm, g) = synth(800, 5);
    c.bench_function("anosim_800_999perm", |b| {
        b.iter(|| anosim(black_box(&dm), black_box(&g), 999, 42).unwrap());
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
