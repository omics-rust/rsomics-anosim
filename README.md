# rsomics-anosim

ANOSIM (Analysis of Similarities, Clarke 1993) test for significant differences
between groups of samples, computed from a square distance matrix. A Rust port
of `skbio.stats.distance.anosim`.

```
rsomics-anosim dm.tsv --grouping groups.tsv [--permutations 999] [--seed S] [-o result.tsv]
```

- `dm.tsv` — square distance-matrix TSV (the `DistanceMatrix` form
  `rsomics-beta-diversity` emits: empty top-left cell, IDs as the header row,
  then one `id<TAB>distances` row per sample). Reads stdin when `-` or omitted.
- `--grouping groups.tsv` — one `id<TAB>group` line per sample. Extra IDs are
  ignored; a matrix ID missing from the file is an error.
- `--permutations` — permutations for the p-value (default 999; `0` skips it and
  reports `nan`).
- `--seed` — permutation RNG seed (default 0).

Output is the skbio result table: method name, test statistic name (R), sample
size, number of groups, the R statistic, number of permutations, and the
p-value.

## Method

The R statistic is

```
R = (r_B − r_W) / (N·(N−1)/4)
```

where `r_W` / `r_B` are the mean ranks of within- / between-group pairwise
distances, ranked over all pairwise distances with ties resolved by average
rank. R ranges from −1 to +1; 0 means random grouping. Significance is assessed
by permuting the group labels and recomputing R; the p-value is the fraction of
permuted R values ≥ the observed R, with the standard `(count + 1)/(perm + 1)`
correction.

## Origin

This crate ports `scikit-bio`'s `skbio.stats.distance.anosim`
(<https://github.com/scikit-bio/scikit-bio>, BSD-3-Clause). scikit-bio's BSD
license permits reading and citing its source, which informed the exact ranking
(`scipy.stats.rankdata(method="average")` over the condensed upper triangle),
the `N·(N−1)/4` divisor, and the result-table layout.

The **R statistic is value-exact** against scikit-bio (agreement to ~1e-9,
including the average-rank tie handling). The **permutation p-value** uses an
independent SplitMix64-seeded Fisher-Yates shuffle (one stream per permutation,
so the p-value is deterministic for a given `--seed` regardless of thread
count), *not* numpy's PCG64 `Generator.permutation`. It is therefore a
seeded-permutation estimate that converges to scikit-bio's value as the
permutation count grows, rather than a bit-for-bit reproduction of numpy's
draws.

Primary reference: Clarke, K.R. (1993). "Non-parametric multivariate analyses of
changes in community structure." *Australian Journal of Ecology* 18(1): 117-143.
doi:10.1111/j.1442-9993.1993.tb00438.x

License: MIT OR Apache-2.0.
