use std::io::{BufRead, Write};

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

mod fmt;
mod rng;

use fmt::push_pyrepr;
use rng::SplitMix64;
pub use rsomics_distance::DistanceMatrix;

#[derive(Serialize)]
pub struct AnosimResult {
    pub sample_size: usize,
    pub num_groups: usize,
    pub r: f64,
    pub permutations: usize,
    pub p_value: f64,
}

/// ANOSIM (Clarke 1993). The R statistic is value-exact against
/// `skbio.stats.distance.anosim`; the p-value is a seeded permutation estimate
/// (its own RNG, not numpy's, so it differs in the last digits from skbio for a
/// given seed but converges to the same value as `permutations` grows).
///
/// # Errors
/// All grouping labels identical (no between distances) or all distinct (no
/// within distances) — both forbidden by the test, matching skbio.
pub fn anosim(
    dm: &DistanceMatrix,
    grouping: &[String],
    permutations: usize,
    seed: u64,
) -> Result<AnosimResult> {
    let n = dm.n();
    if grouping.len() != n {
        return Err(RsomicsError::InvalidInput(format!(
            "grouping has {} entries but the matrix has {n} samples",
            grouping.len()
        )));
    }
    let codes = factorize(grouping);
    let num_groups = codes.iter().copied().max().map_or(0, |m| m + 1);
    if num_groups == 1 {
        return Err(RsomicsError::InvalidInput(
            "all samples are in a single group — no between-group distances".into(),
        ));
    }
    if num_groups == n {
        return Err(RsomicsError::InvalidInput(
            "every sample is in its own group — no within-group distances".into(),
        ));
    }

    let ranked = ranked_condensed(&dm.data, n);
    let divisor = n as f64 * ((n as f64 - 1.0) / 4.0);
    let r = r_stat(&ranked, &codes, n, divisor);

    let p_value = if permutations == 0 {
        f64::NAN
    } else {
        let ge = (0..permutations)
            .into_par_iter()
            .filter(|&i| {
                let mut rng =
                    SplitMix64::new(seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
                let mut perm = codes.clone();
                shuffle(&mut perm, &mut rng);
                r_stat(&ranked, &perm, n, divisor) >= r
            })
            .count();
        (ge + 1) as f64 / (permutations as f64 + 1.0)
    };

    Ok(AnosimResult {
        sample_size: n,
        num_groups,
        r,
        permutations,
        p_value,
    })
}

/// Map labels to dense integer codes. The codes only carry equality (R depends
/// solely on whether a pair shares a code), so any stable assignment matches
/// skbio's `np.unique(..., return_inverse=True)` for R purposes.
fn factorize(grouping: &[String]) -> Vec<usize> {
    let mut order: Vec<&String> = grouping.iter().collect();
    order.sort_unstable();
    order.dedup();
    grouping
        .iter()
        .map(|g| order.binary_search(&g).unwrap())
        .collect()
}

/// Tie-averaged ranks (scipy `rankdata(method="average")`, 1-based) of the
/// upper-triangle distances in triu row-major order (i<j).
fn ranked_condensed(data: &[f64], n: usize) -> Vec<f64> {
    let m = n * (n - 1) / 2;
    let mut vals = Vec::with_capacity(m);
    for i in 0..n {
        for j in (i + 1)..n {
            vals.push(data[i * n + j]);
        }
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| vals[a].total_cmp(&vals[b]));

    let mut ranks = vec![0.0_f64; m];
    let mut k = 0usize;
    while k < m {
        let mut j = k + 1;
        while j < m && vals[order[j]] == vals[order[k]] {
            j += 1;
        }
        // ranks k..j (0-based) are 1-based k+1..j; average = (k+1 + j) / 2
        let avg = (k + 1 + j) as f64 / 2.0;
        for &idx in &order[k..j] {
            ranks[idx] = avg;
        }
        k = j;
    }
    ranks
}

/// R = (mean between-rank − mean within-rank) / divisor, over the triu pairs.
fn r_stat(ranked: &[f64], codes: &[usize], n: usize, divisor: f64) -> f64 {
    let mut sum_w = 0.0_f64;
    let mut cnt_w = 0usize;
    let mut sum_b = 0.0_f64;
    let mut cnt_b = 0usize;
    let mut p = 0usize;
    for i in 0..n {
        let ci = codes[i];
        for &cj in &codes[i + 1..n] {
            let rk = ranked[p];
            if ci == cj {
                sum_w += rk;
                cnt_w += 1;
            } else {
                sum_b += rk;
                cnt_b += 1;
            }
            p += 1;
        }
    }
    let r_w = sum_w / cnt_w as f64;
    let r_b = sum_b / cnt_b as f64;
    (r_b - r_w) / divisor
}

fn shuffle(v: &mut [usize], rng: &mut SplitMix64) {
    for i in (1..v.len()).rev() {
        let j = rng.bounded(i as u64 + 1) as usize;
        v.swap(i, j);
    }
}

impl AnosimResult {
    /// Write the skbio-style result table (one `field<TAB>value` line each).
    /// Floats use Python `repr` so R is byte-comparable against the oracle.
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W) -> Result<()> {
        let mut line = String::new();
        let row = |out: &mut W, k: &str, v: &str| -> Result<()> {
            writeln!(out, "{k}\t{v}").map_err(RsomicsError::Io)
        };
        row(&mut out, "method name", "ANOSIM")?;
        row(&mut out, "test statistic name", "R")?;
        row(&mut out, "sample size", &self.sample_size.to_string())?;
        row(&mut out, "number of groups", &self.num_groups.to_string())?;
        line.clear();
        push_pyrepr(&mut line, self.r);
        row(&mut out, "test statistic", &line)?;
        line.clear();
        push_pyrepr(&mut line, self.p_value);
        row(&mut out, "p-value", &line)?;
        row(
            &mut out,
            "number of permutations",
            &self.permutations.to_string(),
        )
    }
}

/// Parse the grouping TSV: one `id<TAB>group` line per sample (`#` comment /
/// header lines skipped), then return the group label for each matrix ID in
/// matrix order. Extra IDs in the grouping file are ignored; an ID missing from
/// it is an error (matching skbio's superset rule).
///
/// # Errors
/// A matrix ID absent from the grouping file, or a malformed line.
pub fn read_grouping<R: BufRead>(reader: R, ids: &[String]) -> Result<Vec<String>> {
    use std::collections::HashMap;
    let mut map: HashMap<String, String> = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut it = line.splitn(2, '\t');
        let id = it.next().unwrap().trim();
        let grp = it
            .next()
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!("grouping line lacks a group column: '{line}'"))
            })?
            .trim();
        map.insert(id.to_string(), grp.to_string());
    }
    ids.iter()
        .map(|id| {
            map.get(id).cloned().ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "sample '{id}' has no entry in the grouping file"
                ))
            })
        })
        .collect()
}

/// # Errors
/// Propagates parse, grouping, and write errors.
pub fn run<W: Write>(
    dm_reader: impl BufRead,
    grouping_reader: impl BufRead,
    out: W,
    delim: char,
    permutations: usize,
    seed: u64,
) -> Result<()> {
    let dm = DistanceMatrix::parse(dm_reader, delim)?;
    let grouping = read_grouping(grouping_reader, &dm.ids)?;
    let res = anosim(&dm, &grouping, permutations, seed)?;
    res.write_tsv(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc_dm() -> DistanceMatrix {
        DistanceMatrix::parse(
            "\ts1\ts2\ts3\ts4\n\
             s1\t0\t1\t1\t4\n\
             s2\t1\t0\t3\t2\n\
             s3\t1\t3\t0\t3\n\
             s4\t4\t2\t3\t0\n"
                .as_bytes(),
            '\t',
        )
        .unwrap()
    }

    #[test]
    fn doc_example_r_is_quarter() {
        let g: Vec<String> = ["Group1", "Group1", "Group2", "Group2"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let res = anosim(&doc_dm(), &g, 0, 42).unwrap();
        assert!((res.r - 0.25).abs() < 1e-12, "R = {}", res.r);
        assert_eq!(res.num_groups, 2);
        assert_eq!(res.sample_size, 4);
        assert!(res.p_value.is_nan());
    }

    #[test]
    fn tie_averaged_ranks() {
        // ranks of [10, 20, 20, 30] are [1, 2.5, 2.5, 4]
        let r = ranked_condensed(&[0.0, 10.0, 20.0, 10.0, 0.0, 20.0, 20.0, 20.0, 0.0], 3);
        // condensed order: (0,1)=10, (0,2)=20, (1,2)=20  -> values [10,20,20]
        assert_eq!(r, vec![1.0, 2.5, 2.5]);
    }

    #[test]
    fn perfect_separation_r_one() {
        // three tight clusters, all within < all between -> R = 1
        let d = "\ta\tb\tc\td\te\tf\n\
                 a\t0\t1\t1\t5\t5\t5\n\
                 b\t1\t0\t1\t5\t5\t5\n\
                 c\t1\t1\t0\t5\t5\t5\n\
                 d\t5\t5\t5\t0\t1\t1\n\
                 e\t5\t5\t5\t1\t0\t1\n\
                 f\t5\t5\t5\t1\t1\t0\n";
        let dm = DistanceMatrix::parse(d.as_bytes(), '\t').unwrap();
        let g: Vec<String> = ["A", "A", "A", "B", "B", "B"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let res = anosim(&dm, &g, 0, 1).unwrap();
        assert!((res.r - 1.0).abs() < 1e-12, "R = {}", res.r);
    }

    #[test]
    fn single_group_errors() {
        let g: Vec<String> = vec!["X".into(); 4];
        assert!(anosim(&doc_dm(), &g, 0, 1).is_err());
    }

    #[test]
    fn all_distinct_errors() {
        let g: Vec<String> = ["a", "b", "c", "d"].iter().map(|s| s.to_string()).collect();
        assert!(anosim(&doc_dm(), &g, 0, 1).is_err());
    }

    #[test]
    fn p_value_deterministic_across_thread_counts() {
        let g: Vec<String> = ["Group1", "Group1", "Group2", "Group2"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let a = anosim(&doc_dm(), &g, 99, 7).unwrap().p_value;
        let b = anosim(&doc_dm(), &g, 99, 7).unwrap().p_value;
        assert_eq!(a, b);
    }
}
