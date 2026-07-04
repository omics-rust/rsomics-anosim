use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn ours_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-anosim"))
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn oracle_script() -> String {
    format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"))
}

fn parse_result(text: &str) -> HashMap<String, String> {
    text.lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()
}

fn ours(table: &str, perms: usize, seed: u64) -> HashMap<String, String> {
    let out = Command::new(ours_bin())
        .arg(golden(&format!("{table}_dm.tsv")))
        .args(["--grouping", &golden(&format!("{table}_groups.tsv"))])
        .args(["--permutations", &perms.to_string()])
        .args(["--seed", &seed.to_string()])
        .output()
        .expect("run rsomics-anosim");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_result(&String::from_utf8(out.stdout).unwrap())
}

/// Run ours against a distance-matrix golden with a shared grouping, returning
/// the raw process output so degenerate-input tests can assert a clean failure.
fn ours_raw(dm: &str, groups: &str) -> std::process::Output {
    Command::new(ours_bin())
        .arg(golden(dm))
        .args(["--grouping", &golden(groups)])
        .args(["--permutations", "0"])
        .output()
        .expect("run rsomics-anosim")
}

/// skbio's `DistanceMatrix` constructor rejects NaN cells and asymmetric input
/// (both raise `DistanceMatrixError`, exit 1). Ours must reject them the same
/// way — non-zero exit, a clear stderr message, and crucially no panic (a panic
/// would print `thread 'main' panicked` and abort, not fail cleanly). Committed
/// goldens, no live oracle; skbio parity is confirmed in `skbio_rejects_degenerate`.
#[test]
fn rejects_degenerate_matrices() {
    for dm in ["nan_dm.tsv", "asym_dm.tsv"] {
        let out = ours_raw(dm, "degenerate_groups.tsv");
        assert!(
            !out.status.success(),
            "{dm}: ours accepted a matrix skbio rejects (exit 0)"
        );
        assert!(
            out.status.code().is_some(),
            "{dm}: ours did not exit cleanly (killed by signal / abort)"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            !stderr.contains("panicked"),
            "{dm}: ours panicked instead of failing loudly:\n{stderr}"
        );
        assert!(
            !stderr.trim().is_empty(),
            "{dm}: ours failed without an error message"
        );
    }
}

/// Committed skbio-captured R values (permutations=0, deterministic). Runs with
/// no scikit-bio present — this is the always-on regression gate.
#[test]
fn matches_committed_golden_r() {
    let expected = std::fs::read_to_string(golden("expected_R.tsv")).unwrap();
    for line in expected.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (name, ss, ng, r) = (f[0], f[1], f[2], f[3].parse::<f64>().unwrap());
        let got = ours(name, 0, 42);
        assert_eq!(got["sample size"], ss, "{name} sample size");
        assert_eq!(got["number of groups"], ng, "{name} num groups");
        let got_r: f64 = got["test statistic"].parse().unwrap();
        assert!(
            (got_r - r).abs() <= 1e-9,
            "{name} R: ours {got_r} vs skbio {r}"
        );
        assert_eq!(got["method name"], "ANOSIM");
        assert_eq!(got["test statistic name"], "R");
        assert_eq!(got["p-value"], "nan", "{name} p must be nan at 0 perms");
    }
}

/// scikit-bio is the named oracle; skip loudly if it (or python) is unavailable.
/// `RSOMICS_SKBIO_PYTHON` overrides the interpreter (e.g. an isolated venv).
fn skbio_python() -> Option<String> {
    let mut candidates = Vec::new();
    if let Ok(p) = std::env::var("RSOMICS_SKBIO_PYTHON") {
        candidates.push(p);
    }
    candidates.push("python3".into());
    candidates.push("python".into());
    for py in candidates {
        let probe = Command::new(&py)
            .args(["-c", "import skbio.stats.distance"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py);
        }
    }
    eprintln!("SKIP: scikit-bio not importable — install `scikit-bio` to run the differential");
    None
}

fn oracle(py: &str, table: &str, perms: usize, seed: u64) -> HashMap<String, String> {
    let out = Command::new(py)
        .arg(oracle_script())
        .arg(golden(&format!("{table}_dm.tsv")))
        .arg(golden(&format!("{table}_groups.tsv")))
        .arg(perms.to_string())
        .arg(seed.to_string())
        .output()
        .expect("run scikit-bio oracle");
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_result(&String::from_utf8(out.stdout).unwrap())
}

/// Live differential: R value-exact (1e-9) vs skbio across all goldens. The
/// permutation p-value uses a different RNG than numpy's, so it is checked for
/// plausibility (same direction, both in (0, 1]) rather than bit-equality —
/// documented in README ## Origin.
fn differential(table: &str) {
    let Some(py) = skbio_python() else { return };
    let o = ours(table, 999, 7);
    let t = oracle(&py, table, 999, 7);

    assert_eq!(o["sample size"], t["sample size"], "{table} N");
    assert_eq!(
        o["number of groups"], t["number of groups"],
        "{table} groups"
    );

    let or: f64 = o["test statistic"].parse().unwrap();
    let tr: f64 = t["test statistic"].parse().unwrap();
    assert!(
        (or - tr).abs() <= 1e-9,
        "{table} R: ours {or} vs skbio {tr}"
    );

    let op: f64 = o["p-value"].parse().unwrap();
    let tp: f64 = t["p-value"].parse().unwrap();
    assert!(op > 0.0 && op <= 1.0, "{table} ours p out of range: {op}");
    assert!(tp > 0.0 && tp <= 1.0, "{table} skbio p out of range: {tp}");
    // both should land on the same side of 0.05 for these well-separated goldens
    assert_eq!(
        op < 0.05,
        tp < 0.05,
        "{table} p significance disagrees: ours {op} vs skbio {tp}"
    );
}

/// Live parity for the rejection path: feed the same degenerate goldens to
/// skbio and confirm its `DistanceMatrix` constructor raises (exit 1). Gated on
/// scikit-bio being importable, like the numeric differentials.
#[test]
fn skbio_rejects_degenerate() {
    let Some(py) = skbio_python() else { return };
    for dm in ["nan_dm.tsv", "asym_dm.tsv"] {
        let out = Command::new(&py)
            .arg(oracle_script())
            .arg(golden(dm))
            .arg(golden("degenerate_groups.tsv"))
            .arg("0")
            .arg("42")
            .output()
            .expect("run scikit-bio oracle");
        assert!(
            !out.status.success(),
            "{dm}: skbio accepted a matrix we now reject"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("DistanceMatrixError"),
            "{dm}: skbio failed for an unexpected reason:\n{stderr}"
        );
    }
}

#[test]
fn matches_skbio_doc() {
    differential("doc");
}

#[test]
fn matches_skbio_tie() {
    differential("tie");
}

#[test]
fn matches_skbio_med() {
    differential("med");
}
