use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A square symmetric distance matrix in scikit-bio's `DistanceMatrix` TSV form:
/// an empty top-left cell then the IDs as the header row, then one row per ID
/// (ID + delimited distances). This is what `rsomics-beta-diversity` emits.
///
/// Parsing enforces the same invariants skbio's `DistanceMatrix` constructor
/// does — no NaN cells, symmetric, hollow (zero) diagonal — so a matrix that
/// skbio would reject is rejected here too rather than yielding a meaningless R.
pub struct DistanceMatrix {
    pub ids: Vec<String>,
    /// Row-major `n × n`.
    pub data: Vec<f64>,
}

impl DistanceMatrix {
    /// # Errors
    /// Missing header, ragged / non-square body, a row count that disagrees with
    /// the header, a mislabelled row, a non-numeric cell, or a matrix that
    /// violates skbio's invariants (NaN cell, asymmetry, non-zero diagonal).
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<DistanceMatrix> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty distance matrix".into())),
            }
        };
        let ids: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        let n = ids.len();
        if n == 0 {
            return Err(RsomicsError::InvalidInput(
                "header has no ID columns (need an empty top-left cell + ≥1 ID)".into(),
            ));
        }

        let mut data = vec![0.0_f64; n * n];
        let mut row = 0usize;
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            if row >= n {
                return Err(RsomicsError::InvalidInput(format!(
                    "more data rows than the {n} IDs in the header"
                )));
            }
            let mut fields = line.split(delim);
            let label = fields.next().unwrap_or("").trim();
            if label != ids[row] {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} label '{label}' does not match header ID '{}'",
                    row + 1,
                    ids[row]
                )));
            }
            let mut col = 0usize;
            for field in fields {
                if col >= n {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} ('{label}') has more values than the {n} IDs",
                        row + 1
                    )));
                }
                let v: f64 = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row {} ('{label}'), column {}: '{}' is not numeric",
                        row + 1,
                        col + 1,
                        field.trim()
                    ))
                })?;
                if v.is_nan() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} ('{label}'), column {}: NaN cell — distance matrix must contain no NaNs",
                        row + 1,
                        col + 1
                    )));
                }
                data[row * n + col] = v;
                col += 1;
            }
            if col != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} ('{label}') has {col} values, expected {n}",
                    row + 1
                )));
            }
            row += 1;
        }
        if row != n {
            return Err(RsomicsError::InvalidInput(format!(
                "{row} data rows but {n} IDs in the header"
            )));
        }
        for i in 0..n {
            if data[i * n + i] != 0.0 {
                return Err(RsomicsError::InvalidInput(format!(
                    "diagonal must be zero: [{}][{}] = {} (matrix must be hollow)",
                    i + 1,
                    i + 1,
                    data[i * n + i]
                )));
            }
            for j in (i + 1)..n {
                let (a, b) = (data[i * n + j], data[j * n + i]);
                if a != b {
                    return Err(RsomicsError::InvalidInput(format!(
                        "matrix is not symmetric: [{}][{}] = {a} but [{}][{}] = {b}",
                        i + 1,
                        j + 1,
                        j + 1,
                        i + 1
                    )));
                }
            }
        }
        Ok(DistanceMatrix { ids, data })
    }

    #[must_use]
    pub fn n(&self) -> usize {
        self.ids.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Result<DistanceMatrix> {
        DistanceMatrix::parse(s.as_bytes(), '\t')
    }

    fn err_of(s: &str) -> String {
        match parse(s) {
            Ok(_) => panic!("expected parse to reject: {s:?}"),
            Err(e) => e.to_string(),
        }
    }

    #[test]
    fn nan_cell_rejected() {
        let err = err_of("\ta\tb\na\t0\tnan\nb\tnan\t0\n");
        assert!(err.contains("NaN"), "{err}");
    }

    #[test]
    fn asymmetric_rejected() {
        let err = err_of("\ta\tb\na\t0\t1\nb\t2\t0\n");
        assert!(err.contains("not symmetric"), "{err}");
    }

    #[test]
    fn nonzero_diagonal_rejected() {
        let err = err_of("\ta\tb\na\t9\t1\nb\t1\t0\n");
        assert!(err.contains("hollow"), "{err}");
    }

    #[test]
    fn symmetric_hollow_accepted() {
        let dm = parse("\ta\tb\na\t0\t1\nb\t1\t0\n").unwrap();
        assert_eq!(dm.n(), 2);
        assert_eq!(dm.data, vec![0.0, 1.0, 1.0, 0.0]);
    }
}
