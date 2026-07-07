use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A square symmetric distance matrix in scikit-bio's `DistanceMatrix` TSV form:
/// an empty top-left cell then the IDs as the header row, then one row per ID
/// (ID + delimited distances). This is what `rsomics-beta-diversity` emits.
///
/// Parsing enforces the same invariants skbio's `DistanceMatrix` constructor
/// does — unique IDs, hollow (zero) diagonal, symmetric, no NaN cells — so a
/// matrix skbio would reject is rejected here too rather than yielding a
/// meaningless R.
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
        let dm = DistanceMatrix { ids, data };
        dm.validate()?;
        Ok(dm)
    }

    /// Enforce the invariants skbio's `DistanceMatrix` constructor checks, in its
    /// order — unique IDs, hollow diagonal, then symmetry — with skbio's exact
    /// wording. The symmetry pass also rejects any NaN cell: skbio compares cells
    /// with `!=`, and `NaN != NaN`, so a NaN never equals its mirror.
    ///
    /// # Errors
    /// Duplicate IDs, a non-zero diagonal cell, or an asymmetric / NaN cell.
    fn validate(&self) -> Result<()> {
        let n = self.ids.len();

        let mut seen = std::collections::HashSet::with_capacity(n);
        let mut dups = Vec::new();
        for id in &self.ids {
            if !seen.insert(id.as_str()) && !dups.contains(&id.as_str()) {
                dups.push(id.as_str());
            }
        }
        if !dups.is_empty() {
            let listed = dups
                .iter()
                .map(|d| format!("'{d}'"))
                .collect::<Vec<_>>()
                .join(", ");
            return Err(RsomicsError::InvalidInput(format!(
                "IDs must be unique. Found the following duplicate IDs: {listed}"
            )));
        }

        for i in 0..n {
            if self.data[i * n + i] != 0.0 {
                return Err(RsomicsError::InvalidInput(
                    "Data must be hollow (i.e., the diagonal can only contain zeros).".into(),
                ));
            }
        }

        for i in 0..n {
            for j in (i + 1)..n {
                if self.data[i * n + j] != self.data[j * n + i] {
                    return Err(RsomicsError::InvalidInput(
                        "Data must be symmetric and cannot contain NaNs.".into(),
                    ));
                }
            }
        }

        Ok(())
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
        assert!(err.contains("symmetric"), "{err}");
    }

    #[test]
    fn nonzero_diagonal_rejected() {
        let err = err_of("\ta\tb\na\t9\t1\nb\t1\t0\n");
        assert!(err.contains("hollow"), "{err}");
    }

    #[test]
    fn duplicate_ids_rejected() {
        let err = err_of("\ta\ta\na\t0\t1\na\t1\t0\n");
        assert!(err.contains("unique") && err.contains("'a'"), "{err}");
    }

    #[test]
    fn symmetric_hollow_accepted() {
        let dm = parse("\ta\tb\na\t0\t1\nb\t1\t0\n").unwrap();
        assert_eq!(dm.n(), 2);
        assert_eq!(dm.data, vec![0.0, 1.0, 1.0, 0.0]);
    }
}
