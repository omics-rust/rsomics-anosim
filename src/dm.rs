use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A square symmetric distance matrix in scikit-bio's `DistanceMatrix` TSV form:
/// an empty top-left cell then the IDs as the header row, then one row per ID
/// (ID + delimited distances). This is what `rsomics-beta-diversity` emits.
pub struct DistanceMatrix {
    pub ids: Vec<String>,
    /// Row-major `n × n`.
    pub data: Vec<f64>,
}

impl DistanceMatrix {
    /// # Errors
    /// Missing header, ragged / non-square body, a row count that disagrees with
    /// the header, a mislabelled row, or a non-numeric cell.
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
                data[row * n + col] = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row {} ('{label}'), column {}: '{}' is not numeric",
                        row + 1,
                        col + 1,
                        field.trim()
                    ))
                })?;
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
        Ok(DistanceMatrix { ids, data })
    }

    #[must_use]
    pub fn n(&self) -> usize {
        self.ids.len()
    }
}
