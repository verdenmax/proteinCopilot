//! Loader for generic TSV/TXT search result files into [`UnifiedPsm`] records.
//!
//! Uses a configurable [`TsvColumnMap`] to map column headers to the fields
//! of [`UnifiedPsm`]. The only essential column is `peptide`; all others are
//! optional and default to `None` when the column is missing or the value
//! cannot be parsed.

use std::path::Path;

use serde::{Deserialize, Serialize};
use tracing;

use crate::error::EntrapmentError;
use crate::types::UnifiedPsm;

/// Column name mapping for generic TSV/TXT result files.
///
/// Each field specifies the header string in the file that corresponds to the
/// given [`UnifiedPsm`] field. Defaults are lowercase underscore names.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TsvColumnMap {
    /// Header for the peptide sequence column.
    #[serde(default = "default_peptide")]
    pub peptide: String,
    /// Header for the charge state column.
    #[serde(default = "default_charge")]
    pub charge: String,
    /// Header for the precursor m/z column.
    #[serde(default = "default_precursor_mz")]
    pub precursor_mz: String,
    /// Header for the retention time column.
    #[serde(default = "default_retention_time")]
    pub retention_time: String,
    /// Header for the scan number column.
    #[serde(default = "default_scan_number")]
    pub scan_number: String,
    /// Header for the spectrum file column.
    #[serde(default = "default_spectrum_file")]
    pub spectrum_file: String,
    /// Header for the protein IDs column.
    #[serde(default = "default_protein_ids")]
    pub protein_ids: String,
    /// Header for the q-value column.
    #[serde(default = "default_q_value")]
    pub q_value: String,
    /// Field delimiter character (default: tab).
    #[serde(default = "default_delimiter")]
    pub delimiter: char,
}

fn default_peptide() -> String {
    "peptide".to_string()
}
fn default_charge() -> String {
    "charge".to_string()
}
fn default_precursor_mz() -> String {
    "precursor_mz".to_string()
}
fn default_retention_time() -> String {
    "retention_time".to_string()
}
fn default_scan_number() -> String {
    "scan_number".to_string()
}
fn default_spectrum_file() -> String {
    "spectrum_file".to_string()
}
fn default_protein_ids() -> String {
    "protein_ids".to_string()
}
fn default_q_value() -> String {
    "q_value".to_string()
}
fn default_delimiter() -> char {
    '\t'
}

impl Default for TsvColumnMap {
    fn default() -> Self {
        Self {
            peptide: default_peptide(),
            charge: default_charge(),
            precursor_mz: default_precursor_mz(),
            retention_time: default_retention_time(),
            scan_number: default_scan_number(),
            spectrum_file: default_spectrum_file(),
            protein_ids: default_protein_ids(),
            q_value: default_q_value(),
            delimiter: default_delimiter(),
        }
    }
}

/// Load a generic TSV/TXT result file into a vector of [`UnifiedPsm`].
///
/// The `column_map` controls which headers map to which PSM fields.
/// The `peptide` column is essential; all other columns are optional.
pub fn load_generic_tsv(
    path: &Path,
    column_map: &TsvColumnMap,
) -> Result<Vec<UnifiedPsm>, EntrapmentError> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(column_map.delimiter as u8)
        .from_path(path)
        .map_err(|e| EntrapmentError::IoError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;

    let headers = rdr.headers().map_err(|e| EntrapmentError::LoaderError {
        format: "generic TSV".to_string(),
        detail: format!("failed to read headers from '{}': {e}", path.display()),
    })?;

    // Find column indices
    let peptide_idx =
        find_column(headers, &column_map.peptide).ok_or_else(|| EntrapmentError::LoaderError {
            format: "generic TSV".to_string(),
            detail: format!(
                "missing essential column '{}' in '{}'",
                column_map.peptide,
                path.display()
            ),
        })?;

    let charge_idx = find_column(headers, &column_map.charge);
    let precursor_mz_idx = find_column(headers, &column_map.precursor_mz);
    let retention_time_idx = find_column(headers, &column_map.retention_time);
    let scan_number_idx = find_column(headers, &column_map.scan_number);
    let spectrum_file_idx = find_column(headers, &column_map.spectrum_file);
    let protein_ids_idx = find_column(headers, &column_map.protein_ids);
    let q_value_idx = find_column(headers, &column_map.q_value);

    // Log unfound optional columns to aid debugging when data appears missing
    for (name, found) in [
        (&column_map.charge, charge_idx.is_some()),
        (&column_map.precursor_mz, precursor_mz_idx.is_some()),
        (&column_map.retention_time, retention_time_idx.is_some()),
        (&column_map.scan_number, scan_number_idx.is_some()),
        (&column_map.spectrum_file, spectrum_file_idx.is_some()),
        (&column_map.protein_ids, protein_ids_idx.is_some()),
        (&column_map.q_value, q_value_idx.is_some()),
    ] {
        if !found {
            tracing::debug!(
                column = %name,
                file = %path.display(),
                "optional TSV column not found in headers"
            );
        }
    }

    let mut psms = Vec::new();

    for result in rdr.records() {
        let record = result.map_err(|e| EntrapmentError::LoaderError {
            format: "generic TSV".to_string(),
            detail: format!("failed to parse row in '{}': {e}", path.display()),
        })?;

        let peptide = record.get(peptide_idx).unwrap_or("").trim().to_string();
        if peptide.is_empty() {
            continue;
        }

        let charge = charge_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<i32>().ok());

        let precursor_mz = precursor_mz_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        let retention_time = retention_time_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        let scan_number = scan_number_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<u32>().ok());

        let spectrum_file = spectrum_file_idx
            .and_then(|i| record.get(i))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());

        let protein_ids = protein_ids_idx
            .and_then(|i| record.get(i))
            .map(|s| s.trim().to_string())
            .unwrap_or_default();

        let q_value = q_value_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        psms.push(UnifiedPsm {
            peptide,
            charge,
            precursor_mz,
            retention_time,
            rt_start: None,
            rt_stop: None,
            scan_number,
            spectrum_file,
            protein_ids,
            q_value,
            modifications: Vec::new(),
        });
    }

    tracing::info!(
        "loaded {} PSMs from generic TSV: {}",
        psms.len(),
        path.display()
    );
    Ok(psms)
}

/// Find the index of a column header (case-sensitive exact match).
fn find_column(headers: &csv::StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|h| h == name)
}
