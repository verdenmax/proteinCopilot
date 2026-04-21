//! Loader for DIA-NN `report.parquet` files into [`UnifiedPsm`] records.
//!
//! Reads the following columns:
//! - `Stripped.Sequence` — stripped peptide sequence (essential)
//! - `Protein.Ids` — semicolon-separated protein accessions (essential)
//! - `Modified.Sequence` — modified sequence with UniMod annotations (optional)
//! - `Precursor.Charge` — charge state (optional)
//! - `Precursor.Mz` — precursor m/z (optional)
//! - `RT` — retention time in minutes (optional)
//! - `Q.Value` — FDR q-value (optional)
//! - `Run` — raw/spectrum file name (optional)

use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, Float32Array, Float64Array, Int32Array, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use tracing;

use crate::error::EntrapmentError;
use crate::mod_parser::parse_modified_sequence;
use crate::types::UnifiedPsm;

/// Load DIA-NN parquet results into a vector of [`UnifiedPsm`].
///
/// Essential columns (`Stripped.Sequence`, `Protein.Ids`) must be present;
/// all other columns are optional and default to `None` when missing.
pub fn load_diann_parquet(path: &Path) -> Result<Vec<UnifiedPsm>, EntrapmentError> {
    let file = std::fs::File::open(path).map_err(|e| EntrapmentError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;

    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
        EntrapmentError::LoaderError {
            format: "DIA-NN parquet".to_string(),
            detail: format!("failed to open parquet file '{}': {e}", path.display()),
        }
    })?;

    let reader = builder.build().map_err(|e| EntrapmentError::LoaderError {
        format: "DIA-NN parquet".to_string(),
        detail: format!(
            "failed to build parquet reader for '{}': {e}",
            path.display()
        ),
    })?;

    let mut psms = Vec::new();

    for batch_result in reader {
        let batch = batch_result.map_err(|e| EntrapmentError::LoaderError {
            format: "DIA-NN parquet".to_string(),
            detail: format!("failed to read record batch: {e}"),
        })?;

        let schema = batch.schema();

        // Essential columns
        let peptide_col = get_string_column(&batch, &schema, "Stripped.Sequence")?;
        let protein_col = get_string_column(&batch, &schema, "Protein.Ids")?;

        // Optional columns
        let charge_col = get_int_column_optional(&batch, &schema, "Precursor.Charge");
        let mz_col = get_float_column_optional(&batch, &schema, "Precursor.Mz");
        let rt_col = get_float_column_optional(&batch, &schema, "RT");
        let qvalue_col = get_float_column_optional(&batch, &schema, "Q.Value");
        let run_col = get_string_column_optional(&batch, &schema, "Run");
        let modified_seq_col =
            get_string_column_optional(&batch, &schema, "Modified.Sequence");

        for row in 0..batch.num_rows() {
            let peptide = get_str(&peptide_col, row).to_string();
            if peptide.is_empty() {
                continue;
            }

            let protein_ids = get_str(&protein_col, row).to_string();

            let charge = charge_col.as_ref().and_then(|c| get_i32(c, row));
            let precursor_mz = mz_col.as_ref().and_then(|c| get_f64(c, row));
            let retention_time = rt_col.as_ref().and_then(|c| get_f64(c, row));
            let q_value = qvalue_col.as_ref().and_then(|c| get_f64(c, row));
            let spectrum_file = run_col
                .as_ref()
                .map(|c| get_str(c, row))
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string());

            // Parse modifications from Modified.Sequence (if present)
            let modifications = modified_seq_col
                .as_ref()
                .map(|c| get_str(c, row))
                .filter(|s| !s.is_empty())
                .map(|s| {
                    let (_stripped, mods) = parse_modified_sequence(s);
                    mods.iter()
                        .map(|m| (m.position, m.delta_mass))
                        .collect::<Vec<(usize, f64)>>()
                })
                .unwrap_or_default();

            psms.push(UnifiedPsm {
                peptide,
                charge,
                precursor_mz,
                retention_time,
                scan_number: None, // DIA-NN doesn't have scan numbers
                spectrum_file,
                protein_ids,
                q_value,
                modifications,
            });
        }
    }

    tracing::info!(
        "loaded {} PSMs from DIA-NN parquet: {}",
        psms.len(),
        path.display()
    );
    Ok(psms)
}

// ── Column helper functions ─────────────────────────────────────────

fn get_string_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<StringArray>, EntrapmentError> {
    let idx = schema
        .index_of(name)
        .map_err(|_| EntrapmentError::LoaderError {
            format: "DIA-NN parquet".to_string(),
            detail: format!("missing essential column '{name}'"),
        })?;
    let arr = batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| EntrapmentError::LoaderError {
            format: "DIA-NN parquet".to_string(),
            detail: format!("column '{name}' is not a String type"),
        })?;
    Ok(Arc::new(arr.clone()))
}

fn get_string_column_optional(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Option<Arc<StringArray>> {
    let idx = schema.index_of(name).ok()?;
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|a| Arc::new(a.clone()))
}

fn get_int_column_optional(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Option<Arc<dyn Array>> {
    let idx = schema.index_of(name).ok()?;
    Some(batch.column(idx).clone())
}

fn get_float_column_optional(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Option<Arc<dyn Array>> {
    let idx = schema.index_of(name).ok()?;
    Some(batch.column(idx).clone())
}

fn get_str(col: &Arc<StringArray>, row: usize) -> &str {
    if col.is_null(row) {
        ""
    } else {
        col.value(row)
    }
}

fn get_i32(col: &Arc<dyn Array>, row: usize) -> Option<i32> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<Int32Array>() {
        return Some(a.value(row));
    }
    if let Some(a) = col.as_any().downcast_ref::<arrow::array::Int64Array>() {
        return i32::try_from(a.value(row)).ok();
    }
    if let Some(a) = col.as_any().downcast_ref::<arrow::array::Int16Array>() {
        return Some(i32::from(a.value(row)));
    }
    None
}

fn get_f64(col: &Arc<dyn Array>, row: usize) -> Option<f64> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<Float64Array>() {
        return Some(a.value(row));
    }
    if let Some(a) = col.as_any().downcast_ref::<Float32Array>() {
        return Some(a.value(row) as f64);
    }
    None
}

#[cfg(test)]
mod tests {
    use crate::mod_parser::parse_modified_sequence;

    #[test]
    fn test_modifications_from_parsed() {
        let (stripped, mods) = parse_modified_sequence("AAAC(UniMod:4)DFK");
        assert_eq!(stripped, "AAACDFK");
        let modifications: Vec<(usize, f64)> =
            mods.iter().map(|m| (m.position, m.delta_mass)).collect();
        assert_eq!(modifications.len(), 1);
        assert_eq!(modifications[0].0, 3);
        assert!((modifications[0].1 - 57.021464).abs() < 1e-6);
    }

    #[test]
    fn test_modifications_empty_when_unmodified() {
        let (stripped, mods) = parse_modified_sequence("PEPTIDE");
        assert_eq!(stripped, "PEPTIDE");
        let modifications: Vec<(usize, f64)> =
            mods.iter().map(|m| (m.position, m.delta_mass)).collect();
        assert!(modifications.is_empty());
    }

    #[test]
    fn test_modifications_multiple_mods() {
        let (stripped, mods) = parse_modified_sequence("AC(UniMod:4)DEFM(UniMod:35)GK");
        assert_eq!(stripped, "ACDEFMGK");
        let modifications: Vec<(usize, f64)> =
            mods.iter().map(|m| (m.position, m.delta_mass)).collect();
        assert_eq!(modifications.len(), 2);
        assert_eq!(modifications[0].0, 1); // C at position 1
        assert!((modifications[0].1 - 57.021464).abs() < 1e-6);
        assert_eq!(modifications[1].0, 5); // M at position 5
        assert!((modifications[1].1 - 15.994915).abs() < 1e-6);
    }
}
