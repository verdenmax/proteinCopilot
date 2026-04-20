//! Loader for DIA-NN `report.parquet` files into [`UnifiedPsm`] records.
//!
//! Reads the following columns:
//! - `Stripped.Sequence` — stripped peptide sequence (essential)
//! - `Protein.Ids` — semicolon-separated protein accessions (essential)
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

            psms.push(UnifiedPsm {
                peptide,
                charge,
                precursor_mz,
                retention_time,
                scan_number: None, // DIA-NN doesn't have scan numbers
                spectrum_file,
                protein_ids,
                q_value,
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
