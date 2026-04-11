//! Parser for DIA-NN report.parquet files.
//!
//! Key columns:
//! - `Modified.Sequence`: `_AAAC(UniMod:4)DM(UniMod:35)K_`
//! - `Precursor.Charge`: charge state
//! - `Precursor.Mz`: precursor m/z
//! - `RT`: retention time in minutes
//! - `Q.Value`: FDR q-value
//! - `Run`: raw file name
//! - `Protein.Names`: protein accessions (semicolon-separated)

use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, Float32Array, Float64Array, Int32Array, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use regex::Regex;

use protein_copilot_core::search_params::{ModPosition, Modification};

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

pub struct DiannParser {
    /// Maximum Q.Value to include. PSMs above this threshold are filtered out.
    pub filter_qvalue: Option<f64>,
    /// Optional: only import PSMs from this specific Run.
    pub run_filter: Option<String>,
}

impl Default for DiannParser {
    fn default() -> Self {
        Self {
            filter_qvalue: Some(0.01),
            run_filter: None,
        }
    }
}

impl DiannParser {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Parse DIA-NN `Modified.Sequence` like `_AAAC(UniMod:4)DM(UniMod:35)K_`
///
/// Returns `(clean_sequence, Vec<(position_1based, unimod_id)>)`.
/// N-terminal modifications (before any residue) get position 0.
fn parse_modified_sequence(modified_seq: &str) -> (String, Vec<(usize, u32)>) {
    let re = Regex::new(r"\(UniMod:(\d+)\)").unwrap();
    let mut clean = String::new();
    let mut mods = Vec::new();

    // Strip leading/trailing underscores
    let trimmed = modified_seq.trim_matches('_');

    let mut pos = 0usize; // 1-based position in clean sequence
    let bytes = trimmed.as_bytes();

    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'(' {
            // Find closing paren
            if let Some(j) = trimmed[i..].find(')') {
                let inside = &trimmed[i..i + j + 1];
                if let Some(cap) = re.captures(inside) {
                    if let Ok(id) = cap[1].parse::<u32>() {
                        mods.push((pos, id));
                    }
                }
                i += j + 1;
                continue;
            }
        }

        if bytes[i].is_ascii_uppercase() {
            pos += 1;
            clean.push(bytes[i] as char);
        }
        i += 1;
    }

    (clean, mods)
}

/// Convert a (position, unimod_id) pair to a core Modification.
///
/// Position 0 is treated as N-terminal; position > 0 delegates to `UnimodDb`.
fn resolve_modification(
    unimod: &UnimodDb,
    pos: usize,
    id: u32,
    sequence: &str,
) -> Result<Modification, ResultImportError> {
    if pos == 0 {
        // N-terminal modification: no residue to attach to
        let entry = unimod
            .get(id)
            .ok_or(ResultImportError::UnknownUnimodId(id))?;
        Ok(Modification {
            name: entry.title.clone(),
            mass_delta: entry.mono_mass,
            residues: vec![],
            position: ModPosition::AnyNTerm,
        })
    } else {
        unimod.to_modification(id, pos, sequence)
    }
}

impl ResultParser for DiannParser {
    fn parse(
        &self,
        path: &Path,
        unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let file = std::fs::File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let reader = builder.build()?;

        let mut psms = Vec::new();
        let mut filtered_count = 0u64;

        let mut null_skip_count = 0u64;

        for batch_result in reader {
            let batch = batch_result?;
            let schema = batch.schema();

            let mod_seq_col = get_string_column(&batch, &schema, "Modified.Sequence")?;
            let charge_col = get_int_column(&batch, &schema, "Precursor.Charge")?;
            let mz_col = get_float_column(&batch, &schema, "Precursor.Mz")?;
            let rt_col = get_float_column(&batch, &schema, "RT")?;
            let qvalue_col = get_float_column(&batch, &schema, "Q.Value")?;
            let run_col = get_string_column(&batch, &schema, "Run")?;
            let protein_col = get_string_column_optional(&batch, &schema, "Protein.Names");

            for row in 0..batch.num_rows() {
                let qvalue = match get_f64(&qvalue_col, row) {
                    Some(v) => v,
                    None => { null_skip_count += 1; continue; }
                };
                if let Some(max_q) = self.filter_qvalue {
                    if qvalue > max_q {
                        filtered_count += 1;
                        continue;
                    }
                }

                let run = get_str(&run_col, row);
                if run.is_empty() {
                    null_skip_count += 1;
                    continue;
                }
                let run = run.to_string();

                if let Some(ref filter) = self.run_filter {
                    if &run != filter {
                        continue;
                    }
                }

                let mod_seq = get_str(&mod_seq_col, row);
                if mod_seq.is_empty() {
                    null_skip_count += 1;
                    continue;
                }
                let (sequence, mod_positions) = parse_modified_sequence(mod_seq);

                let charge = match get_i32(&charge_col, row) {
                    Some(c) if c > 0 => c,
                    _ => { null_skip_count += 1; continue; }
                };
                let precursor_mz = match get_f64(&mz_col, row) {
                    Some(v) if v > 0.0 => v,
                    _ => { null_skip_count += 1; continue; }
                };
                let rt_min = match get_f64(&rt_col, row) {
                    Some(v) => v,
                    None => { null_skip_count += 1; continue; }
                };

                let mut modifications = Vec::new();
                for (pos, id) in &mod_positions {
                    match resolve_modification(unimod, *pos, *id, &sequence) {
                        Ok(m) => modifications.push(m),
                        Err(e) => {
                            tracing::debug!("DIA-NN mod conversion: {e}");
                        }
                    }
                }

                let proteins: Vec<String> = protein_col
                    .as_ref()
                    .map(|col| {
                        get_str(col, row)
                            .split(';')
                            .map(|s| s.trim().to_string())
                            .filter(|s| !s.is_empty())
                            .collect()
                    })
                    .unwrap_or_default();

                psms.push(ImportedPsm {
                    sequence,
                    charge,
                    precursor_mz,
                    rt_min: rt_min, // DIA-NN reports in minutes; store as-is
                    modifications,
                    score: Some(1.0 - qvalue), // invert: Q.Value (lower=better) → score (higher=better)
                    q_value: Some(qvalue),
                    protein_accessions: proteins,
                    raw_name: run,
                    matched_scan: None,
                    rt_delta_min: None,
                });
            }
        }

        if null_skip_count > 0 {
            tracing::warn!("skipped {null_skip_count} rows with null/invalid required fields");
        }

        if filtered_count > 0 {
            tracing::info!("filtered {filtered_count} PSMs above Q.Value threshold");
        }
        tracing::info!(
            "parsed {} PSMs from DIA-NN parquet: {}",
            psms.len(),
            path.display()
        );
        Ok(psms)
    }
}

// ── Column helper functions ─────────────────────────────────────────

fn get_string_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<StringArray>, ResultImportError> {
    let idx = schema
        .index_of(name)
        .map_err(|_| ResultImportError::MissingColumn {
            column: name.to_string(),
            expected: "Modified.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run"
                .to_string(),
        })?;
    Ok(Arc::new(
        batch
            .column(idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| {
                ResultImportError::Other(format!("column '{name}' is not String type"))
            })?
            .clone(),
    ))
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

fn get_int_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<dyn arrow::array::Array>, ResultImportError> {
    let idx = schema
        .index_of(name)
        .map_err(|_| ResultImportError::MissingColumn {
            column: name.to_string(),
            expected: "Modified.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run"
                .to_string(),
        })?;
    Ok(batch.column(idx).clone())
}

fn get_float_column(
    batch: &arrow::record_batch::RecordBatch,
    schema: &arrow::datatypes::SchemaRef,
    name: &str,
) -> Result<Arc<dyn arrow::array::Array>, ResultImportError> {
    let idx = schema
        .index_of(name)
        .map_err(|_| ResultImportError::MissingColumn {
            column: name.to_string(),
            expected: "Modified.Sequence, Precursor.Charge, Precursor.Mz, RT, Q.Value, Run"
                .to_string(),
        })?;
    Ok(batch.column(idx).clone())
}

fn get_str(col: &Arc<StringArray>, row: usize) -> &str {
    if col.is_null(row) {
        ""
    } else {
        col.value(row)
    }
}

fn get_i32(col: &Arc<dyn arrow::array::Array>, row: usize) -> Option<i32> {
    if col.is_null(row) {
        return None;
    }
    if let Some(a) = col.as_any().downcast_ref::<Int32Array>() {
        return Some(a.value(row));
    }
    if let Some(a) = col.as_any().downcast_ref::<arrow::array::Int64Array>() {
        return Some(a.value(row) as i32);
    }
    if let Some(a) = col.as_any().downcast_ref::<arrow::array::Int16Array>() {
        return Some(a.value(row) as i32);
    }
    None
}

fn get_f64(col: &Arc<dyn arrow::array::Array>, row: usize) -> Option<f64> {
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
    use super::*;

    #[test]
    fn parse_modified_sequence_no_mods() {
        let (seq, mods) = parse_modified_sequence("_PEPTIDE_");
        assert_eq!(seq, "PEPTIDE");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_modified_sequence_single_mod() {
        let (seq, mods) = parse_modified_sequence("_PEPTM(UniMod:35)IDE_");
        assert_eq!(seq, "PEPTMIDE");
        assert_eq!(mods, vec![(5, 35)]);
    }

    #[test]
    fn parse_modified_sequence_multiple_mods() {
        let (seq, mods) = parse_modified_sequence("_AAAC(UniMod:4)DM(UniMod:35)K_");
        assert_eq!(seq, "AAACDMK");
        assert_eq!(mods, vec![(4, 4), (6, 35)]);
    }

    #[test]
    fn parse_modified_sequence_nterm_mod() {
        let (seq, mods) = parse_modified_sequence("_(UniMod:1)PEPTIDE_");
        assert_eq!(seq, "PEPTIDE");
        assert_eq!(mods, vec![(0, 1)]);
    }

    #[test]
    fn parse_modified_sequence_bare() {
        let (seq, mods) = parse_modified_sequence("PEPTIDE");
        assert_eq!(seq, "PEPTIDE");
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_real_diann_parquet() {
        let path = Path::new("/home/verden/pfind/2025-fall/code/report.parquet");
        if !path.exists() {
            eprintln!("skipping parquet test: report.parquet not found");
            return;
        }
        let db = UnimodDb::builtin();
        let mut parser = DiannParser::new();
        parser.filter_qvalue = Some(0.001);
        let psms = parser.parse(path, &db).unwrap();
        assert!(!psms.is_empty(), "should parse some PSMs");
        for psm in &psms[..5.min(psms.len())] {
            assert!(
                psm.rt_min > 1.0,
                "RT should be in minutes (>1), got {}",
                psm.rt_min
            );
        }
        tracing::info!("parsed {} PSMs from real DIA-NN parquet", psms.len());
    }
}
