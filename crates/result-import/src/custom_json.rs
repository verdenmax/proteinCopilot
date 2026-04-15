//! Parser for custom JSON format (hela.json style).
//!
//! Input format: JSON array of objects with fields:
//! - `sequence`: peptide sequence
//! - `charge`: charge state
//! - `modify`: `[[position, unimod_id], ...]` (1-based position)
//! - `rt`: retention time in minutes
//! - `precursor_mz`: precursor m/z
//! - `raw_title`: raw file name (without extension)
//! - `protein_names`: array of protein accessions

use std::path::Path;

use serde::Deserialize;

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

/// Raw JSON record matching the hela.json schema.
#[derive(Debug, Deserialize)]
struct RawJsonPsm {
    sequence: String,
    charge: i32,
    /// Modifications as [[position, unimod_id], ...]. May be empty or absent.
    #[serde(default)]
    modify: Vec<Vec<u32>>,
    /// Retention time in minutes.
    rt: f64,
    precursor_mz: f64,
    raw_title: String,
    #[serde(default)]
    protein_names: Vec<String>,
}

/// Parser for the custom JSON format.
pub struct CustomJsonParser;

impl ResultParser for CustomJsonParser {
    fn parse(&self, path: &Path, unimod: &UnimodDb) -> Result<Vec<ImportedPsm>, ResultImportError> {
        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let data = std::fs::read_to_string(path)?;
        let raw_psms: Vec<RawJsonPsm> = serde_json::from_str(&data)?;

        let mut psms = Vec::with_capacity(raw_psms.len());
        let mut mod_errors = 0u64;

        for (i, raw) in raw_psms.into_iter().enumerate() {
            let mut modifications = Vec::new();
            for pair in &raw.modify {
                if pair.len() != 2 {
                    tracing::warn!("PSM #{i}: invalid modify entry (expected [pos, id]): {pair:?}");
                    continue;
                }
                let position = pair[0] as usize;
                let unimod_id = pair[1];
                match unimod.to_modification(unimod_id, position, &raw.sequence) {
                    Ok(m) => modifications.push(m),
                    Err(e) => {
                        mod_errors += 1;
                        if mod_errors <= 5 {
                            tracing::warn!("PSM #{i} ({}): modification error: {e}", raw.sequence);
                        }
                    }
                }
            }

            psms.push(ImportedPsm {
                sequence: raw.sequence,
                charge: raw.charge,
                precursor_mz: raw.precursor_mz,
                rt_min: raw.rt, // already in minutes; store as-is
                modifications,
                score: None,
                q_value: None,
                protein_accessions: raw.protein_names,
                raw_name: raw.raw_title,
                matched_scan: None,
                rt_delta_min: None,
            });
        }

        if mod_errors > 5 {
            tracing::warn!(
                "... and {} more modification errors suppressed",
                mod_errors - 5
            );
        }
        tracing::info!(
            "parsed {} PSMs from custom JSON: {}",
            psms.len(),
            path.display()
        );
        Ok(psms)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn sample_json() -> &'static str {
        r#"[
            {
                "sequence": "AADLLDDVSQK",
                "charge": 2,
                "modify": [],
                "rt": 12.648,
                "precursor_mz": 588.2993,
                "raw_title": "hela_Rep1",
                "protein_names": ["sp|P12345|TEST_HUMAN"]
            },
            {
                "sequence": "PEPTMIDER",
                "charge": 3,
                "modify": [[5, 35]],
                "rt": 25.5,
                "precursor_mz": 372.185,
                "raw_title": "hela_Rep1",
                "protein_names": ["sp|P67890|TEST2_HUMAN"]
            },
            {
                "sequence": "ACDEFGHIK",
                "charge": 2,
                "modify": [[2, 4]],
                "rt": 40.0,
                "precursor_mz": 510.234,
                "raw_title": "hela_Rep2",
                "protein_names": ["sp|P11111|TEST3_HUMAN"]
            }
        ]"#
    }

    #[test]
    fn parse_custom_json_basic() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.json");
        std::fs::write(&path, sample_json()).unwrap();

        let db = UnimodDb::builtin();
        let parser = CustomJsonParser;
        let psms = parser.parse(&path, &db).unwrap();

        assert_eq!(psms.len(), 3);

        // First PSM: no modifications, RT converted to seconds
        assert_eq!(psms[0].sequence, "AADLLDDVSQK");
        assert_eq!(psms[0].charge, 2);
        assert!((psms[0].rt_min - 12.648).abs() < 0.01);
        assert!(psms[0].modifications.is_empty());
        assert_eq!(psms[0].raw_name, "hela_Rep1");

        // Second PSM: Oxidation at position 5 (M)
        assert_eq!(psms[1].sequence, "PEPTMIDER");
        assert_eq!(psms[1].modifications.len(), 1);
        assert_eq!(psms[1].modifications[0].name, "Oxidation");
        assert!((psms[1].modifications[0].mass_delta - 15.994915).abs() < 0.001);
        assert_eq!(psms[1].modifications[0].residues, vec!['M']);

        // Third PSM: Carbamidomethyl at position 2 (C)
        assert_eq!(psms[2].sequence, "ACDEFGHIK");
        assert_eq!(psms[2].modifications.len(), 1);
        assert_eq!(psms[2].modifications[0].name, "Carbamidomethyl");
        assert_eq!(psms[2].modifications[0].residues, vec!['C']);
        assert_eq!(psms[2].raw_name, "hela_Rep2");
    }

    #[test]
    fn parse_custom_json_rt_conversion() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.json");
        std::fs::write(&path, r#"[{"sequence":"PEPTIDE","charge":2,"modify":[],"rt":1.0,"precursor_mz":400.0,"raw_title":"test","protein_names":[]}]"#).unwrap();

        let db = UnimodDb::builtin();
        let psms = CustomJsonParser.parse(&path, &db).unwrap();
        // rt=1.0 means 1 minute, stored as 1.0 minutes
        assert!((psms[0].rt_min - 1.0).abs() < 0.01);
    }

    #[test]
    fn parse_custom_json_file_not_found() {
        let db = UnimodDb::builtin();
        let result = CustomJsonParser.parse(Path::new("/nonexistent.json"), &db);
        assert!(result.is_err());
    }
}
