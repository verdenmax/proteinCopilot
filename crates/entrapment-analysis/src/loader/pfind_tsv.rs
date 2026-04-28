//! Loader for pFind result TSV files into [`UnifiedPsm`] records.
//!
//! Parses pFind's 16-column TSV format (tab-separated, with header row).
//! Auto-detects pFind format via [`detect`] by probing the header for
//! `PeptideSequence`, `ScanNo`, and `FileName` columns.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use crate::error::EntrapmentError;
use crate::types::UnifiedPsm;

/// Proton mass constant used for MH+ → m/z conversion.
const PROTON_MASS: f64 = 1.007276;

/// Check whether the file at `path` is a pFind TSV by inspecting its header.
///
/// Returns `true` if the first line contains all three marker columns:
/// `PeptideSequence`, `ScanNo`, and `FileName`.
pub fn detect(path: &Path) -> bool {
    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return false,
    };
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();
    let header = match lines.next() {
        Some(Ok(line)) => line,
        _ => return false,
    };
    let fields: Vec<&str> = header.split('\t').collect();
    fields.contains(&"PeptideSequence")
        && fields.contains(&"ScanNo")
        && fields.contains(&"FileName")
}

/// Build the modification name → delta mass lookup table.
fn mod_mass_table() -> HashMap<&'static str, f64> {
    let mut m = HashMap::new();
    m.insert("Carbamidomethyl", 57.021464);
    m.insert("Oxidation", 15.994915);
    m.insert("Acetyl", 42.010565);
    m.insert("Phospho", 79.966331);
    m.insert("Deamidated", 0.984016);
    m.insert("Methyl", 14.015650);
    m.insert("Dimethyl", 28.031300);
    m.insert("Trimethyl", 42.046950);
    m.insert("GlyGly", 114.042927);
    m.insert("Formyl", 27.994915);
    m.insert("Dehydrated", -18.010565);
    m.insert("Glu->pyro-Glu", -18.010565);
    m.insert("Gln->pyro-Glu", -17.026549);
    m.insert("Carbamyl", 43.005814);
    m
}

/// Parse pFind modification string into `Vec<(usize, f64)>`.
///
/// pFind format: `pos,Name[Residue];pos,Name[Residue];`
/// - `pos` is 1-based for residue mods, 0 for N-term
/// - Converts to 0-based: `pos - 1` for pos > 0, keep 0 for N-term
///
/// Unknown modification names are skipped with a tracing warning.
fn parse_modifications(mods_str: &str) -> Vec<(usize, f64)> {
    let mods_str = mods_str.trim();
    if mods_str.is_empty() {
        return Vec::new();
    }

    let table = mod_mass_table();
    let mut result = Vec::new();

    for entry in mods_str.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        // Format: "pos,Name[Residue]"
        let comma_pos = match entry.find(',') {
            Some(p) => p,
            None => {
                tracing::warn!(entry = %entry, "unexpected pFind mod entry format (no comma)");
                continue;
            }
        };
        let pos_str = &entry[..comma_pos];
        let rest = &entry[comma_pos + 1..];

        let pos: usize = match pos_str.parse() {
            Ok(p) => p,
            Err(_) => {
                tracing::warn!(entry = %entry, "cannot parse modification position");
                continue;
            }
        };

        // Extract mod name: everything before '['
        let name = match rest.find('[') {
            Some(bracket) => &rest[..bracket],
            None => rest,
        };

        let delta = match table.get(name) {
            Some(&d) => d,
            None => {
                tracing::warn!(name = %name, entry = %entry, "unknown pFind modification name");
                continue;
            }
        };

        // Convert 1-based to 0-based; N-term (0) stays 0
        let zero_based = if pos > 0 { pos - 1 } else { 0 };
        result.push((zero_based, delta));
    }

    result
}

/// Parse pFind protein string: split by `/`, filter empty, rejoin with `;`.
fn parse_proteins(proteins_str: &str) -> String {
    proteins_str
        .split('/')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join(";")
}

/// Load a pFind TSV result file into a vector of [`UnifiedPsm`].
///
/// Reads columns: `PeptideSequence`, `Charge`, `MH+`, `PredRT`, `ScanNo`,
/// `FileName`, `Proteins`, `QValue`, `Modifications`.
///
/// MH+ is converted to precursor m/z via:
/// `(MH+ + (charge - 1) * 1.007276) / charge`
pub fn load_pfind_tsv(path: &Path) -> Result<Vec<UnifiedPsm>, EntrapmentError> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(|e| EntrapmentError::IoError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;

    let headers = rdr.headers().map_err(|e| EntrapmentError::LoaderError {
        format: "pFind TSV".to_string(),
        detail: format!("failed to read headers from '{}': {e}", path.display()),
    })?;

    let find_col = |name: &str| -> Option<usize> { headers.iter().position(|h| h == name) };

    let peptide_idx =
        find_col("PeptideSequence").ok_or_else(|| EntrapmentError::LoaderError {
            format: "pFind TSV".to_string(),
            detail: format!(
                "missing essential column 'PeptideSequence' in '{}'",
                path.display()
            ),
        })?;

    let charge_idx = find_col("Charge");
    let mh_plus_idx = find_col("MH+");
    let pred_rt_idx = find_col("PredRT");
    let scan_no_idx = find_col("ScanNo");
    let file_name_idx = find_col("FileName");
    let proteins_idx = find_col("Proteins");
    let q_value_idx = find_col("QValue");
    let mods_idx = find_col("Modifications");

    let mut psms = Vec::new();

    for result in rdr.records() {
        let record = result.map_err(|e| EntrapmentError::LoaderError {
            format: "pFind TSV".to_string(),
            detail: format!("failed to parse row in '{}': {e}", path.display()),
        })?;

        let peptide = record.get(peptide_idx).unwrap_or("").trim().to_string();
        if peptide.is_empty() {
            continue;
        }

        let charge = charge_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<i32>().ok());

        let mh_plus = mh_plus_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        // Convert MH+ to precursor m/z: (MH+ + (charge-1) * proton) / charge
        let precursor_mz = match (mh_plus, charge) {
            (Some(mh), Some(z)) if z > 0 => {
                Some((mh + (z as f64 - 1.0) * PROTON_MASS) / z as f64)
            }
            _ => None,
        };

        let retention_time = pred_rt_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        let scan_number = scan_no_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<u32>().ok());

        let spectrum_file = file_name_idx
            .and_then(|i| record.get(i))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());

        let protein_ids = proteins_idx
            .and_then(|i| record.get(i))
            .map(parse_proteins)
            .unwrap_or_default();

        let q_value = q_value_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.trim().parse::<f64>().ok());

        let modifications = mods_idx
            .and_then(|i| record.get(i))
            .map(parse_modifications)
            .unwrap_or_default();

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
            modifications,
        });
    }

    tracing::info!(
        "loaded {} PSMs from pFind TSV: {}",
        psms.len(),
        path.display()
    );
    Ok(psms)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_mods_empty() {
        let mods = parse_modifications("");
        assert!(mods.is_empty());

        let mods = parse_modifications("   ");
        assert!(mods.is_empty());
    }

    #[test]
    fn test_parse_mods_single() {
        // "10,Carbamidomethyl[C];" → (9, 57.021464)
        let mods = parse_modifications("10,Carbamidomethyl[C];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].0, 9); // 10 - 1 = 9 (0-based)
        assert!((mods[0].1 - 57.021464).abs() < 1e-6);
    }

    #[test]
    fn test_parse_mods_double() {
        // "3,Carbamidomethyl[C];5,Carbamidomethyl[C];"
        let mods = parse_modifications("3,Carbamidomethyl[C];5,Carbamidomethyl[C];");
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].0, 2); // 3 - 1
        assert!((mods[0].1 - 57.021464).abs() < 1e-6);
        assert_eq!(mods[1].0, 4); // 5 - 1
        assert!((mods[1].1 - 57.021464).abs() < 1e-6);
    }

    #[test]
    fn test_parse_mods_nterm() {
        // "0,Acetyl[ProteinN-term];" → (0, 42.010565)
        let mods = parse_modifications("0,Acetyl[ProteinN-term];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].0, 0); // N-term stays 0
        assert!((mods[0].1 - 42.010565).abs() < 1e-6);
    }

    #[test]
    fn test_parse_mods_oxidation() {
        let mods = parse_modifications("2,Oxidation[M];");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].0, 1); // 2 - 1
        assert!((mods[0].1 - 15.994915).abs() < 1e-6);
    }

    #[test]
    fn test_parse_proteins() {
        assert_eq!(parse_proteins("sp|P50475|SYAC_RAT/"), "sp|P50475|SYAC_RAT");
        assert_eq!(
            parse_proteins("sp|P1|A/sp|P2|B/"),
            "sp|P1|A;sp|P2|B"
        );
        assert_eq!(parse_proteins(""), "");
        assert_eq!(parse_proteins("/"), "");
    }

    #[test]
    fn test_detect_pfind_format() {
        let dir = tempfile::tempdir().unwrap(); // OK in test code
        let path = dir.path().join("pfind.tsv");

        // pFind header
        {
            let mut f = std::fs::File::create(&path).unwrap(); // OK in test code
            use std::io::Write;
            writeln!(f, "FileName\tPeptideSequence\tModifications\tPepMass\tPredRT\tCleavageType\tProNCTerm\tProteins\tMH+\tCharge\tScanNo\tRawScore\tDeltaMassPPM\tDeltaRT(Min)\tFinalScore\tQValue").unwrap();
        }
        assert!(detect(&path));

        // Generic header
        let path2 = dir.path().join("generic.tsv");
        {
            let mut f = std::fs::File::create(&path2).unwrap(); // OK in test code
            use std::io::Write;
            writeln!(f, "peptide\tcharge\tprotein_ids\tq_value").unwrap();
        }
        assert!(!detect(&path2));
    }

    #[test]
    fn test_detect_nonexistent_file() {
        assert!(!detect(Path::new("/nonexistent/file.tsv")));
    }

    #[test]
    fn test_mh_plus_to_mz_conversion() {
        // MH+ = 1012.470123, charge = 2
        // precursor_mz = (1012.470123 + 1 * 1.007276) / 2 = 506.738700
        let mh_plus = 1012.470123;
        let charge = 2;
        let mz = (mh_plus + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64;
        assert!((mz - 506.738700).abs() < 0.001);
    }

    #[test]
    fn test_load_fixture() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            // Skip if fixture not available
            return;
        }
        let psms = load_pfind_tsv(&fixture).unwrap(); // OK in test code
        assert_eq!(psms.len(), 10);

        // First row: HNDLDDVGK, charge 2, MH+ 1012.470123
        assert_eq!(psms[0].peptide, "HNDLDDVGK");
        assert_eq!(psms[0].charge, Some(2));
        assert!(psms[0].precursor_mz.is_some());
        let expected_mz = (1012.470123 + 1.0 * PROTON_MASS) / 2.0;
        assert!((psms[0].precursor_mz.unwrap() - expected_mz).abs() < 0.001); // OK in test
        assert_eq!(psms[0].scan_number, Some(9911));
        assert_eq!(psms[0].protein_ids, "sp|P50475|SYAC_RAT");
        assert_eq!(psms[0].q_value, Some(0.0));
        assert!(psms[0].modifications.is_empty());

        // Row with single mod: LWADHGVQACFGR, 10,Carbamidomethyl[C];
        assert_eq!(psms[2].peptide, "LWADHGVQACFGR");
        assert_eq!(psms[2].modifications.len(), 1);
        assert_eq!(psms[2].modifications[0].0, 9); // 10 - 1

        // Row with Oxidation: IMTYLFDFPHGPK, 2,Oxidation[M];
        assert_eq!(psms[6].peptide, "IMTYLFDFPHGPK");
        assert_eq!(psms[6].modifications.len(), 1);
        assert_eq!(psms[6].modifications[0].0, 1); // 2 - 1
        assert!((psms[6].modifications[0].1 - 15.994915).abs() < 1e-6);

        // Row with N-term Acetyl: ADQLTEEQIAEFK, 0,Acetyl[ProteinN-term];
        assert_eq!(psms[7].peptide, "ADQLTEEQIAEFK");
        assert_eq!(psms[7].modifications.len(), 1);
        assert_eq!(psms[7].modifications[0].0, 0); // N-term stays 0
        assert!((psms[7].modifications[0].1 - 42.010565).abs() < 1e-6);
        assert!(psms[7].q_value.unwrap() > 0.001); // OK in test

        // Row with double mod: ITCLCQVPQNAANR, 3,Carbamidomethyl[C];5,Carbamidomethyl[C];
        assert_eq!(psms[8].peptide, "ITCLCQVPQNAANR");
        assert_eq!(psms[8].modifications.len(), 2);
        assert_eq!(psms[8].modifications[0].0, 2); // 3 - 1
        assert_eq!(psms[8].modifications[1].0, 4); // 5 - 1

        // Last row: different spectrum file
        assert_eq!(psms[9].peptide, "VLGGLGK");
        assert!(psms[9]
            .spectrum_file
            .as_deref()
            .unwrap()
            .contains("600_650"));

        // Retention time check
        assert!((psms[0].retention_time.unwrap() - 15.122).abs() < 0.001); // OK in test
    }
}
