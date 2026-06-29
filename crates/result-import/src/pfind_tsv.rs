//! Parser for pFind TSV result files (tab-separated, 16 columns).
//!
//! pFind results include scan numbers directly, so scan matching
//! is not required — `matched_scan` is set from the `ScanNo` column.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use protein_copilot_core::search_params::{ModPosition, Modification};

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

/// Builtin lookup table for common modification masses.
fn builtin_mod_masses() -> HashMap<&'static str, f64> {
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

/// Proton mass in Da.
const PROTON_MASS: f64 = 1.007276;

/// Parser for pFind TSV result files.
pub struct PFindTsvParser;

impl ResultParser for PFindTsvParser {
    fn parse(
        &self,
        path: &Path,
        _unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        let _span = tracing::info_span!("parse_pfind_tsv", file = %path.display()).entered();

        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        let mod_masses = builtin_mod_masses();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(path)
            .map_err(|e| ResultImportError::Other(format!("CSV reader error: {e}")))?;

        let headers = reader
            .headers()
            .map_err(|e| ResultImportError::Other(format!("CSV header error: {e}")))?
            .clone();

        let col = |name: &str| -> Result<usize, ResultImportError> {
            headers.iter().position(|h| h == name).ok_or_else(|| {
                ResultImportError::MissingColumn {
                    column: name.to_string(),
                    expected: "FileName, PeptideSequence, Modifications, PepMass, PredRT, CleavageType, ProNCTerm, Proteins, MH+, Charge, ScanNo, RawScore, DeltaMassPPM, DeltaRT(Min), FinalScore, QValue".to_string(),
                }
            })
        };

        let i_file = col("FileName")?;
        let i_seq = col("PeptideSequence")?;
        let i_mods = col("Modifications")?;
        let i_pred_rt = col("PredRT")?;
        let i_proteins = col("Proteins")?;
        let i_mh_plus = col("MH+")?;
        let i_charge = col("Charge")?;
        let i_scan = col("ScanNo")?;
        let i_final_score = col("FinalScore")?;
        let i_qvalue = col("QValue")?;
        let i_delta_rt = col("DeltaRT(Min)")?;

        let mut psms = Vec::new();

        for result in reader.records() {
            let record =
                result.map_err(|e| ResultImportError::Other(format!("CSV record error: {e}")))?;

            let file_name = record.get(i_file).unwrap_or_default();
            let sequence = record.get(i_seq).unwrap_or_default();
            let mod_str = record.get(i_mods).unwrap_or_default();
            let pred_rt: f64 = record.get(i_pred_rt).unwrap_or("0").parse().unwrap_or(0.0);
            let proteins_str = record.get(i_proteins).unwrap_or_default();
            let mh_plus: f64 = record.get(i_mh_plus).unwrap_or("0").parse().unwrap_or(0.0);
            let charge: i32 = record.get(i_charge).unwrap_or("0").parse().unwrap_or(0);
            let scan_no: u32 = record.get(i_scan).unwrap_or("0").parse().unwrap_or(0);
            let final_score: f64 = record
                .get(i_final_score)
                .unwrap_or("0")
                .parse()
                .unwrap_or(0.0);
            let q_value: f64 = record.get(i_qvalue).unwrap_or("0").parse().unwrap_or(0.0);
            let delta_rt: f64 = record.get(i_delta_rt).unwrap_or("0").parse().unwrap_or(0.0);

            if charge <= 0 {
                tracing::warn!(
                    scan = scan_no,
                    "skipping pFind PSM with non-positive charge"
                );
                continue;
            }

            let precursor_mz = mh_plus_to_precursor_mz(mh_plus, charge);
            let modifications = parse_modifications(mod_str, &mod_masses);
            let protein_accessions = parse_proteins(proteins_str);

            psms.push(ImportedPsm {
                sequence: sequence.to_string(),
                charge,
                precursor_mz,
                rt_min: pred_rt,
                modifications,
                score: Some(final_score),
                q_value: Some(q_value),
                protein_accessions,
                raw_name: file_name.to_string(),
                matched_scan: Some(scan_no),
                rt_delta_min: Some(delta_rt),
            });
        }

        tracing::info!(count = psms.len(), "parsed pFind TSV PSMs");
        Ok(psms)
    }
}

/// Convert MH+ to precursor m/z.
///
/// Formula: `(MH+ + (charge - 1) * proton_mass) / charge`
fn mh_plus_to_precursor_mz(mh_plus: f64, charge: i32) -> f64 {
    (mh_plus + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64
}

/// Parse pFind modification string.
///
/// Format: `pos,Name[Residue];pos,Name[Residue];`
/// - Empty string → empty vec
/// - `10,Carbamidomethyl[C];` → single mod at position 10
/// - `0,Acetyl[ProteinN-term];` → N-term mod
fn parse_modifications(mod_str: &str, mod_masses: &HashMap<&str, f64>) -> Vec<Modification> {
    let trimmed = mod_str.trim();
    if trimmed.is_empty() {
        return Vec::new();
    }

    let mut mods = Vec::new();

    for entry in trimmed.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }

        // Format: "pos,Name[Residue]"
        let Some((pos_str, rest)) = entry.split_once(',') else {
            continue;
        };

        let _pos: usize = pos_str.parse().unwrap_or(0);

        // Extract name and residue from "Name[Residue]"
        let Some((name, residue_bracket)) = rest.split_once('[') else {
            continue;
        };
        let residue_str = residue_bracket.trim_end_matches(']');

        let mass_delta = mod_masses.get(name).copied().unwrap_or(0.0);

        let (residues, position) = if residue_str == "ProteinN-term" {
            (Vec::new(), ModPosition::ProteinNTerm)
        } else if residue_str == "ProteinC-term" {
            (Vec::new(), ModPosition::ProteinCTerm)
        } else if residue_str == "N-term" || residue_str == "AnyN-term" {
            (Vec::new(), ModPosition::AnyNTerm)
        } else if residue_str == "C-term" || residue_str == "AnyC-term" {
            (Vec::new(), ModPosition::AnyCTerm)
        } else {
            let residues: Vec<char> = residue_str.chars().collect();
            (residues, ModPosition::Anywhere)
        };

        mods.push(Modification {
            name: name.to_string(),
            mass_delta,
            residues,
            position,
        });
    }

    mods
}

/// Parse pFind protein string.
///
/// Format: `sp|P50475|SYAC_RAT/sp|Q9JMG1|EDF1_MOUSE/`
/// Split by `/`, remove trailing empty entries.
fn parse_proteins(proteins_str: &str) -> Vec<String> {
    proteins_str
        .split('/')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

/// Detect whether a file is a pFind TSV file by checking its header.
///
/// Returns `true` if the first line contains `PeptideSequence`, `ScanNo`, and `FileName`.
pub fn detect(path: &Path) -> bool {
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let reader = std::io::BufReader::new(file);
    let Some(Ok(header)) = reader.lines().next() else {
        return false;
    };
    header.contains("PeptideSequence") && header.contains("ScanNo") && header.contains("FileName")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mod_masses() -> HashMap<&'static str, f64> {
        builtin_mod_masses()
    }

    #[test]
    fn parse_empty_mods() {
        let mods = parse_modifications("", &mod_masses());
        assert!(mods.is_empty());
    }

    #[test]
    fn parse_single_mod() {
        let mods = parse_modifications("10,Carbamidomethyl[C];", &mod_masses());
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].name, "Carbamidomethyl");
        assert!((mods[0].mass_delta - 57.021464).abs() < 1e-6);
        assert_eq!(mods[0].residues, vec!['C']);
        assert_eq!(mods[0].position, ModPosition::Anywhere);
    }

    #[test]
    fn parse_double_mod() {
        let mods = parse_modifications("3,Carbamidomethyl[C];5,Carbamidomethyl[C];", &mod_masses());
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].name, "Carbamidomethyl");
        assert_eq!(mods[1].name, "Carbamidomethyl");
    }

    #[test]
    fn parse_oxidation_mod() {
        let mods = parse_modifications("2,Oxidation[M];", &mod_masses());
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].name, "Oxidation");
        assert!((mods[0].mass_delta - 15.994915).abs() < 1e-6);
        assert_eq!(mods[0].residues, vec!['M']);
    }

    #[test]
    fn parse_nterm_mod() {
        let mods = parse_modifications("0,Acetyl[ProteinN-term];", &mod_masses());
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].name, "Acetyl");
        assert!((mods[0].mass_delta - 42.010565).abs() < 1e-6);
        assert!(mods[0].residues.is_empty());
        assert_eq!(mods[0].position, ModPosition::ProteinNTerm);
    }

    #[test]
    fn parse_proteins_single() {
        let proteins = parse_proteins("sp|P50475|SYAC_RAT/");
        assert_eq!(proteins, vec!["sp|P50475|SYAC_RAT"]);
    }

    #[test]
    fn parse_proteins_multiple() {
        let proteins = parse_proteins("sp|P50475|SYAC_RAT/sp|Q9JMG1|EDF1_MOUSE/");
        assert_eq!(proteins, vec!["sp|P50475|SYAC_RAT", "sp|Q9JMG1|EDF1_MOUSE"]);
    }

    #[test]
    fn parse_proteins_empty() {
        let proteins = parse_proteins("");
        assert!(proteins.is_empty());
    }

    #[test]
    fn precursor_mz_calculation() {
        // MH+ = 1012.470123, charge = 2
        // (1012.470123 + 1 * 1.007276) / 2 = 506.738700
        let mz = mh_plus_to_precursor_mz(1012.470123, 2);
        assert!((mz - 506.738700).abs() < 0.001);
    }

    #[test]
    fn precursor_mz_charge3() {
        // MH+ = 1592.692795, charge = 3
        // (1592.692795 + 2 * 1.007276) / 3 = 531.569116
        let mz = mh_plus_to_precursor_mz(1592.692795, 3);
        assert!((mz - 531.569116).abs() < 0.001);
    }

    #[test]
    fn parse_skips_non_positive_charge() {
        // One charge-0 row (must be skipped) followed by one valid charge-2 row.
        let header = "FileName\tPeptideSequence\tModifications\tPepMass\tPredRT\tCleavageType\tProNCTerm\tProteins\tMH+\tCharge\tScanNo\tRawScore\tDeltaMassPPM\tDeltaRT(Min)\tFinalScore\tQValue";
        let row_charge0 =
            "raw1\tBADPEPK\t\t900.0\t10.0\tfull\t-\tsp|P1|A/\t901.0\t0\t100\t1.0\t0.5\t0.1\t5.0\t0.01";
        let row_valid =
            "raw1\tHNDLDDVGK\t\t1000.0\t12.0\tfull\t-\tsp|P2|B/\t1012.470123\t2\t200\t1.0\t0.5\t0.1\t6.0\t0.01";
        let content = format!("{header}\n{row_charge0}\n{row_valid}\n");

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pfind_charge0.tsv");
        std::fs::write(&path, content).unwrap();

        let db = UnimodDb::builtin();
        let parser = PFindTsvParser;
        let psms = parser.parse(&path, &db).expect("should parse pFind TSV");

        assert_eq!(psms.len(), 1, "the charge-0 PSM must be skipped");
        assert_eq!(psms[0].sequence, "HNDLDDVGK");
        assert_eq!(psms[0].charge, 2);
        assert!(
            psms[0].precursor_mz.is_finite() && psms[0].precursor_mz > 0.0,
            "valid PSM must have a finite, positive precursor m/z, got {}",
            psms[0].precursor_mz
        );
    }

    #[test]
    fn parse_fixture_file() {
        let fixture =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            // Skip if fixture not available
            return;
        }

        let db = UnimodDb::builtin();
        let parser = PFindTsvParser;
        let psms = parser.parse(&fixture, &db).expect("should parse fixture");

        assert_eq!(psms.len(), 10);

        // First row: HNDLDDVGK, no mods, charge 2, scan 9911
        assert_eq!(psms[0].sequence, "HNDLDDVGK");
        assert!(psms[0].modifications.is_empty());
        assert_eq!(psms[0].charge, 2);
        assert_eq!(psms[0].matched_scan, Some(9911));

        // Third row: LWADHGVQACFGR, Carbamidomethyl[C] at pos 10
        assert_eq!(psms[2].sequence, "LWADHGVQACFGR");
        assert_eq!(psms[2].modifications.len(), 1);
        assert_eq!(psms[2].modifications[0].name, "Carbamidomethyl");

        // Row 7 (index 6): Oxidation[M]
        assert_eq!(psms[6].sequence, "IMTYLFDFPHGPK");
        assert_eq!(psms[6].modifications.len(), 1);
        assert_eq!(psms[6].modifications[0].name, "Oxidation");

        // Row 8 (index 7): Acetyl[ProteinN-term]
        assert_eq!(psms[7].sequence, "ADQLTEEQIAEFK");
        assert_eq!(psms[7].modifications[0].name, "Acetyl");
        assert_eq!(psms[7].modifications[0].position, ModPosition::ProteinNTerm);

        // Row 9 (index 8): double Carbamidomethyl
        assert_eq!(psms[8].sequence, "ITCLCQVPQNAANR");
        assert_eq!(psms[8].modifications.len(), 2);

        // Last row: charge 1, different raw file
        assert_eq!(psms[9].charge, 1);
        assert!(psms[9].raw_name.contains("600_650"));
    }

    #[test]
    fn detect_fixture_file() {
        let fixture =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/pfind_sample.tsv");
        if !fixture.exists() {
            return;
        }
        assert!(detect(&fixture));
    }

    #[test]
    fn detect_non_pfind_file() {
        // A non-existent file should return false
        assert!(!detect(Path::new("/nonexistent/file.tsv")));
    }
}
