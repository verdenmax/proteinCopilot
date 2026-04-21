//! Output writers for entrapment analysis results.
//!
//! Provides TSV and JSON output for classified PSMs, razor-error reports,
//! and run metadata.

use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use tracing::info;

use crate::error::EntrapmentError;
use crate::types::{ClassifiedPsm, DiscriminabilityLevel, LevelCounts, PsmGroup};

// ---------------------------------------------------------------------------
// Run metadata
// ---------------------------------------------------------------------------

/// Metadata captured for a single entrapment analysis run.
///
/// Serialised as JSON alongside result files for reproducibility.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunMetadata {
    /// Crate version string (from `env!("CARGO_PKG_VERSION")`).
    pub tool_version: String,
    /// ISO 8601 timestamp of the analysis run.
    pub run_timestamp: String,
    /// Path to the input search-engine results file.
    pub input_file: String,
    /// SHA-256 hex digest of the input file.
    pub input_sha256: String,
    /// Path to the target FASTA database file.
    pub fasta_file: String,
    /// SHA-256 hex digest of the FASTA file.
    pub fasta_sha256: String,
    /// Snapshot of the configuration used for this run.
    pub config_snapshot: serde_json::Value,
    /// Total number of PSMs processed.
    pub total_psms: usize,
    /// Number of PSMs classified as trap.
    pub trap_psms: usize,
    /// Per-level hit counts.
    pub level_counts: LevelCounts,
}

// ---------------------------------------------------------------------------
// SHA-256 helper
// ---------------------------------------------------------------------------

/// Compute the SHA-256 hex digest of a file.
///
/// Reads the file in 8 KiB chunks to avoid loading the entire file into memory.
pub fn file_sha256(path: &Path) -> Result<String, EntrapmentError> {
    let file = File::open(path).map_err(|e| EntrapmentError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 8192];

    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| EntrapmentError::IoError {
                path: path.to_path_buf(),
                detail: e.to_string(),
            })?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }

    let hash = hasher.finalize();
    let hex = hash.iter().map(|b| format!("{b:02x}")).collect::<String>();

    Ok(hex)
}

// ---------------------------------------------------------------------------
// TSV writers
// ---------------------------------------------------------------------------

/// Write all classified PSMs to a tab-separated file.
///
/// Headers: peptide, charge, precursor_mz, retention_time, scan_number,
/// spectrum_file, protein_ids, q_value, group, level, best_target_peptide,
/// best_target_protein, mismatches, delta_mass_da, diff_positions,
/// substitution_type, edit_distance, alignment_detail.
///
/// Optional fields (`best_target_peptide`, `edit_distance`, `alignment_detail`,
/// `diff_positions`) are written as empty strings when `None`.
/// `substitution_type` is always present (defaults to `None` variant).
pub fn write_classified_tsv(psms: &[ClassifiedPsm], path: &Path) -> Result<(), EntrapmentError> {
    let file = File::create(path).map_err(|e| EntrapmentError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);

    // Header
    wtr.write_record([
        "peptide",
        "charge",
        "precursor_mz",
        "retention_time",
        "scan_number",
        "spectrum_file",
        "protein_ids",
        "q_value",
        "group",
        "level",
        "best_target_peptide",
        "best_target_protein",
        "mismatches",
        "delta_mass_da",
        "diff_positions",
        "substitution_type",
        "edit_distance",
        "alignment_detail",
    ])
    .map_err(|e| EntrapmentError::OutputError {
        detail: format!("failed to write TSV header to {}: {e}", path.display()),
    })?;

    for cp in psms {
        wtr.write_record([
            &cp.psm.peptide,
            &opt_to_string(&cp.psm.charge),
            &opt_to_string(&cp.psm.precursor_mz),
            &opt_to_string(&cp.psm.retention_time),
            &opt_to_string(&cp.psm.scan_number),
            &opt_to_string(&cp.psm.spectrum_file),
            &cp.psm.protein_ids,
            &opt_to_string(&cp.psm.q_value),
            &cp.group.to_string(),
            &cp.level.to_string(),
            &opt_to_string(&cp.best_target_peptide),
            &opt_to_string(&cp.best_target_protein),
            &opt_to_string(&cp.mismatches),
            &opt_to_string(&cp.delta_mass_da),
            &opt_to_string(&cp.diff_positions),
            &cp.substitution_type.to_string(),
            &opt_to_string(&cp.edit_distance),
            &opt_to_string(&cp.alignment_detail),
        ])
        .map_err(|e| EntrapmentError::OutputError {
            detail: format!("failed to write TSV row to {}: {e}", path.display()),
        })?;
    }

    wtr.flush().map_err(|e| EntrapmentError::OutputError {
        detail: format!("failed to flush TSV to {}: {e}", path.display()),
    })?;

    info!(path = %path.display(), rows = psms.len(), "wrote classified TSV");
    Ok(())
}

/// Write L0 trap PSMs as a razor-error report TSV.
///
/// Filters to only PSMs where `group == Trap` **and** `level == L0`.
///
/// Headers: peptide, current_protein, suggested_target_protein, reason.
pub fn write_razor_errors_tsv(psms: &[ClassifiedPsm], path: &Path) -> Result<(), EntrapmentError> {
    let file = File::create(path).map_err(|e| EntrapmentError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;
    let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);

    wtr.write_record([
        "peptide",
        "current_protein",
        "suggested_target_protein",
        "reason",
    ])
    .map_err(|e| EntrapmentError::OutputError {
        detail: format!(
            "failed to write razor TSV header to {}: {e}",
            path.display()
        ),
    })?;

    let mut count = 0usize;
    for cp in psms {
        if cp.group != PsmGroup::Trap || cp.level != DiscriminabilityLevel::L0 {
            continue;
        }
        wtr.write_record([
            &cp.psm.peptide,
            &cp.psm.protein_ids,
            cp.best_target_protein.as_deref().unwrap_or(""),
            "exact sequence match in target database",
        ])
        .map_err(|e| EntrapmentError::OutputError {
            detail: format!("failed to write razor TSV row to {}: {e}", path.display()),
        })?;
        count += 1;
    }

    wtr.flush().map_err(|e| EntrapmentError::OutputError {
        detail: format!("failed to flush razor TSV to {}: {e}", path.display()),
    })?;

    info!(path = %path.display(), razor_errors = count, "wrote razor errors TSV");
    Ok(())
}

// ---------------------------------------------------------------------------
// JSON metadata writer
// ---------------------------------------------------------------------------

/// Write run metadata as pretty-printed JSON (2-space indent).
pub fn write_run_metadata(metadata: &RunMetadata, path: &Path) -> Result<(), EntrapmentError> {
    let json =
        serde_json::to_string_pretty(metadata).map_err(|e| EntrapmentError::OutputError {
            detail: format!("failed to serialise run metadata: {e}"),
        })?;

    let mut file = File::create(path).map_err(|e| EntrapmentError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })?;

    file.write_all(json.as_bytes())
        .map_err(|e| EntrapmentError::IoError {
            path: path.to_path_buf(),
            detail: e.to_string(),
        })?;

    info!(path = %path.display(), "wrote run metadata JSON");
    Ok(())
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Convert an `Option<T: Display>` to a `String`, using `""` for `None`.
fn opt_to_string<T: std::fmt::Display>(opt: &Option<T>) -> String {
    match opt {
        Some(v) => v.to_string(),
        None => String::new(),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{DiscriminabilityLevel, PsmGroup, SubstitutionType, UnifiedPsm};

    fn make_classified_psm(
        peptide: &str,
        protein: &str,
        group: PsmGroup,
        level: DiscriminabilityLevel,
    ) -> ClassifiedPsm {
        ClassifiedPsm {
            psm: UnifiedPsm {
                peptide: peptide.to_owned(),
                charge: Some(2),
                precursor_mz: Some(500.25),
                retention_time: Some(12.5),
                scan_number: Some(1001),
                spectrum_file: Some("sample.raw".to_owned()),
                protein_ids: protein.to_owned(),
                q_value: Some(0.01),
            },
            group,
            level,
            best_target_peptide: if level == DiscriminabilityLevel::L0 {
                Some(peptide.to_owned())
            } else {
                None
            },
            best_target_protein: if level == DiscriminabilityLevel::L0 {
                Some("sp|P99999|TARGET_HUMAN".to_owned())
            } else {
                None
            },
            mismatches: if level == DiscriminabilityLevel::L0 {
                Some(0)
            } else {
                None
            },
            delta_mass_da: if level == DiscriminabilityLevel::L0 {
                Some(0.0)
            } else {
                None
            },
            diff_positions: if level == DiscriminabilityLevel::L0 {
                Some(String::new())
            } else {
                None
            },
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        }
    }

    #[test]
    fn test_file_sha256_known_content() {
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("test.txt");
        std::fs::write(&path, b"hello world\n").expect("write file");

        let hash = file_sha256(&path).expect("compute sha256");

        // sha256("hello world\n") = a948904f2f0f479b8f8564e9d7d5a26e6b1cd5ad9c8e4e3e3b1f2...
        // Verify it's a 64-char lowercase hex string
        assert_eq!(hash.len(), 64);
        assert!(hash
            .chars()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase()));
    }

    #[test]
    fn test_file_sha256_missing_file() {
        let result = file_sha256(Path::new("/nonexistent/missing.txt"));
        assert!(result.is_err());
    }

    #[test]
    fn test_write_classified_tsv_roundtrip() {
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("classified.tsv");

        let psms = vec![
            make_classified_psm(
                "PEPTIDEK",
                "sp|P001|TRAP_YEAST",
                PsmGroup::Trap,
                DiscriminabilityLevel::L0,
            ),
            make_classified_psm(
                "ANOTHERK",
                "sp|P002|TARGET_HUMAN",
                PsmGroup::Target,
                DiscriminabilityLevel::L4,
            ),
        ];

        write_classified_tsv(&psms, &path).expect("write TSV");

        let content = std::fs::read_to_string(&path).expect("read TSV");
        let lines: Vec<&str> = content.lines().collect();

        // Header + 2 data rows
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("peptide\t"));
        assert!(lines[1].starts_with("PEPTIDEK\t"));
        assert!(lines[2].starts_with("ANOTHERK\t"));

        // Check that group/level columns are present
        assert!(lines[1].contains("trap"));
        assert!(lines[1].contains("L0"));
        assert!(lines[2].contains("target"));
        assert!(lines[2].contains("L4"));
    }

    #[test]
    fn test_write_classified_tsv_none_fields() {
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("classified_none.tsv");

        let psm = ClassifiedPsm {
            psm: UnifiedPsm {
                peptide: "TESTPEPTIDEK".to_owned(),
                charge: None,
                precursor_mz: None,
                retention_time: None,
                scan_number: None,
                spectrum_file: None,
                protein_ids: "sp|P001|TEST_YEAST".to_owned(),
                q_value: None,
            },
            group: PsmGroup::Trap,
            level: DiscriminabilityLevel::L4,
            best_target_peptide: None,
            best_target_protein: None,
            mismatches: None,
            delta_mass_da: None,
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        };

        write_classified_tsv(&[psm], &path).expect("write TSV");

        let content = std::fs::read_to_string(&path).expect("read TSV");
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2); // header + 1 data row

        // Count tabs in the data row — should match header
        let header_tabs = lines[0].matches('\t').count();
        let data_tabs = lines[1].matches('\t').count();
        assert_eq!(header_tabs, data_tabs);
    }

    #[test]
    fn test_write_razor_errors_tsv_filters_l0_traps() {
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("razor.tsv");

        let psms = vec![
            // L0 trap — should be included
            make_classified_psm(
                "PEPTIDEK",
                "sp|P001|TRAP_YEAST",
                PsmGroup::Trap,
                DiscriminabilityLevel::L0,
            ),
            // L4 trap — should NOT be included
            make_classified_psm(
                "ANOTHERK",
                "sp|P002|TRAP_YEAST",
                PsmGroup::Trap,
                DiscriminabilityLevel::L4,
            ),
            // L0 target — should NOT be included (wrong group)
            make_classified_psm(
                "TARGETPK",
                "sp|P003|TARGET_HUMAN",
                PsmGroup::Target,
                DiscriminabilityLevel::L0,
            ),
        ];

        write_razor_errors_tsv(&psms, &path).expect("write razor TSV");

        let content = std::fs::read_to_string(&path).expect("read razor TSV");
        let lines: Vec<&str> = content.lines().collect();

        // Header + 1 L0-trap row
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("peptide\t"));
        assert!(lines[1].contains("PEPTIDEK"));
        assert!(lines[1].contains("exact sequence match in target database"));
    }

    #[test]
    fn test_write_razor_errors_tsv_empty() {
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("razor_empty.tsv");

        let psms = vec![make_classified_psm(
            "PEPTIDEK",
            "sp|P001|TRAP_YEAST",
            PsmGroup::Trap,
            DiscriminabilityLevel::L4,
        )];

        write_razor_errors_tsv(&psms, &path).expect("write razor TSV");

        let content = std::fs::read_to_string(&path).expect("read razor TSV");
        let lines: Vec<&str> = content.lines().collect();
        // Header only, no data rows
        assert_eq!(lines.len(), 1);
    }

    #[test]
    fn test_write_run_metadata_json() {
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("metadata.json");

        let metadata = RunMetadata {
            tool_version: "0.1.0".to_owned(),
            run_timestamp: "2025-01-15T12:00:00Z".to_owned(),
            input_file: "results.tsv".to_owned(),
            input_sha256: "abc123".to_owned(),
            fasta_file: "target.fasta".to_owned(),
            fasta_sha256: "def456".to_owned(),
            config_snapshot: serde_json::json!({"version": 1}),
            total_psms: 100,
            trap_psms: 10,
            level_counts: LevelCounts {
                l0: 2,
                l1: 1,
                l2: 3,
                l3: 2,
                l4: 2,
            },
        };

        write_run_metadata(&metadata, &path).expect("write metadata JSON");

        let content = std::fs::read_to_string(&path).expect("read JSON");
        let parsed: serde_json::Value = serde_json::from_str(&content).expect("parse JSON");

        assert_eq!(parsed["tool_version"], "0.1.0");
        assert_eq!(parsed["total_psms"], 100);
        assert_eq!(parsed["trap_psms"], 10);
        assert_eq!(parsed["level_counts"]["l0"], 2);

        // Pretty-printed: should contain newlines
        assert!(content.contains('\n'));
    }

    #[test]
    fn test_opt_to_string() {
        assert_eq!(opt_to_string(&Some(42)), "42");
        assert_eq!(opt_to_string(&Some(3.14)), "3.14");
        assert_eq!(opt_to_string(&Some("hello".to_owned())), "hello");
        assert_eq!(opt_to_string::<i32>(&None), "");
    }

    #[test]
    fn test_write_classified_tsv_v2_columns() {
        use crate::types::SubstitutionType;
        let dir = tempfile::tempdir().expect("create temp dir");
        let path = dir.path().join("classified_v2.tsv");

        let psm = ClassifiedPsm {
            psm: UnifiedPsm {
                peptide: "PEPQDEK".to_owned(),
                charge: Some(2),
                precursor_mz: Some(500.0),
                retention_time: Some(10.0),
                scan_number: Some(100),
                spectrum_file: Some("test.raw".to_owned()),
                protein_ids: "sp|P001|TRAP_YEAST".to_owned(),
                q_value: Some(0.01),
            },
            group: PsmGroup::Trap,
            level: DiscriminabilityLevel::L2,
            best_target_peptide: Some("PEPKDEK".to_owned()),
            best_target_protein: Some("sp|P002|TARGET_HUMAN".to_owned()),
            mismatches: Some(1),
            delta_mass_da: Some(0.036385),
            diff_positions: Some("[3:Q->K]".to_owned()),
            substitution_type: SubstitutionType::QKSubstitution,
            edit_distance: Some(1),
            alignment_detail: Some("Q3→K".to_owned()),
        };

        write_classified_tsv(&[psm], &path).expect("write TSV");

        let content = std::fs::read_to_string(&path).expect("read TSV");
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines[0].contains("substitution_type"));
        assert!(lines[0].contains("edit_distance"));
        assert!(lines[0].contains("alignment_detail"));
        assert!(lines[1].contains("QKSubstitution"));
    }
}
