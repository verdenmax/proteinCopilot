//! HTML report generation for entrapment analysis results.
//!
//! Renders classified PSMs into a self-contained interactive HTML dashboard
//! using Plotly.js for charts and vanilla JS for table search/sort.

use std::path::Path;

use serde::Serialize;
use tracing::info;

use crate::error::EntrapmentError;
use crate::types::{ClassifiedPsm, EntrapmentSummary, PsmGroup};

/// The HTML template is compiled into the binary at build time.
const TEMPLATE: &str = include_str!("../templates/entrapment_report.html");

/// Placeholder token inside the template that is replaced with JSON data.
const DATA_PLACEHOLDER: &str = "/*__REPORT_DATA__*/{}";

/// Top-level JSON structure embedded into the HTML template.
#[derive(Debug, Serialize)]
struct ReportData {
    summary: EntrapmentSummary,
    psms: Vec<PsmRow>,
}

/// A flat row for the PSM table in the HTML report.
#[derive(Debug, Serialize)]
struct PsmRow {
    peptide: String,
    group: String,
    level: String,
    best_target_peptide: String,
    best_target_protein: String,
    mismatches: String,
    delta_mass_da: String,
    diff_positions: String,
    charge: String,
    precursor_mz: String,
    retention_time: String,
    scan_number: String,
    spectrum_file: String,
    protein_ids: String,
    q_value: String,
}

impl PsmRow {
    /// Convert a [`ClassifiedPsm`] into a flat row suitable for JSON serialisation.
    fn from_classified(cp: &ClassifiedPsm) -> Self {
        Self {
            peptide: cp.psm.peptide.clone(),
            group: cp.group.to_string(),
            level: cp.level.as_str().to_string(),
            best_target_peptide: cp.best_target_peptide.clone().unwrap_or_default(),
            best_target_protein: cp.best_target_protein.clone().unwrap_or_default(),
            mismatches: cp.mismatches.map(|m| m.to_string()).unwrap_or_default(),
            delta_mass_da: cp
                .delta_mass_da
                .map(|d| format!("{:.4}", d))
                .unwrap_or_default(),
            diff_positions: cp.diff_positions.clone().unwrap_or_default(),
            charge: cp.psm.charge.map(|c| c.to_string()).unwrap_or_default(),
            precursor_mz: cp
                .psm
                .precursor_mz
                .map(|m| format!("{:.4}", m))
                .unwrap_or_default(),
            retention_time: cp
                .psm
                .retention_time
                .map(|r| format!("{:.2}", r))
                .unwrap_or_default(),
            scan_number: cp
                .psm
                .scan_number
                .map(|s| s.to_string())
                .unwrap_or_default(),
            spectrum_file: cp.psm.spectrum_file.clone().unwrap_or_default(),
            protein_ids: cp.psm.protein_ids.clone(),
            q_value: cp
                .psm
                .q_value
                .map(|q| format!("{:.6}", q))
                .unwrap_or_default(),
        }
    }
}

/// Render an interactive HTML report and write it to `output_path`.
///
/// Only trap PSMs are included in the detail table; the summary section
/// covers all groups (target, trap, ambiguous).
///
/// # Errors
///
/// Returns [`EntrapmentError::ReportError`] if JSON serialisation or file I/O fails.
pub fn render_report(
    summary: &EntrapmentSummary,
    classified: &[ClassifiedPsm],
    output_path: &Path,
) -> Result<(), EntrapmentError> {
    // Only include trap PSMs in the detail table
    let psm_rows: Vec<PsmRow> = classified
        .iter()
        .filter(|cp| cp.group == PsmGroup::Trap)
        .map(PsmRow::from_classified)
        .collect();

    info!(
        total = classified.len(),
        trap_rows = psm_rows.len(),
        output = %output_path.display(),
        "rendering entrapment HTML report"
    );

    let report_data = ReportData {
        summary: summary.clone(),
        psms: psm_rows,
    };

    let json = serde_json::to_string(&report_data).map_err(|e| EntrapmentError::ReportError {
        detail: format!("JSON serialization failed: {e}"),
    })?;

    // Escape "</script>" sequences in JSON to prevent HTML injection (I5)
    let safe_json = json.replace("</", "<\\/");
    let html = TEMPLATE.replace(DATA_PLACEHOLDER, &safe_json);

    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| EntrapmentError::ReportError {
            detail: format!("failed to create directory {}: {e}", parent.display()),
        })?;
    }

    std::fs::write(output_path, html).map_err(|e| EntrapmentError::ReportError {
        detail: format!("failed to write report to {}: {e}", output_path.display()),
    })?;

    info!(path = %output_path.display(), "entrapment report written");
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{
        ClassifiedPsm, DiscriminabilityLevel, EntrapmentSummary, LevelCounts, PsmGroup,
        RazorFamily, SubstitutionType, UnifiedPsm,
    };

    fn make_psm(peptide: &str, group: PsmGroup, level: DiscriminabilityLevel) -> ClassifiedPsm {
        ClassifiedPsm {
            psm: UnifiedPsm {
                peptide: peptide.to_string(),
                charge: Some(2),
                precursor_mz: Some(500.25),
                retention_time: Some(12.5),
                scan_number: Some(100),
                spectrum_file: Some("test.raw".to_string()),
                protein_ids: "sp|P12345|TEST_HUMAN".to_string(),
                q_value: Some(0.01),
            },
            group,
            level,
            best_target_peptide: Some("PEPTIDER".to_string()),
            best_target_protein: Some("sp|Q99999|TARG_HUMAN".to_string()),
            mismatches: Some(1),
            delta_mass_da: Some(0.0364),
            diff_positions: Some("[3:A->G]".to_string()),
            substitution_type: SubstitutionType::None,
            edit_distance: None,
            alignment_detail: None,
        }
    }

    fn make_summary() -> EntrapmentSummary {
        EntrapmentSummary {
            total_psms: 100,
            target_psms: 80,
            trap_psms: 15,
            ambiguous_psms: 5,
            level_counts: LevelCounts {
                l0: 3,
                l1: 2,
                l2: 4,
                l3: 3,
                l4: 3,
            },
            top_razor_families: vec![RazorFamily {
                family: "EF1A1".to_string(),
                count: 3,
                example_peptide: "PEPTIDEK".to_string(),
                example_trap_protein: "sp|P12345|EF1A1_BOVIN".to_string(),
                example_target_protein: "sp|Q99999|EF1A1_HUMAN".to_string(),
            }],
        }
    }

    #[test]
    fn test_psm_row_from_classified() {
        let cp = make_psm("TESTPEPTIDE", PsmGroup::Trap, DiscriminabilityLevel::L2);
        let row = PsmRow::from_classified(&cp);
        assert_eq!(row.peptide, "TESTPEPTIDE");
        assert_eq!(row.group, "trap");
        assert_eq!(row.level, "L2");
        assert_eq!(row.best_target_peptide, "PEPTIDER");
        assert_eq!(row.mismatches, "1");
        assert_eq!(row.delta_mass_da, "0.0364");
        assert_eq!(row.charge, "2");
    }

    #[test]
    fn test_psm_row_empty_optionals() {
        let mut cp = make_psm("SEQ", PsmGroup::Trap, DiscriminabilityLevel::L4);
        cp.best_target_peptide = None;
        cp.best_target_protein = None;
        cp.mismatches = None;
        cp.delta_mass_da = None;
        cp.diff_positions = None;
        cp.psm.charge = None;
        cp.psm.precursor_mz = None;
        let row = PsmRow::from_classified(&cp);
        assert_eq!(row.best_target_peptide, "");
        assert_eq!(row.mismatches, "");
        assert_eq!(row.delta_mass_da, "");
        assert_eq!(row.charge, "");
    }

    #[test]
    fn test_render_report_creates_file() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("report.html");
        let summary = make_summary();
        let classified = vec![
            make_psm("AAAK", PsmGroup::Target, DiscriminabilityLevel::L4),
            make_psm("BBBK", PsmGroup::Trap, DiscriminabilityLevel::L0),
            make_psm("CCCK", PsmGroup::Trap, DiscriminabilityLevel::L2),
        ];
        render_report(&summary, &classified, &out).expect("render_report");
        let content = std::fs::read_to_string(&out).expect("read");
        // Should contain the template structure
        assert!(content.contains("Entrapment Analysis Report"));
        // Should contain serialised PSM data
        assert!(content.contains("BBBK"));
        assert!(content.contains("CCCK"));
        // Target PSMs should NOT be in the detail table data
        // (the JSON psms array only has trap rows)
        // Check that "AAAK" only appears if it slipped through - it shouldn't be in psms
        assert!(content.contains("Plotly"));
    }

    #[test]
    fn test_render_report_only_trap_psms_in_table() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("report2.html");
        let summary = make_summary();
        let classified = vec![
            make_psm("TARGET_ONLY", PsmGroup::Target, DiscriminabilityLevel::L4),
            make_psm("TRAP_ONE", PsmGroup::Trap, DiscriminabilityLevel::L1),
        ];
        render_report(&summary, &classified, &out).expect("render_report");
        let content = std::fs::read_to_string(&out).expect("read");
        assert!(content.contains("TRAP_ONE"));
        // TARGET_ONLY should not appear in the psms JSON array
        // We can verify by checking it doesn't appear as a peptide value
        assert!(!content.contains("TARGET_ONLY"));
    }

    #[test]
    fn test_template_contains_placeholder_structure() {
        // Verify the template can be loaded and has the expected placeholder
        assert!(TEMPLATE.contains(DATA_PLACEHOLDER));
        assert!(TEMPLATE.contains("Entrapment Analysis Report"));
        assert!(TEMPLATE.contains("plotly"));
    }

    #[test]
    fn test_render_report_creates_parent_dirs() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("deep").join("nested").join("report.html");
        let summary = make_summary();
        render_report(&summary, &[], &out).expect("render_report");
        assert!(out.exists());
    }
}
