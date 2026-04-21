//! Integration test for entrapment analysis with known peptides.
//!
//! Tests L0-L4 classification using 6 peptides from the HeLa SILAC DIA analysis.
//! Requires human_swissprot.fasta to be present; gracefully skips if not found.

use std::path::Path;

use protein_copilot_entrapment_analysis::{
    config::EntrapmentConfig, DiscriminabilityLevel, EntrapmentAnalyzer, SubstitutionType,
    UnifiedPsm,
};

/// Relative path from the workspace root to the human SwissProt FASTA file.
const FASTA_RELATIVE: &str = ".proteincopilot/databases/human_swissprot.fasta";

/// Resolve the FASTA path relative to the workspace root.
///
/// `CARGO_MANIFEST_DIR` points to `crates/integration-tests/`; the workspace
/// root is two levels up.
fn fasta_path() -> std::path::PathBuf {
    let manifest_dir =
        std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR should be set by cargo");
    Path::new(&manifest_dir).join("../..").join(FASTA_RELATIVE)
}

fn make_psm(peptide: &str, protein_ids: &str) -> UnifiedPsm {
    UnifiedPsm {
        peptide: peptide.to_string(),
        charge: Some(2),
        precursor_mz: None,
        retention_time: None,
        scan_number: None,
        spectrum_file: None,
        protein_ids: protein_ids.to_string(),
        q_value: Some(0.01),
    }
}

fn test_config_yaml() -> &'static str {
    r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI", "_DICDI", "_YARLI", "_YERPE"]
conflict_resolution: prefer_target
similarity:
  max_mismatches: 2
  delta_mz_threshold_da: 1.0
"#
}

#[test]
fn test_known_peptide_classifications() {
    let fasta_path = fasta_path();
    let fasta = fasta_path.as_path();
    if !fasta.exists() {
        eprintln!(
            "SKIP: {} not found — run `download_database human_swissprot` first",
            fasta.display()
        );
        return;
    }

    let config = EntrapmentConfig::from_yaml_str(test_config_yaml()).expect("config should parse");

    let analyzer = EntrapmentAnalyzer::new(config, fasta).expect("analyzer should build");

    // Test cases: (peptide, protein_ids, expected_level, description)
    let cases = vec![
        (
            "STTTGHLIYK",
            "sp|Q6CFZ6|EF1A_YARLI",
            DiscriminabilityLevel::L0,
            "exact match → L0",
        ),
        (
            "GYSFTTTAER",
            "sp|P18616|ACT1_DICDI",
            DiscriminabilityLevel::L0,
            "exact match → L0",
        ),
        (
            "ELTALAPSTMK",
            "sp|P18616|ACT1_DICDI",
            DiscriminabilityLevel::L1,
            "L/I isomer → L1",
        ),
        (
            "HPFPGPGIAIR",
            "sp|P07926|GUAA_YEAST",
            DiscriminabilityLevel::L1,
            "L/I isomer → L1",
        ),
        (
            "DGFLLDGFPR",
            "sp|P69428|KAD_YERPE",
            DiscriminabilityLevel::L2,
            "1mm near-isobaric → L2",
        ),
        (
            "IGSEVYHNLK",
            "sp|P00924|ENO2_YEAST",
            DiscriminabilityLevel::L3,
            "1mm large Δm → L3",
        ),
    ];

    for (peptide, protein, expected_level, desc) in &cases {
        let psm = make_psm(peptide, protein);
        let result = analyzer
            .classify(&psm)
            .unwrap_or_else(|e| panic!("classify failed for {peptide}: {e}"));

        assert_eq!(
            result.level,
            *expected_level,
            "Peptide {peptide} ({desc}): expected {expected_level}, got {}. \
             Best match: {:?} in {:?}, mismatches={:?}, Δm={:?}",
            result.level,
            result.best_target_peptide,
            result.best_target_protein,
            result.mismatches,
            result.delta_mass_da,
        );
    }

    // Also verify summary statistics
    let psms: Vec<UnifiedPsm> = cases
        .iter()
        .map(|(pep, prot, _, _)| make_psm(pep, prot))
        .collect();
    let classified = analyzer
        .classify_all(&psms)
        .expect("classify_all should work");
    let summary = analyzer.summary(&classified);

    assert_eq!(summary.trap_psms, 6);
    assert_eq!(summary.level_counts.l0, 2);
    assert_eq!(summary.level_counts.l1, 2);
    assert_eq!(summary.level_counts.l2, 1);
    assert_eq!(summary.level_counts.l3, 1);
    assert_eq!(summary.level_counts.l4, 0);
}

/// V2 end-to-end test: verify substitution_type, edit_distance, alignment_detail
/// are correctly populated through the full pipeline, including TSV output.
#[test]
fn test_v2_fields_end_to_end() {
    let fasta_path = fasta_path();
    let fasta = fasta_path.as_path();
    if !fasta.exists() {
        eprintln!(
            "SKIP: {} not found — run `download_database human_swissprot` first",
            fasta.display()
        );
        return;
    }

    let config_yaml = r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI", "_DICDI", "_YARLI", "_YERPE"]
conflict_resolution: prefer_target
similarity:
  max_mismatches: 2
  delta_mass_threshold_da: 1.0
  len_tolerance: 2
"#;
    let config = EntrapmentConfig::from_yaml_str(config_yaml).expect("config should parse");
    let analyzer = EntrapmentAnalyzer::new(config, fasta).expect("analyzer should build");

    // L2 near-isobaric: DGFLLDGFPR (trap) → NGFLLDGFPR (human), D→N, Δm ≈ -0.984 Da
    let psm_l2 = make_psm("DGFLLDGFPR", "sp|P69428|KAD_YERPE");
    let r_l2 = analyzer.classify(&psm_l2).unwrap();
    assert_eq!(r_l2.level, DiscriminabilityLevel::L2);
    assert!(
        r_l2.edit_distance.is_some(),
        "L2 should have edit_distance, got None"
    );
    assert!(
        r_l2.alignment_detail.is_some(),
        "L2 should have alignment_detail, got None"
    );
    assert_ne!(
        r_l2.substitution_type,
        SubstitutionType::None,
        "L2 near-isobaric should have a non-None substitution_type"
    );

    // L0 exact match: fields should be None/SubstitutionType::None
    let psm_l0 = make_psm("STTTGHLIYK", "sp|Q6CFZ6|EF1A_YARLI");
    let r_l0 = analyzer.classify(&psm_l0).unwrap();
    assert_eq!(r_l0.level, DiscriminabilityLevel::L0);
    assert_eq!(r_l0.substitution_type, SubstitutionType::None);
    assert!(r_l0.edit_distance.is_none());
    assert!(r_l0.alignment_detail.is_none());

    // L3 large delta-mass: should have v2 fields populated
    let psm_l3 = make_psm("IGSEVYHNLK", "sp|P00924|ENO2_YEAST");
    let r_l3 = analyzer.classify(&psm_l3).unwrap();
    assert_eq!(r_l3.level, DiscriminabilityLevel::L3);
    assert!(r_l3.edit_distance.is_some());

    // TSV output: verify v2 columns present and populated
    let psms = vec![psm_l2, psm_l0, psm_l3];
    let classified = analyzer.classify_all(&psms).unwrap();

    let tmp_dir = tempfile::tempdir().unwrap();
    let tsv_path = tmp_dir.path().join("v2_e2e.tsv");
    protein_copilot_entrapment_analysis::output::write_classified_tsv(&classified, &tsv_path)
        .unwrap();

    let content = std::fs::read_to_string(&tsv_path).unwrap();
    let header = content.lines().next().unwrap();
    assert!(
        header.contains("substitution_type"),
        "TSV missing substitution_type column"
    );
    assert!(
        header.contains("edit_distance"),
        "TSV missing edit_distance column"
    );
    assert!(
        header.contains("alignment_detail"),
        "TSV missing alignment_detail column"
    );

    // Verify data rows have values (not all empty)
    let data_lines: Vec<&str> = content.lines().skip(1).collect();
    assert_eq!(data_lines.len(), 3, "should have 3 data rows");

    // The L2 row should have non-empty substitution_type and edit_distance
    let l2_row = data_lines
        .iter()
        .find(|l| l.contains("DGFLLDGFPR"))
        .expect("should find DGFLLDGFPR row");
    let cols: Vec<&str> = l2_row.split('\t').collect();
    // substitution_type is second-to-last-2, edit_distance is second-to-last-1
    let sub_type_col = cols.iter().rev().nth(2).unwrap();
    let edit_dist_col = cols.iter().rev().nth(1).unwrap();
    assert!(
        !sub_type_col.is_empty(),
        "L2 row substitution_type should not be empty"
    );
    assert!(
        !edit_dist_col.is_empty(),
        "L2 row edit_distance should not be empty"
    );
}
