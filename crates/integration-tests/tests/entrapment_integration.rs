//! Integration test for entrapment analysis with known peptides.
//!
//! Tests L0-L4 classification using 6 peptides from the HeLa SILAC DIA analysis.
//! Requires human_swissprot.fasta to be present; gracefully skips if not found.

use std::path::Path;

use protein_copilot_entrapment_analysis::{
    config::EntrapmentConfig, DiscriminabilityLevel, EntrapmentAnalyzer, UnifiedPsm,
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
