//! Integration tests for entrapment v2: edit distance + substitution type.

use std::io::Write;

use protein_copilot_entrapment_analysis::config::EntrapmentConfig;
use protein_copilot_entrapment_analysis::{
    DiscriminabilityLevel, EntrapmentAnalyzer, PsmGroup, SubstitutionType, UnifiedPsm,
};

fn make_psm(peptide: &str, protein: &str) -> UnifiedPsm {
    UnifiedPsm {
        peptide: peptide.to_string(),
        charge: Some(2),
        precursor_mz: None,
        retention_time: None,
        scan_number: None,
        spectrum_file: None,
        protein_ids: protein.to_string(),
        q_value: Some(0.01),
    }
}

fn make_config_yaml() -> String {
    r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST", "_ECOLI"]
similarity:
  max_mismatches: 2
  delta_mass_threshold_da: 1.0
  len_tolerance: 2
"#
    .to_string()
}

#[test]
fn test_cross_length_indel_gets_l2_not_l4() {
    // Target has "PEPGGDEK" (8 aa), trap PSM has "PEPNDEK" (7 aa)
    // N↔GG isobaric dipeptide, delta_mass ≈ 0
    // v1 would classify as L4 (different length), v2 should classify as L2
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(
        fasta_file,
        ">sp|P001|TEST_HUMAN Test protein\nPEPGGDEKLASTR\n"
    )
    .unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("PEPNDEK", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.group, PsmGroup::Trap);
    // This is the key v2 assertion: should NOT be L4 anymore
    assert!(
        result.level == DiscriminabilityLevel::L2 || result.level == DiscriminabilityLevel::L3,
        "cross-length indel should be L2 or L3, got {:?}",
        result.level
    );
    assert!(result.edit_distance.is_some());
}

#[test]
fn test_qk_substitution_detected() {
    // Target: "PEPTKDEK", trap: "PEPTQDEK" — Q↔K substitution
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(fasta_file, ">sp|P001|TEST_HUMAN Test\nPEPTKDEKLASTR\n").unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("PEPTQDEK", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.level, DiscriminabilityLevel::L2);
    assert_eq!(result.substitution_type, SubstitutionType::QKSubstitution);
}

#[test]
fn test_backward_compatible_same_length_matching() {
    // Same-length D→N substitution should still work exactly as v1
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(fasta_file, ">sp|P001|TEST_HUMAN Test\nNGFLLDGFPR\n").unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("DGFLLDGFPR", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.level, DiscriminabilityLevel::L2);
    assert_eq!(result.mismatches, Some(1));
    assert!(result.delta_mass_da.unwrap().abs() < 1.0);
}

#[test]
fn test_true_trap_still_l4() {
    // Completely unrelated peptide should still be L4
    let mut fasta_file = tempfile::NamedTempFile::new().unwrap();
    write!(fasta_file, ">sp|P001|TEST_HUMAN Test\nAAAAAAAAR\n").unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config, fasta_file.path()).unwrap();

    let psm = make_psm("WWWWWWWWR", "sp|Q001|TRAP_YEAST");
    let result = analyzer.classify(&psm).unwrap();

    assert_eq!(result.level, DiscriminabilityLevel::L4);
    assert_eq!(result.substitution_type, SubstitutionType::None);
}
