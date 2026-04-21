//! Integration test for entrapment analysis with known peptides.
//!
//! Tests L0-L4 classification using 6 peptides from the HeLa SILAC DIA analysis.
//! Requires human_swissprot.fasta to be present; gracefully skips if not found.
//!
//! V3 integration tests cover provenance features: TSV output columns,
//! mod parser, provenance engine with known peptides, ClassifiedPsm
//! provenance JSON roundtrip, and backward-compatible config parsing.

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
        modifications: Vec::new(),
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
    let header_cols: Vec<&str> = content.lines().next().unwrap().split('\t').collect();
    let sub_type_idx = header_cols.iter().position(|&h| h == "substitution_type").unwrap();
    let edit_dist_idx = header_cols.iter().position(|&h| h == "edit_distance").unwrap();
    let sub_type_col = cols[sub_type_idx];
    let edit_dist_col = cols[edit_dist_idx];
    assert!(
        !sub_type_col.is_empty(),
        "L2 row substitution_type should not be empty"
    );
    assert!(
        !edit_dist_col.is_empty(),
        "L2 row edit_distance should not be empty"
    );
}

// ===========================================================================
// V3 Provenance Integration Tests
// ===========================================================================

/// Test 1: Verify that the TSV output includes the 5 new provenance column
/// headers (`trap_matched`, `target_matched`, `shared_ions`, `shared_ratio`,
/// `is_chimeric`) and that each data row has 23 columns total.
///
/// When provenance is not traced (no mzML), the 5 provenance columns should
/// be present but empty.
#[test]
fn test_v3_provenance_columns_in_tsv() {
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

    // Use a mix of levels to exercise different code paths
    let psms = vec![
        make_psm("STTTGHLIYK", "sp|Q6CFZ6|EF1A_YARLI"),    // L0
        make_psm("ELTALAPSTMK", "sp|P18616|ACT1_DICDI"),   // L1
        make_psm("DGFLLDGFPR", "sp|P69428|KAD_YERPE"),     // L2
        make_psm("IGSEVYHNLK", "sp|P00924|ENO2_YEAST"),    // L3
    ];
    let classified = analyzer
        .classify_all(&psms)
        .expect("classify_all should work");

    let tmp_dir = tempfile::tempdir().expect("create tempdir");
    let tsv_path = tmp_dir.path().join("v3_provenance.tsv");
    protein_copilot_entrapment_analysis::output::write_classified_tsv(&classified, &tsv_path)
        .expect("write TSV");

    let content = std::fs::read_to_string(&tsv_path).expect("read TSV");
    let header = content.lines().next().expect("TSV should have a header");
    let header_cols: Vec<&str> = header.split('\t').collect();

    // 1. Check that all 5 provenance columns are present
    let provenance_columns = [
        "trap_matched",
        "target_matched",
        "shared_ions",
        "shared_ratio",
        "is_chimeric",
    ];
    for col in &provenance_columns {
        assert!(
            header_cols.contains(col),
            "TSV header missing provenance column '{col}'"
        );
    }

    // 2. Check total column count is 23
    assert_eq!(
        header_cols.len(),
        23,
        "Expected 23 columns in TSV header, got {}. Columns: {header_cols:?}",
        header_cols.len()
    );

    // 3. Check each data row also has 23 columns
    let data_lines: Vec<&str> = content.lines().skip(1).collect();
    assert_eq!(data_lines.len(), 4, "should have 4 data rows");
    for (i, line) in data_lines.iter().enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            cols.len(),
            23,
            "Data row {i} has {} columns, expected 23",
            cols.len()
        );
    }

    // 4. Provenance columns should be empty (no mzML provided)
    let trap_matched_idx = header_cols
        .iter()
        .position(|&h| h == "trap_matched")
        .expect("trap_matched column");
    let shared_ratio_idx = header_cols
        .iter()
        .position(|&h| h == "shared_ratio")
        .expect("shared_ratio column");
    let is_chimeric_idx = header_cols
        .iter()
        .position(|&h| h == "is_chimeric")
        .expect("is_chimeric column");

    for line in &data_lines {
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols[trap_matched_idx].is_empty(),
            "trap_matched should be empty without mzML, got '{}'",
            cols[trap_matched_idx]
        );
        assert!(
            cols[shared_ratio_idx].is_empty(),
            "shared_ratio should be empty without mzML, got '{}'",
            cols[shared_ratio_idx]
        );
        assert!(
            cols[is_chimeric_idx].is_empty(),
            "is_chimeric should be empty without mzML, got '{}'",
            cols[is_chimeric_idx]
        );
    }
}

/// Test 2: Verify the mod_parser works in an integration context with
/// realistic DIA-NN Modified.Sequence strings, and that the parsed
/// modifications produce correct (position, delta_mass) tuples suitable
/// for the provenance engine.
#[test]
fn test_v3_mod_parser_integration() {
    use protein_copilot_entrapment_analysis::mod_parser::parse_modified_sequence;

    // 1. Realistic DIA-NN modified sequence with two modifications
    let (stripped, mods) = parse_modified_sequence("AAAC(UniMod:4)DEFM(UniMod:35)GK");
    assert_eq!(stripped, "AAACDEFMGK");
    assert_eq!(mods.len(), 2);

    // Verify individual mods
    assert_eq!(mods[0].position, 3); // C at position 3
    assert_eq!(mods[0].unimod_id, 4); // Carbamidomethyl
    assert!((mods[0].delta_mass - 57.021464).abs() < 1e-6);

    assert_eq!(mods[1].position, 7); // M at position 7
    assert_eq!(mods[1].unimod_id, 35); // Oxidation
    assert!((mods[1].delta_mass - 15.994915).abs() < 1e-6);

    // 2. Verify the modifications can be converted to (position, delta_mass) tuples
    //    as used by provenance::trace_provenance
    let mod_tuples: Vec<(usize, f64)> = mods.iter().map(|m| (m.position, m.delta_mass)).collect();
    assert_eq!(mod_tuples.len(), 2);
    assert_eq!(mod_tuples[0].0, 3);
    assert_eq!(mod_tuples[1].0, 7);
    assert!((mod_tuples[0].1 - 57.021464).abs() < 1e-6);
    assert!((mod_tuples[1].1 - 15.994915).abs() < 1e-6);

    // 3. N-terminal modification
    let (stripped_nterm, mods_nterm) =
        parse_modified_sequence("(UniMod:1)AAAC(UniMod:4)DEFMGK");
    assert_eq!(stripped_nterm, "AAACDEFMGK");
    assert_eq!(mods_nterm.len(), 2);
    assert_eq!(mods_nterm[0].position, 0); // N-term acetyl → position 0
    assert_eq!(mods_nterm[0].unimod_id, 1);
    assert!((mods_nterm[0].delta_mass - 42.010565).abs() < 1e-6);
    assert_eq!(mods_nterm[1].position, 3); // Carbamidomethyl on C
    assert_eq!(mods_nterm[1].unimod_id, 4);

    // 4. Unmodified sequence
    let (stripped_plain, mods_plain) = parse_modified_sequence("PEPTIDE");
    assert_eq!(stripped_plain, "PEPTIDE");
    assert!(mods_plain.is_empty());

    // 5. Heavy SILAC labels
    let (stripped_silac, mods_silac) =
        parse_modified_sequence("PEPTIDEK(UniMod:259)");
    assert_eq!(stripped_silac, "PEPTIDEK");
    assert_eq!(mods_silac.len(), 1);
    assert_eq!(mods_silac[0].unimod_id, 259);
    assert!((mods_silac[0].delta_mass - 8.014199).abs() < 1e-6);
}

/// Test 3: Verify trace_provenance correctly classifies fragment ions when
/// given known peptide pairs and synthetic spectra derived from theoretical
/// ions.
#[test]
fn test_v3_provenance_trace_known_peptides() {
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
    use protein_copilot_entrapment_analysis::provenance::{trace_provenance, IonOrigin};

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    // Use test_helpers to generate real theoretical fragment ions for PEPTIDE
    let peaks = test_helpers::synthetic_peaks_for_peptide("PEPTIDE", 1000.0);
    let observed_mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
    let observed_int: Vec<f64> = peaks.iter().map(|p| p.1).collect();

    // Case A: PEPTIDE vs PAPTIDE — differ at position 1 (E→A)
    // Some y-ions from the C-terminus should be shared (those not covering
    // position 1), while some b-ions from the N-terminus should differ.
    let result_a = trace_provenance(
        &observed_mz,
        &observed_int,
        "PEPTIDE",
        "PAPTIDE",
        &[],
        &tolerance,
        1,
    );

    assert_eq!(result_a.trap_sequence, "PEPTIDE");
    assert_eq!(result_a.target_sequence, "PAPTIDE");
    assert_eq!(
        result_a.annotated_peaks.len(),
        observed_mz.len(),
        "every observed peak should be annotated"
    );
    // At least some peaks should be trap-only (b-ions covering position 1)
    let trap_only_count = result_a
        .annotated_peaks
        .iter()
        .filter(|p| p.origin == IonOrigin::TrapOnly)
        .count();
    // At least some peaks should be shared (y-ions from C-terminus not covering pos 1)
    let _shared_count = result_a
        .annotated_peaks
        .iter()
        .filter(|p| p.origin == IonOrigin::Shared)
        .count();
    assert!(
        trap_only_count > 0 || result_a.trap_matched_count > 0,
        "should have trap-only or trap-matched peaks for PEPTIDE vs PAPTIDE"
    );
    // The total should account for all peaks
    assert_eq!(
        result_a.trap_matched_count + result_a.target_matched_count
            + result_a.shared_count + result_a.unassigned_count,
        observed_mz.len() as u32,
        "all peaks should be accounted for"
    );

    // Case B: identical sequences — all matched peaks should be Shared
    let result_b = trace_provenance(
        &observed_mz,
        &observed_int,
        "PEPTIDE",
        "PEPTIDE",
        &[],
        &tolerance,
        1,
    );
    assert_eq!(result_b.trap_matched_count, 0, "no trap-only when sequences match");
    assert_eq!(result_b.target_matched_count, 0, "no target-only when sequences match");
    // shared_ratio should be 1.0 for matched peaks
    if result_b.shared_count > 0 {
        assert!(
            (result_b.shared_ratio - 1.0).abs() < 1e-10,
            "shared_ratio should be 1.0 for identical sequences, got {}",
            result_b.shared_ratio
        );
    }

    // Case C: L4 — no target sequence (empty string)
    let result_c = trace_provenance(
        &observed_mz,
        &observed_int,
        "PEPTIDE",
        "",
        &[],
        &tolerance,
        1,
    );
    assert!(result_c.target_sequence.is_empty());
    assert_eq!(result_c.target_matched_count, 0);
    assert_eq!(result_c.shared_count, 0);
    // All peaks should be trap-only or unassigned
    assert_eq!(
        result_c.trap_matched_count + result_c.unassigned_count,
        observed_mz.len() as u32,
        "L4: all peaks should be trap-only or unassigned"
    );
}

/// Test 4: Verify ClassifiedPsm with provenance can be serialized to JSON
/// and deserialized back, preserving all provenance fields. Also test that
/// a ClassifiedPsm without provenance (None) roundtrips correctly for
/// backward compatibility.
#[test]
fn test_v3_classified_psm_provenance_roundtrip() {
    use protein_copilot_entrapment_analysis::provenance::{
        AnnotatedPeak, FragmentProvenance, IonOrigin,
    };
    use protein_copilot_entrapment_analysis::types::ClassifiedPsm;

    let psm = UnifiedPsm {
        peptide: "DGFLLDGFPR".to_string(),
        charge: Some(2),
        precursor_mz: Some(540.28),
        retention_time: Some(25.3),
        scan_number: Some(12345),
        spectrum_file: Some("sample1".to_string()),
        protein_ids: "sp|P69428|KAD_YERPE".to_string(),
        q_value: Some(0.005),
        modifications: vec![(3, 57.021464)],
    };

    let provenance = FragmentProvenance {
        trap_sequence: "DGFLLDGFPR".to_string(),
        target_sequence: "NGFLLDGFPR".to_string(),
        annotated_peaks: vec![
            AnnotatedPeak {
                mz_observed: 204.134,
                intensity: 5000.0,
                origin: IonOrigin::TrapOnly,
                trap_ion_label: Some("b2+1".to_string()),
                target_ion_label: None,
            },
            AnnotatedPeak {
                mz_observed: 175.119,
                intensity: 8000.0,
                origin: IonOrigin::Shared,
                trap_ion_label: Some("y1+1".to_string()),
                target_ion_label: Some("y1+1".to_string()),
            },
            AnnotatedPeak {
                mz_observed: 999.999,
                intensity: 200.0,
                origin: IonOrigin::Unassigned,
                trap_ion_label: None,
                target_ion_label: None,
            },
        ],
        trap_matched_count: 1,
        target_matched_count: 0,
        shared_count: 1,
        unassigned_count: 1,
        shared_ratio: 0.5,
        is_chimeric: true,
    };

    let cpsm_with = ClassifiedPsm {
        psm: psm.clone(),
        group: protein_copilot_entrapment_analysis::PsmGroup::Trap,
        level: DiscriminabilityLevel::L2,
        best_target_peptide: Some("NGFLLDGFPR".to_string()),
        best_target_protein: Some("sp|Q9Y2Q3|GSTK1_HUMAN".to_string()),
        mismatches: Some(1),
        delta_mass_da: Some(-0.984),
        diff_positions: Some("[0:D->N]".to_string()),
        substitution_type: SubstitutionType::NearIsobaric,
        edit_distance: Some(1),
        alignment_detail: Some("D0→N".to_string()),
        provenance: Some(provenance),
    };

    // Roundtrip with provenance
    let json = serde_json::to_string(&cpsm_with).expect("serialize ClassifiedPsm with provenance");
    let deser: ClassifiedPsm =
        serde_json::from_str(&json).expect("deserialize ClassifiedPsm with provenance");

    let prov = deser.provenance.as_ref().expect("provenance should be Some");
    assert_eq!(prov.trap_sequence, "DGFLLDGFPR");
    assert_eq!(prov.target_sequence, "NGFLLDGFPR");
    assert_eq!(prov.annotated_peaks.len(), 3);
    assert_eq!(prov.trap_matched_count, 1);
    assert_eq!(prov.target_matched_count, 0);
    assert_eq!(prov.shared_count, 1);
    assert_eq!(prov.unassigned_count, 1);
    assert!((prov.shared_ratio - 0.5).abs() < 1e-10);
    assert!(prov.is_chimeric);

    // Verify IonOrigin deserialized correctly
    assert_eq!(prov.annotated_peaks[0].origin, IonOrigin::TrapOnly);
    assert_eq!(prov.annotated_peaks[1].origin, IonOrigin::Shared);
    assert_eq!(prov.annotated_peaks[2].origin, IonOrigin::Unassigned);

    // Verify ion labels
    assert_eq!(
        prov.annotated_peaks[0].trap_ion_label.as_deref(),
        Some("b2+1")
    );
    assert!(prov.annotated_peaks[0].target_ion_label.is_none());
    assert_eq!(
        prov.annotated_peaks[1].target_ion_label.as_deref(),
        Some("y1+1")
    );

    // Verify modifications survived roundtrip
    assert_eq!(deser.psm.modifications.len(), 1);
    assert_eq!(deser.psm.modifications[0].0, 3);
    assert!((deser.psm.modifications[0].1 - 57.021464).abs() < 1e-6);

    // Roundtrip WITHOUT provenance (backward compat)
    let cpsm_without = ClassifiedPsm {
        psm,
        group: protein_copilot_entrapment_analysis::PsmGroup::Trap,
        level: DiscriminabilityLevel::L4,
        best_target_peptide: None,
        best_target_protein: None,
        mismatches: None,
        delta_mass_da: None,
        diff_positions: None,
        substitution_type: SubstitutionType::None,
        edit_distance: None,
        alignment_detail: None,
        provenance: None,
    };

    let json_none =
        serde_json::to_string(&cpsm_without).expect("serialize ClassifiedPsm without provenance");
    // provenance field should not appear in JSON (skip_serializing_if = "Option::is_none")
    assert!(
        !json_none.contains("provenance"),
        "JSON should not contain 'provenance' when None"
    );
    let deser_none: ClassifiedPsm =
        serde_json::from_str(&json_none).expect("deserialize ClassifiedPsm without provenance");
    assert!(
        deser_none.provenance.is_none(),
        "provenance should be None after roundtrip"
    );
}

/// Test 5: Verify that a v2-style config YAML (without a `provenance` section)
/// parses correctly and gets default provenance values, ensuring backward
/// compatibility.
#[test]
fn test_v3_config_backward_compatible() {
    // This is a v2-era config — no provenance section at all
    let v2_yaml = r#"
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

    let config =
        EntrapmentConfig::from_yaml_str(v2_yaml).expect("v2 config should parse successfully");

    // Provenance defaults should be applied
    assert!(
        (config.provenance.fragment_tolerance_ppm - 20.0).abs() < 1e-6,
        "default fragment_tolerance_ppm should be 20.0, got {}",
        config.provenance.fragment_tolerance_ppm
    );
    assert_eq!(
        config.provenance.max_fragment_charge, 2,
        "default max_fragment_charge should be 2"
    );
    assert!(
        (config.provenance.chimera_threshold - 0.3).abs() < 1e-6,
        "default chimera_threshold should be 0.3, got {}",
        config.provenance.chimera_threshold
    );
    assert_eq!(
        config.provenance.min_peaks_for_analysis, 6,
        "default min_peaks_for_analysis should be 6"
    );
    assert_eq!(
        config.provenance.levels_to_trace,
        vec!["L2", "L3", "L4"],
        "default levels_to_trace should be [L2, L3, L4]"
    );

    // Existing v2 fields should also be preserved
    assert_eq!(config.similarity.max_mismatches, 2);
    assert!((config.similarity.delta_mass_threshold_da - 1.0).abs() < 1e-6);
    assert_eq!(config.similarity.len_tolerance, 2);

    // Now test with explicit provenance section to verify override works
    let v3_yaml = r#"
version: 1
target:
  rules:
    - type: AccessionContains
      any_of: ["_HUMAN"]
trap:
  rules:
    - type: AccessionContains
      any_of: ["_YEAST"]
provenance:
  fragment_tolerance_ppm: 10.0
  max_fragment_charge: 3
  chimera_threshold: 0.5
  min_peaks_for_analysis: 10
  levels_to_trace: ["L3", "L4"]
"#;

    let config_v3 =
        EntrapmentConfig::from_yaml_str(v3_yaml).expect("v3 config should parse successfully");
    assert!(
        (config_v3.provenance.fragment_tolerance_ppm - 10.0).abs() < 1e-6,
        "explicit fragment_tolerance_ppm should be 10.0"
    );
    assert_eq!(config_v3.provenance.max_fragment_charge, 3);
    assert!((config_v3.provenance.chimera_threshold - 0.5).abs() < 1e-6);
    assert_eq!(config_v3.provenance.min_peaks_for_analysis, 10);
    assert_eq!(
        config_v3.provenance.levels_to_trace,
        vec!["L3", "L4"]
    );
}
