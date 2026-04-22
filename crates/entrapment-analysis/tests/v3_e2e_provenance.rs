//! End-to-end integration test for entrapment v3: provenance tracing.
//!
//! Creates synthetic FASTA, config, PSM, and MGF data with known theoretical
//! fragment ions, then runs the full classification → provenance pipeline and
//! verifies classification levels, provenance counts, chimera flag, and output.

use protein_copilot_entrapment_analysis::config::EntrapmentConfig;
use protein_copilot_entrapment_analysis::output::{write_classified_tsv, write_run_metadata, RunMetadata};
use protein_copilot_entrapment_analysis::provenance::{trace_provenance, IonOrigin};
use protein_copilot_entrapment_analysis::report;
use protein_copilot_entrapment_analysis::{
    trace_provenance_batch, DiscriminabilityLevel, EntrapmentAnalyzer, PsmGroup, UnifiedPsm,
};

use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

// ---------------------------------------------------------------------------
// Constants: amino acid masses (duplicated here for test-data generation)
// ---------------------------------------------------------------------------

const PROTON: f64 = 1.007276;
const WATER: f64 = 18.010565;

fn aa_mass(aa: char) -> f64 {
    match aa {
        'G' => 57.02146,
        'A' => 71.03711,
        'V' => 99.06841,
        'L' => 113.08406,
        'I' => 113.08406,
        'P' => 97.05276,
        'F' => 147.06841,
        'W' => 186.07931,
        'M' => 131.04049,
        'S' => 87.03203,
        'T' => 101.04768,
        'C' => 103.00919,
        'Y' => 163.06333,
        'H' => 137.05891,
        'D' => 115.02694,
        'E' => 129.04259,
        'N' => 114.04293,
        'Q' => 128.05858,
        'K' => 128.09496,
        'R' => 156.10111,
        _ => 0.0,
    }
}

/// Compute singly-charged b-ion m/z values for a peptide sequence.
fn b_ions(seq: &str) -> Vec<f64> {
    let residues: Vec<char> = seq.chars().collect();
    let n = residues.len();
    let mut ions = Vec::new();
    let mut cumulative = 0.0;
    for &r in residues.iter().take(n - 1) {
        cumulative += aa_mass(r);
        ions.push(cumulative + PROTON);
    }
    ions
}

/// Compute singly-charged y-ion m/z values for a peptide sequence.
fn y_ions(seq: &str) -> Vec<f64> {
    let residues: Vec<char> = seq.chars().collect();
    let n = residues.len();
    let mut ions = Vec::new();
    let mut cumulative = WATER;
    for i in (1..n).rev() {
        cumulative += aa_mass(residues[i]);
        ions.push(cumulative + PROTON);
    }
    ions
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn make_psm(
    peptide: &str,
    protein_ids: &str,
    scan_number: Option<u32>,
    spectrum_file: Option<&str>,
) -> UnifiedPsm {
    UnifiedPsm {
        peptide: peptide.to_string(),
        charge: Some(2),
        precursor_mz: None,
        retention_time: Some(300.0),
        scan_number,
        spectrum_file: spectrum_file.map(|s| s.to_string()),
        protein_ids: protein_ids.to_string(),
        q_value: Some(0.005),
        modifications: Vec::new(),
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
unmatched: trap
similarity:
  max_mismatches: 2
  delta_mass_threshold_da: 1.0
  len_tolerance: 2
provenance:
  levels_to_trace: ["L2", "L3", "L4"]
  fragment_tolerance_ppm: 20.0
  max_fragment_charge: 1
  chimera_threshold: 0.3
  min_peaks_for_analysis: 3
"#
    .to_string()
}

/// Build a minimal MGF with synthetic spectra whose peaks correspond to
/// known theoretical b/y-ions.
fn build_test_mgf() -> String {
    // Trap peptide: DGFLLDGFPR, Target peptide: NGFLLDGFPR
    // Difference: position 0, D→N (delta ≈ 0.984 Da, near-isobaric)
    //
    // b-ions differ at all positions (cumulative includes position 0).
    // y-ions y1..y9 are SHARED (suffix GFLLDGFPR is identical).

    let trap_b = b_ions("DGFLLDGFPR");
    let target_b = b_ions("NGFLLDGFPR");
    let shared_y = y_ions("DGFLLDGFPR"); // same as y_ions("NGFLLDGFPR")

    // Scan 1: spectrum for trap PSM "DGFLLDGFPR" (L2 classification)
    //   Peaks: trap b2, trap b4 (TrapOnly), target b3 (TargetOnly),
    //          shared y2, y4, y6 (Shared), noise 500.0, 800.0 (Unassigned)
    let scan1_peaks = vec![
        (trap_b[1], 600.0),     // trap b2 → TrapOnly
        (target_b[2], 500.0),   // target b3 → TargetOnly
        (shared_y[1], 700.0),   // y2 → Shared
        (trap_b[3], 800.0),     // trap b4 → TrapOnly
        (shared_y[3], 750.0),   // y4 → Shared
        (500.0, 300.0),         // noise → Unassigned
        (shared_y[5], 650.0),   // y6 → Shared
        (800.0, 200.0),         // noise → Unassigned
    ];

    // Compute precursor m/z for DGFLLDGFPR, charge 2
    let precursor_mass: f64 = "DGFLLDGFPR".chars().map(aa_mass).sum::<f64>() + WATER;
    let precursor_mz = (precursor_mass + 2.0 * PROTON) / 2.0;

    // Scan 2: spectrum for target PSM "PEPTIDEK" (will be Target, not traced)
    let scan2_peaks = vec![
        (200.0, 400.0),
        (350.0, 500.0),
        (500.0, 600.0),
    ];

    // Scan 3: spectrum for L4 trap PSM "WWWWWWWK" (no target match)
    // Include some of its b-ions so provenance can be traced (trap-only ions)
    let l4_b = b_ions("WWWWWWWK");
    let l4_precursor_mass: f64 = "WWWWWWWK".chars().map(aa_mass).sum::<f64>() + WATER;
    let l4_precursor_mz = (l4_precursor_mass + 2.0 * PROTON) / 2.0;
    let scan3_peaks = vec![
        (l4_b[0], 500.0),  // b1 of WWWWWWWK
        (l4_b[1], 600.0),  // b2
        (l4_b[2], 700.0),  // b3
        (l4_b[3], 800.0),  // b4
        (400.0, 200.0),    // noise
    ];

    let mut mgf = String::new();

    // Scan 1
    mgf.push_str("BEGIN IONS\n");
    mgf.push_str("TITLE=scan=1\n");
    mgf.push_str(&format!("PEPMASS={:.4} 10000.0\n", precursor_mz));
    mgf.push_str("CHARGE=2+\n");
    mgf.push_str("RTINSECONDS=300.0\n");
    for (mz, intensity) in &scan1_peaks {
        mgf.push_str(&format!("{:.4} {:.0}\n", mz, intensity));
    }
    mgf.push_str("END IONS\n");

    // Scan 2
    mgf.push_str("BEGIN IONS\n");
    mgf.push_str("TITLE=scan=2\n");
    mgf.push_str("PEPMASS=500.0 8000.0\n");
    mgf.push_str("CHARGE=2+\n");
    mgf.push_str("RTINSECONDS=350.0\n");
    for (mz, intensity) in &scan2_peaks {
        mgf.push_str(&format!("{:.4} {:.0}\n", mz, intensity));
    }
    mgf.push_str("END IONS\n");

    // Scan 3
    mgf.push_str("BEGIN IONS\n");
    mgf.push_str("TITLE=scan=3\n");
    mgf.push_str(&format!("PEPMASS={:.4} 6000.0\n", l4_precursor_mz));
    mgf.push_str("CHARGE=2+\n");
    mgf.push_str("RTINSECONDS=400.0\n");
    for (mz, intensity) in &scan3_peaks {
        mgf.push_str(&format!("{:.4} {:.0}\n", mz, intensity));
    }
    mgf.push_str("END IONS\n");

    mgf
}

// ===========================================================================
// Test 1: Unit-level provenance tracing with known ions
// ===========================================================================

#[test]
fn test_provenance_trace_known_ions() {
    let trap_seq = "DGFLLDGFPR";
    let target_seq = "NGFLLDGFPR";

    let trap_b = b_ions(trap_seq);
    let target_b = b_ions(target_seq);
    let shared_y = y_ions(trap_seq);

    // Build observed spectrum: 2 trap b-ions, 1 target b-ion, 3 shared y-ions, 2 noise
    let observed_mz = vec![
        trap_b[1],     // trap b2 → TrapOnly
        target_b[2],   // target b3 → TargetOnly
        shared_y[1],   // y2 → Shared
        trap_b[3],     // trap b4 → TrapOnly
        shared_y[3],   // y4 → Shared
        500.0,         // noise → Unassigned
        shared_y[5],   // y6 → Shared
        800.0,         // noise → Unassigned
    ];
    let observed_int = vec![600.0, 500.0, 700.0, 800.0, 750.0, 300.0, 650.0, 200.0];

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    let prov = trace_provenance(
        &observed_mz,
        &observed_int,
        trap_seq,
        target_seq,
        &[],
        &tolerance,
        1, // max_fragment_charge=1
    );

    // Verify counts
    assert_eq!(prov.trap_matched_count, 2, "expected 2 TrapOnly peaks");
    assert_eq!(prov.target_matched_count, 1, "expected 1 TargetOnly peak");
    assert_eq!(prov.shared_count, 3, "expected 3 Shared peaks");
    assert_eq!(prov.unassigned_count, 2, "expected 2 Unassigned peaks");

    // Verify shared_ratio = 3 / (2 + 1 + 3) = 0.5
    let expected_ratio = 3.0 / 6.0;
    assert!(
        (prov.shared_ratio - expected_ratio).abs() < 1e-6,
        "shared_ratio: expected {expected_ratio}, got {}",
        prov.shared_ratio
    );

    // Verify per-peak origins
    assert_eq!(prov.annotated_peaks[0].origin, IonOrigin::TrapOnly);
    assert_eq!(prov.annotated_peaks[1].origin, IonOrigin::TargetOnly);
    assert_eq!(prov.annotated_peaks[2].origin, IonOrigin::Shared);
    assert_eq!(prov.annotated_peaks[5].origin, IonOrigin::Unassigned);

    // Verify ion labels are populated
    assert!(prov.annotated_peaks[0].trap_ion_label.is_some());
    assert!(prov.annotated_peaks[0].target_ion_label.is_none());
    assert!(prov.annotated_peaks[1].trap_ion_label.is_none());
    assert!(prov.annotated_peaks[1].target_ion_label.is_some());
    assert!(prov.annotated_peaks[2].trap_ion_label.is_some());
    assert!(prov.annotated_peaks[2].target_ion_label.is_some());
}

// ===========================================================================
// Test 2: L4 provenance (no target match — all peaks are TrapOnly/Unassigned)
// ===========================================================================

#[test]
fn test_provenance_l4_no_target() {
    let trap_seq = "WWWWWWWK";
    let target_seq = ""; // L4 — no target

    let trap_b = b_ions(trap_seq);

    let observed_mz = vec![trap_b[0], trap_b[1], trap_b[2], trap_b[3], 400.0];
    let observed_int = vec![500.0, 600.0, 700.0, 800.0, 200.0];

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    let prov = trace_provenance(
        &observed_mz,
        &observed_int,
        trap_seq,
        target_seq,
        &[],
        &tolerance,
        1,
    );

    assert_eq!(prov.trap_matched_count, 4, "4 b-ions should be TrapOnly");
    assert_eq!(prov.target_matched_count, 0, "no target ions");
    assert_eq!(prov.shared_count, 0, "no shared ions");
    assert_eq!(prov.unassigned_count, 1, "1 noise peak");
    assert!((prov.shared_ratio - 0.0).abs() < 1e-6);
}

// ===========================================================================
// Test 3: Full pipeline — classify + provenance batch
// ===========================================================================

#[test]
fn test_e2e_classify_and_provenance_batch() {
    // --- 1. Create temp files ---
    let tmp_dir = tempfile::tempdir().unwrap();

    // FASTA: target proteins containing NGFLLDGFPR and PEPTIDEK
    let fasta_path = tmp_dir.path().join("target.fasta");
    std::fs::write(
        &fasta_path,
        ">sp|P001|PRTN_HUMAN Test human protein\n\
         NGFLLDGFPRLASTK\n\
         >sp|P002|PEPT_HUMAN Another human protein\n\
         PEPTIDEKLASTK\n",
    )
    .unwrap();

    // Config YAML
    let config_yaml = make_config_yaml();
    let config_path = tmp_dir.path().join("config.yaml");
    std::fs::write(&config_path, &config_yaml).unwrap();

    // MGF file with synthetic spectra
    let mgf_content = build_test_mgf();
    let mgf_path = tmp_dir.path().join("test_spectra.mgf");
    std::fs::write(&mgf_path, &mgf_content).unwrap();

    // --- 2. Build analyser ---
    let config = EntrapmentConfig::from_yaml(&config_path).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config.clone(), &fasta_path).unwrap();

    // --- 3. Create PSMs ---
    let psms = vec![
        // Target PSM (should be classified as Target, not provenance-traced)
        make_psm(
            "PEPTIDEK",
            "sp|P002|PEPT_HUMAN",
            Some(2),
            Some("test_spectra.mgf"),
        ),
        // Trap L2 PSM: DGFLLDGFPR similar to target NGFLLDGFPR (D→N, delta ≈ 0.98 Da)
        make_psm(
            "DGFLLDGFPR",
            "sp|Q001|TRAP_YEAST",
            Some(1),
            Some("test_spectra.mgf"),
        ),
        // Trap L4 PSM: completely different peptide
        make_psm(
            "WWWWWWWK",
            "sp|Q002|TRAP_ECOLI",
            Some(3),
            Some("test_spectra.mgf"),
        ),
    ];

    // --- 4. Classify all PSMs ---
    let mut classified = analyzer.classify_all(&psms).unwrap();
    assert_eq!(classified.len(), 3);

    // Verify classification groups and levels
    // PSM 0: Target
    assert_eq!(classified[0].group, PsmGroup::Target);

    // PSM 1: Trap, L2 (near-isobaric substitution D→N)
    assert_eq!(classified[1].group, PsmGroup::Trap);
    assert!(
        classified[1].level == DiscriminabilityLevel::L2
            || classified[1].level == DiscriminabilityLevel::L3,
        "DGFLLDGFPR should be L2 or L3 (D→N sub), got {:?}",
        classified[1].level
    );

    // PSM 2: Trap, L4 (no match)
    assert_eq!(classified[2].group, PsmGroup::Trap);
    assert_eq!(classified[2].level, DiscriminabilityLevel::L4);

    // --- 5. Run provenance tracing ---
    let traced = trace_provenance_batch(
        &mut classified,
        tmp_dir.path(), // mzml_dir — MGF is found via bare-name fallback
        &config,
    )
    .unwrap();

    // PSM 0 (Target) should NOT be traced
    assert!(
        classified[0].provenance.is_none(),
        "target PSMs should not have provenance"
    );

    // PSM 1 (Trap L2) should be traced
    assert!(
        classified[1].provenance.is_some(),
        "L2 trap PSM should have provenance traced"
    );
    let prov1 = classified[1].provenance.as_ref().unwrap();
    assert_eq!(prov1.trap_sequence, "DGFLLDGFPR");
    assert_eq!(prov1.target_sequence, "NGFLLDGFPR");
    assert!(prov1.trap_matched_count > 0, "should have trap-only ions");
    assert!(prov1.shared_count > 0, "should have shared ions (y-ions)");
    // shared_ratio > 0.3 → is_chimeric should be true
    assert!(
        prov1.is_chimeric,
        "shared_ratio={:.2} should exceed chimera threshold 0.3",
        prov1.shared_ratio
    );

    // PSM 2 (Trap L4) should be traced (L4 is in levels_to_trace)
    assert!(
        classified[2].provenance.is_some(),
        "L4 trap PSM should have provenance traced"
    );
    let prov2 = classified[2].provenance.as_ref().unwrap();
    assert_eq!(prov2.trap_sequence, "WWWWWWWK");
    assert_eq!(prov2.target_sequence, "");
    assert!(prov2.trap_matched_count > 0, "L4 should have trap-only b-ions");
    assert_eq!(prov2.target_matched_count, 0);
    assert_eq!(prov2.shared_count, 0);
    assert!(!prov2.is_chimeric);

    // Verify total traced count
    assert_eq!(traced, 2, "should have traced 2 PSMs (L2 + L4)");
}

// ===========================================================================
// Test 4: Output files — classified.tsv, report, metadata
// ===========================================================================

#[test]
fn test_e2e_output_files() {
    let tmp_dir = tempfile::tempdir().unwrap();

    // FASTA
    let fasta_path = tmp_dir.path().join("target.fasta");
    std::fs::write(
        &fasta_path,
        ">sp|P001|PRTN_HUMAN Test\nNGFLLDGFPRLASTK\n",
    )
    .unwrap();

    // Config
    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config.clone(), &fasta_path).unwrap();

    // MGF
    let mgf_content = build_test_mgf();
    std::fs::write(tmp_dir.path().join("test_spectra.mgf"), &mgf_content).unwrap();

    // PSMs
    let psms = vec![
        make_psm("DGFLLDGFPR", "sp|Q001|TRAP_YEAST", Some(1), Some("test_spectra.mgf")),
        make_psm("WWWWWWWK", "sp|Q002|TRAP_ECOLI", Some(3), Some("test_spectra.mgf")),
    ];

    let mut classified = analyzer.classify_all(&psms).unwrap();
    trace_provenance_batch(&mut classified, tmp_dir.path(), &config).unwrap();

    // Write classified TSV
    let out_dir = tmp_dir.path().join("output");
    std::fs::create_dir_all(&out_dir).unwrap();

    let classified_path = out_dir.join("classified.tsv");
    write_classified_tsv(&classified, &classified_path).unwrap();
    assert!(classified_path.exists(), "classified.tsv should exist");

    // Read and verify TSV content
    let tsv_content = std::fs::read_to_string(&classified_path).unwrap();
    let lines: Vec<&str> = tsv_content.lines().collect();
    assert_eq!(lines.len(), 3, "header + 2 data rows");

    // Header should include provenance columns
    let header = lines[0];
    assert!(header.contains("trap_matched"));
    assert!(header.contains("target_matched"));
    assert!(header.contains("shared_ions"));
    assert!(header.contains("shared_ratio"));
    assert!(header.contains("is_chimeric"));

    // Data rows should have provenance values (non-empty)
    for line in &lines[1..] {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 23, "expected 23 columns, got {}", fields.len());
        // trap_matched (index 18) should be non-empty
        assert!(!fields[18].is_empty(), "trap_matched should be populated");
    }

    // Write HTML report
    let summary = analyzer.summary(&classified);
    let report_path = out_dir.join("report.html");
    report::render_report(&summary, &classified, &report_path).unwrap();
    assert!(report_path.exists(), "HTML report should exist");

    let html = std::fs::read_to_string(&report_path).unwrap();
    assert!(html.contains("<!DOCTYPE html>") || html.contains("<html"));
    assert!(html.len() > 100, "HTML report should have substantial content");

    // Write run metadata
    let metadata = RunMetadata {
        tool_version: "test".to_string(),
        run_timestamp: "2025-01-01T00:00:00Z".to_string(),
        input_file: "test.tsv".to_string(),
        input_sha256: "0".repeat(64),
        fasta_file: fasta_path.display().to_string(),
        fasta_sha256: "0".repeat(64),
        config_snapshot: serde_json::to_value(&config).unwrap(),
        total_psms: classified.len(),
        trap_psms: 2,
        level_counts: summary.level_counts.clone(),
    };
    let metadata_path = out_dir.join("run_metadata.json");
    write_run_metadata(&metadata, &metadata_path).unwrap();
    assert!(metadata_path.exists(), "run_metadata.json should exist");

    let meta_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&metadata_path).unwrap()).unwrap();
    assert_eq!(meta_json["total_psms"], 2);
}

// ===========================================================================
// Test 5: Provenance with modifications
// ===========================================================================

#[test]
fn test_provenance_with_modifications() {
    // Test that modifications shift theoretical ion masses correctly.
    // Trap: DGFLLDGFPR with Oxidation(M) at position... wait, no M in this peptide.
    // Use a simpler peptide: PEPTMDER with oxidation at position 4 (M)
    let trap_seq = "PEPTMDER";
    let target_seq = ""; // L4

    // Oxidation of M: +15.994915
    let mods = vec![(4_usize, 15.994915_f64)];

    // Compute expected b5 with mod: P+E+P+T+M(ox) = 97.05276+129.04259+97.05276+101.04768+(131.04049+15.994915)
    // = 97.05276 + 129.04259 + 97.05276 + 101.04768 + 147.035405 + 1.007276 = 572.235871
    // b4 without mod: P+E+P+T = 97.05276+129.04259+97.05276+101.04768 + 1.007276 = 425.200466
    // b5 with mod: b4_sum + (131.04049 + 15.994915) + 1.007276 = 572.235871

    let b4_mz = 97.05276 + 129.04259 + 97.05276 + 101.04768 + PROTON;
    let b5_mz = b4_mz - PROTON + (131.04049 + 15.994915) + PROTON;

    let observed_mz = vec![b4_mz, b5_mz, 400.0];
    let observed_int = vec![500.0, 600.0, 200.0];

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    let prov = trace_provenance(
        &observed_mz,
        &observed_int,
        trap_seq,
        target_seq,
        &mods,
        &tolerance,
        1,
    );

    // b4 should match (unmodified prefix), b5 should match (includes mod)
    assert_eq!(prov.trap_matched_count, 2, "b4 and b5 should match with mod");
    assert_eq!(prov.unassigned_count, 1, "noise peak");
}

// ===========================================================================
// Test 6: Provenance batch skips PSMs without scan_number (B2 fix)
// ===========================================================================

#[test]
fn test_provenance_batch_no_scan_number_warns() {
    let tmp_dir = tempfile::tempdir().unwrap();

    let fasta_path = tmp_dir.path().join("target.fasta");
    std::fs::write(
        &fasta_path,
        ">sp|P001|PRTN_HUMAN Test\nNGFLLDGFPRLASTK\n",
    )
    .unwrap();

    let config = EntrapmentConfig::from_yaml_str(&make_config_yaml()).unwrap();
    let analyzer = EntrapmentAnalyzer::new(config.clone(), &fasta_path).unwrap();

    // PSM without scan_number — simulates DIA-NN results
    let psms = vec![make_psm("DGFLLDGFPR", "sp|Q001|TRAP_YEAST", None, None)];

    let mut classified = analyzer.classify_all(&psms).unwrap();
    assert!(classified[0].level != DiscriminabilityLevel::L4);

    // Provenance batch should not crash, should return 0 traced
    let traced = trace_provenance_batch(&mut classified, tmp_dir.path(), &config).unwrap();
    assert_eq!(traced, 0, "no PSMs should be traced without scan_number");
    assert!(
        classified[0].provenance.is_none(),
        "provenance should remain None"
    );
}

// ===========================================================================
// Test 7: Chimera threshold logic (B3 fix)
// ===========================================================================

#[test]
fn test_chimera_threshold_flag() {
    let trap_seq = "DGFLLDGFPR";
    let target_seq = "NGFLLDGFPR";

    // Build a spectrum with many shared y-ions, no trap-only/target-only
    let shared_y = y_ions(trap_seq);
    let observed_mz = vec![shared_y[0], shared_y[1], shared_y[2]];
    let observed_int = vec![500.0; 3];

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    let prov = trace_provenance(
        &observed_mz,
        &observed_int,
        trap_seq,
        target_seq,
        &[],
        &tolerance,
        1,
    );

    // All peaks are shared → shared_ratio = 1.0
    assert_eq!(prov.shared_count, 3);
    assert!((prov.shared_ratio - 1.0).abs() < 1e-6);

    // is_chimeric is set by caller (trace_provenance returns false)
    assert!(!prov.is_chimeric, "trace_provenance always returns false");

    // Caller applies threshold: shared_ratio > 0.3 → chimeric
    let is_chimeric = prov.shared_ratio > 0.3;
    assert!(is_chimeric, "shared_ratio=1.0 > 0.3 → chimeric");
}

// ===========================================================================
// Test 8: Mirror plot generation (HTML output)
// ===========================================================================

#[test]
fn test_mirror_plot_html_generation() {
    use protein_copilot_entrapment_analysis::mirror_plot::render_mirror_plot;

    let trap_seq = "DGFLLDGFPR";
    let target_seq = "NGFLLDGFPR";

    let trap_b = b_ions(trap_seq);
    let shared_y = y_ions(trap_seq);

    let observed_mz = vec![trap_b[1], shared_y[1], 500.0];
    let observed_int = vec![600.0, 700.0, 300.0];

    let tolerance = MassTolerance {
        value: 20.0,
        unit: ToleranceUnit::Ppm,
    };

    let prov = trace_provenance(
        &observed_mz,
        &observed_int,
        trap_seq,
        target_seq,
        &[],
        &tolerance,
        1,
    );

    let tmp_dir = tempfile::tempdir().unwrap();
    let html_path = tmp_dir.path().join("mirror_plot.html");

    render_mirror_plot(&prov, &html_path).unwrap();

    assert!(html_path.exists());
    let html = std::fs::read_to_string(&html_path).unwrap();
    assert!(html.contains("TrapOnly") || html.contains("trap"));
    assert!(html.len() > 500, "mirror plot HTML should be substantial");
}
