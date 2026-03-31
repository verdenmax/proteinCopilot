//! End-to-end integration tests for the complete ProteinCopilot pipeline.
//!
//! Tests the full workflow: spectrum-io → param-recommend → search-engine → report
//! using realistic fixture data (100 spectra + 100 proteins).

use std::path::PathBuf;

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::progress::noop_progress;
use protein_copilot_core::search_params::*;
use protein_copilot_param_recommend::{ParamRecommender, UserHints};
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::SimpleSearchEngine;
use protein_copilot_spectrum_io::{create_reader, detect_format};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures_e2e")
}

fn mgf_path() -> PathBuf {
    fixtures_dir().join("test_100.mgf")
}

fn fasta_path() -> PathBuf {
    fixtures_dir().join("test_100.fasta")
}

// ─────────────────────────────────────────────────────────
// Scenario A: Full automated pipeline
// ─────────────────────────────────────────────────────────

#[tokio::test]
async fn scenario_a_full_pipeline() {
    // Step 1: Read spectra
    let file_info = detect_format(&mgf_path()).unwrap();
    let reader = create_reader(&file_info);
    let summary = reader.read_summary(&mgf_path()).unwrap();

    assert_eq!(summary.total_spectra, 100);
    assert_eq!(summary.ms2_count, 100);
    assert!(summary.mz_range[0] > 0.0);
    assert!(summary.mz_range[1] > summary.mz_range[0]);

    // Step 2: Recommend params
    let recommendation = ParamRecommender.recommend(&summary, None).unwrap();
    assert!(recommendation.confidence > 0.0);
    assert!(!recommendation.explanation.is_empty());

    // Step 3: Set database path
    let mut params = recommendation.decision;
    params.database_path = fasta_path().to_string_lossy().to_string();
    assert!(params.validate().is_ok());

    // Step 4: Run search
    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(&params, &[mgf_path()], noop_progress())
        .await
        .unwrap();

    // Step 5: Verify result structure
    assert_eq!(result.summary.total_spectra_searched, 100);
    assert!(result.summary.search_duration_sec >= 0.0);
    assert!(result.summary.identification_rate >= 0.0);
    assert!(result.summary.identification_rate <= 1.0);
    assert_eq!(result.run_id, result.metadata.run_id);
    assert_eq!(result.engine_info, result.metadata.engine_info);
    assert_eq!(result.params_used, result.metadata.params_used);

    // Should identify some PSMs (fixture data has matching b/y ions)
    assert!(
        !result.psms.is_empty(),
        "Expected some PSMs from matched fixture data"
    );

    // Step 6: Generate FDR-filtered summary
    let filtered_summary = ReportGenerator::generate_summary(&result);
    assert_eq!(
        filtered_summary.total_spectra_searched,
        result.summary.total_spectra_searched
    );

    // Step 7: Export results
    let output_dir = tempfile::tempdir().unwrap();
    ReportGenerator::export_tsv(&result, output_dir.path()).unwrap();
    ReportGenerator::export_json(&result, &output_dir.path().join("result.json")).unwrap();
    ReportGenerator::export_metadata(
        &result.metadata,
        &output_dir.path().join("run_metadata.json"),
    )
    .unwrap();

    // Verify files exist
    assert!(output_dir.path().join("psm.tsv").exists());
    assert!(output_dir.path().join("peptide.tsv").exists());
    assert!(output_dir.path().join("protein.tsv").exists());
    assert!(output_dir.path().join("result.json").exists());
    assert!(output_dir.path().join("run_metadata.json").exists());

    // Verify JSON roundtrip
    let json_content = std::fs::read_to_string(output_dir.path().join("result.json")).unwrap();
    let back: protein_copilot_core::search_result::SearchResult =
        serde_json::from_str(&json_content).unwrap();
    assert_eq!(result.run_id, back.run_id);
    assert_eq!(result.psms.len(), back.psms.len());
}

// ─────────────────────────────────────────────────────────
// Scenario B: Direct parameters (skip recommendation)
// ─────────────────────────────────────────────────────────

#[tokio::test]
async fn scenario_b_direct_params() {
    let params = SearchParams {
        database_path: fasta_path().to_string_lossy().to_string(),
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 1,
        fixed_modifications: vec![Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        }],
        variable_modifications: vec![],
        precursor_tolerance: MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        decoy_strategy: DecoyStrategy::Reverse,
    };

    assert!(params.validate().is_ok());

    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(&params, &[mgf_path()], noop_progress())
        .await
        .unwrap();

    assert_eq!(result.summary.total_spectra_searched, 100);
    assert_eq!(result.params_used.enzyme, Enzyme::Trypsin);
    assert!(!result.psms.is_empty());

    // Verify PSM fields are valid
    for psm in &result.psms {
        assert!(psm.spectrum_scan >= 1);
        assert!(!psm.peptide_sequence.is_empty());
        assert!(psm.charge != 0);
        assert!(psm.precursor_mz > 0.0);
        assert!(psm.score >= 0.0);
        assert!(!psm.protein_accessions.is_empty());
    }
}

// ─────────────────────────────────────────────────────────
// Scenario C: Data exploration only
// ─────────────────────────────────────────────────────────

#[test]
fn scenario_c_data_exploration() {
    // Just read and summarize — no search
    let file_info = detect_format(&mgf_path()).unwrap();
    let reader = create_reader(&file_info);

    // Summary
    let summary = reader.read_summary(&mgf_path()).unwrap();
    assert_eq!(summary.total_spectra, 100);
    assert!(summary.validate().is_ok());
    assert!(!summary.is_empty());

    // Read individual spectrum
    let spectrum = reader.read_spectrum(&mgf_path(), 1).unwrap();
    assert_eq!(spectrum.scan_number, 1);
    assert!(spectrum.validate().is_ok());
    assert!(!spectrum.mz_array.is_empty());
    assert_eq!(spectrum.precursors.len(), 1);

    // Charge distribution should have entries
    assert!(!summary.precursor_charge_distribution.is_empty());
}

// ─────────────────────────────────────────────────────────
// Scenario: Phospho experiment with hints
// ─────────────────────────────────────────────────────────

#[tokio::test]
async fn scenario_phospho_with_hints() {
    let file_info = detect_format(&mgf_path()).unwrap();
    let summary = create_reader(&file_info).read_summary(&mgf_path()).unwrap();

    let hints = UserHints {
        experiment_type: Some("phosphorylation".to_string()),
        instrument_type: Some("Orbitrap".to_string()),
        ..Default::default()
    };
    let recommendation = ParamRecommender.recommend(&summary, Some(&hints)).unwrap();

    // With both experiment_type + instrument_type hints, confidence should be >= 0.90
    assert!(
        recommendation.confidence >= 0.90,
        "Expected confidence >= 0.90 with hints, got {}",
        recommendation.confidence
    );
    assert!(recommendation
        .decision
        .variable_modifications
        .iter()
        .any(|m| m.name == "Phospho"));

    let mut params = recommendation.decision;
    params.database_path = fasta_path().to_string_lossy().to_string();

    let result = SimpleSearchEngine::new()
        .search(&params, &[mgf_path()], noop_progress())
        .await
        .unwrap();

    assert_eq!(result.summary.total_spectra_searched, 100);
}

// ─────────────────────────────────────────────────────────
// Error scenarios
// ─────────────────────────────────────────────────────────

#[test]
fn error_missing_file() {
    let err = detect_format(std::path::Path::new("/nonexistent/file.mgf")).unwrap_err();
    assert!(err.to_string().contains("not found"));
}

#[tokio::test]
async fn error_invalid_params() {
    let params = SearchParams {
        database_path: "".to_string(), // invalid
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 0,
        fixed_modifications: vec![],
        variable_modifications: vec![],
        precursor_tolerance: MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        decoy_strategy: DecoyStrategy::Reverse,
    };

    let engine = SimpleSearchEngine::new();
    let result = engine.search(&params, &[mgf_path()], noop_progress()).await;
    assert!(result.is_err());
}

// ─────────────────────────────────────────────────────────
// Scenario: Progress tracking e2e
// ─────────────────────────────────────────────────────────

#[tokio::test]
async fn scenario_progress_tracking() {
    use protein_copilot_core::progress::{ProgressCallback, SearchProgress};
    use std::sync::{Arc, Mutex};

    let stages: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let stages_clone = Arc::clone(&stages);
    let on_progress: ProgressCallback = Box::new(move |p: SearchProgress| {
        if let Some(ref stage) = p.stage {
            let mut s = stages_clone.lock().unwrap();
            if s.last().map(|l| l != stage).unwrap_or(true) {
                s.push(stage.clone());
            }
        }
        // Verify progress_pct is valid when present
        if let Some(pct) = p.progress_pct {
            assert!(
                (0.0..=1.0).contains(&pct),
                "progress_pct out of range: {pct}"
            );
        }
    });

    let file_info = detect_format(&mgf_path()).unwrap();
    let summary = create_reader(&file_info).read_summary(&mgf_path()).unwrap();
    let mut params = ParamRecommender.recommend(&summary, None).unwrap().decision;
    params.database_path = fasta_path().to_string_lossy().to_string();

    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(&params, &[mgf_path()], on_progress)
        .await
        .unwrap();

    assert_eq!(result.summary.total_spectra_searched, 100);

    let recorded = stages.lock().unwrap();
    assert!(
        recorded.len() >= 4,
        "Expected at least 4 stages, got: {recorded:?}"
    );
    assert!(
        recorded[0].contains("FASTA"),
        "First stage should mention FASTA"
    );
    assert!(
        recorded.iter().any(|s| s.contains("Matching")),
        "Should have a Matching stage"
    );
    assert!(
        recorded.last().unwrap().contains("Aggregating"),
        "Last stage should be Aggregating"
    );
}
