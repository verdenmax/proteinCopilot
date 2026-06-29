//! Integration test: full search pipeline (read → recommend → search → FDR → summary → annotate).
//!
//! Uses the test_100 fixtures from the search-engine crate.

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::progress::noop_progress;
use protein_copilot_core::search_params::MassTolerance;
use protein_copilot_param_recommend::ParamRecommender;
use protein_copilot_report::ReportGenerator;
use protein_copilot_search_engine::annotate::annotate_spectrum;
use protein_copilot_search_engine::SimpleSearchEngine;
use protein_copilot_spectrum_io::{create_reader, detect_format};
use test_helpers::search_engine_fixtures;

fn mgf_path() -> std::path::PathBuf {
    search_engine_fixtures()
        .parent()
        .unwrap()
        .join("fixtures_e2e")
        .join("test_100.mgf")
}

fn fasta_path() -> std::path::PathBuf {
    search_engine_fixtures()
        .parent()
        .unwrap()
        .join("fixtures_e2e")
        .join("test_100.fasta")
}

/// Full pipeline: read → recommend → search → FDR summary → annotate best PSM
#[tokio::test]
async fn full_pipeline_read_search_annotate() {
    let mgf = mgf_path();
    let fasta = fasta_path();

    assert!(
        mgf.exists() && fasta.exists(),
        "required fixtures missing: {mgf:?} / {fasta:?}"
    );

    // Step 1: Read spectra
    let file_info = detect_format(&mgf).unwrap();
    let reader = create_reader(&file_info);
    let summary = reader.read_summary(&mgf).unwrap();
    assert!(summary.ms2_count > 0, "should have MS2 spectra");

    // Step 2: Recommend parameters
    let rec = ParamRecommender.recommend(&summary, None).unwrap();
    let mut params = rec.decision;
    params.database_path = fasta.to_string_lossy().to_string();
    assert!(params.validate().is_ok());

    // Step 3: Search
    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(
            &params,
            std::slice::from_ref(&mgf),
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
        .await
        .unwrap();
    assert!(!result.psms.is_empty(), "should find some PSMs");

    // Step 4: FDR-filtered summary
    let fdr_summary = ReportGenerator::generate_summary(&result);
    assert!(fdr_summary.identification_rate >= 0.0);
    assert!(fdr_summary.identification_rate <= 1.0);

    // Step 5: Annotate the best-scoring PSM
    let best_psm = result
        .psms
        .iter()
        .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
        .unwrap();

    let spectrum = reader.read_spectrum(&mgf, best_psm.spectrum_scan).unwrap();
    let tol = MassTolerance {
        value: 20.0,
        unit: protein_copilot_core::search_params::ToleranceUnit::Ppm,
    };

    let annotation = annotate_spectrum(
        &spectrum,
        &best_psm.peptide_sequence,
        best_psm.charge,
        &tol,
        &best_psm.modifications,
        best_psm.protein_accessions.clone(),
        false,
        false,
    )
    .unwrap();

    assert!(annotation.matched_ions > 0, "best PSM should match ions");
    assert!(!annotation.b_ions.is_empty());
    assert!(!annotation.y_ions.is_empty());

    // Step 6: Export to temp dir
    let dir = tempfile::tempdir().unwrap();
    ReportGenerator::export_tsv(&result, dir.path()).unwrap();
    assert!(dir.path().join("psm.tsv").exists());
    assert!(dir.path().join("peptide.tsv").exists());
    assert!(dir.path().join("protein.tsv").exists());
}
