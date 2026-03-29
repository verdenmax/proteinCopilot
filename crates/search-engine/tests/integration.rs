//! End-to-end integration test: spectrum-io → param-recommend → search-engine.
//!
//! Validates the complete data flow: read spectra → recommend params → search.

use std::io::Write;
use std::path::PathBuf;

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::search_params::Enzyme;
use protein_copilot_param_recommend::{ParamRecommender, UserHints};
use protein_copilot_search_engine::SimpleSearchEngine;
use protein_copilot_spectrum_io::{create_reader, detect_format};

fn mgf_fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("spectrum-io")
        .join("tests")
        .join("fixtures")
        .join("small.mgf")
}

fn create_test_fasta() -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(
        f,
        ">sp|P001|TEST1 Test protein 1\n\
         PEPTIDEKANOTHERRLASTR\n\
         >sp|P002|TEST2 Test protein 2\n\
         AVCDEFGKHIKLMNPQRST\n\
         >sp|P003|TEST3 Test protein 3\n\
         MKWVTFISLLFLFSSAYSRGVFRR\n"
    )
    .unwrap();
    f
}

/// Full pipeline: read file → summary → recommend params → search
#[tokio::test]
async fn full_pipeline_mgf_to_search_result() {
    let mgf_path = mgf_fixture();
    let fasta = create_test_fasta();

    // Step 1: Read spectrum summary
    let file_info = detect_format(&mgf_path).unwrap();
    let reader = create_reader(&file_info);
    let summary = reader.read_summary(&mgf_path).unwrap();
    assert_eq!(summary.total_spectra, 10);

    // Step 2: Recommend parameters
    let recommender = ParamRecommender;
    let recommendation = recommender.recommend(&summary, None).unwrap();
    assert!(recommendation.confidence > 0.0);

    // Step 3: Set real database path on recommended params
    let mut params = recommendation.decision;
    params.database_path = fasta.path().to_string_lossy().to_string();
    assert!(params.validate().is_ok());

    // Step 4: Run search
    let engine = SimpleSearchEngine::new();
    let result = engine.search(&params, &[mgf_path]).await.unwrap();

    // Verify result structure
    assert_eq!(result.summary.total_spectra_searched, 10);
    assert!(result.summary.search_duration_sec >= 0.0);
    assert!(result.summary.identification_rate >= 0.0);
    assert!(result.summary.identification_rate <= 1.0);
    assert_eq!(result.run_id, result.metadata.run_id);
    assert_eq!(result.engine_info, result.metadata.engine_info);
    assert_eq!(result.params_used, result.metadata.params_used);
    assert_eq!(
        result.metadata.status,
        protein_copilot_core::run_metadata::RunStatus::Completed
    );
}

/// Pipeline with user hints: phospho experiment type
#[tokio::test]
async fn full_pipeline_with_phospho_hints() {
    let mgf_path = mgf_fixture();
    let fasta = create_test_fasta();

    let file_info = detect_format(&mgf_path).unwrap();
    let summary = create_reader(&file_info).read_summary(&mgf_path).unwrap();

    let hints = UserHints {
        experiment_type: Some("phosphorylation".to_string()),
        ..Default::default()
    };
    let recommendation = ParamRecommender.recommend(&summary, Some(&hints)).unwrap();

    // Phospho preset should include Phospho modification
    assert!(recommendation
        .decision
        .variable_modifications
        .iter()
        .any(|m| m.name == "Phospho"));

    let mut params = recommendation.decision;
    params.database_path = fasta.path().to_string_lossy().to_string();

    let result = SimpleSearchEngine::new()
        .search(&params, &[mgf_path])
        .await
        .unwrap();
    assert_eq!(result.summary.total_spectra_searched, 10);
    assert_eq!(result.params_used.enzyme, Enzyme::Trypsin);
}

/// Pipeline with custom enzyme hint
#[tokio::test]
async fn full_pipeline_with_enzyme_override() {
    let mgf_path = mgf_fixture();
    let fasta = create_test_fasta();

    let file_info = detect_format(&mgf_path).unwrap();
    let summary = create_reader(&file_info).read_summary(&mgf_path).unwrap();

    let hints = UserHints {
        enzyme: Some(Enzyme::LysC),
        ..Default::default()
    };
    let recommendation = ParamRecommender.recommend(&summary, Some(&hints)).unwrap();
    assert_eq!(recommendation.decision.enzyme, Enzyme::LysC);

    let mut params = recommendation.decision;
    params.database_path = fasta.path().to_string_lossy().to_string();

    let result = SimpleSearchEngine::new()
        .search(&params, &[mgf_path])
        .await
        .unwrap();
    assert_eq!(result.params_used.enzyme, Enzyme::LysC);
}
