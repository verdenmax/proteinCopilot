//! Cross-module integration test: spectrum-io → param-recommend.
//!
//! Reads a real fixture file, generates a summary, then feeds it
//! to the parameter recommender — verifying the full pipeline.

use std::path::PathBuf;

use protein_copilot_param_recommend::{ParamRecommender, UserHints};
use protein_copilot_spectrum_io::{create_reader, detect_format};

fn fixtures_dir() -> PathBuf {
    // spectrum-io fixtures are in a sibling crate
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("spectrum-io")
        .join("tests")
        .join("fixtures")
}

#[test]
fn mgf_to_recommendation() {
    let path = fixtures_dir().join("small.mgf");
    let info = detect_format(&path).unwrap();
    let reader = create_reader(&info);
    let summary = reader.read_summary(&path).unwrap();

    let recommender = ParamRecommender;
    let result = recommender.recommend(&summary, None).unwrap();

    assert!(result.confidence > 0.0);
    assert!(!result.explanation.is_empty());
    assert!(!result.evidence.is_empty());
    assert!(!result.alternatives.is_empty());
    assert!(!result.input_summary.is_empty());

    // Summary has 10 spectra → should produce a valid recommendation
    assert!(result.input_summary.contains("10 spectra"));
}

#[test]
fn mzml_to_recommendation() {
    let path = fixtures_dir().join("small.mzml");
    let info = detect_format(&path).unwrap();
    let reader = create_reader(&info);
    let summary = reader.read_summary(&path).unwrap();

    let recommender = ParamRecommender;
    let result = recommender.recommend(&summary, None).unwrap();

    assert!(result.confidence > 0.0);
    assert!(!result.explanation.is_empty());
}

#[test]
fn mgf_with_phospho_hint() {
    let path = fixtures_dir().join("small.mgf");
    let info = detect_format(&path).unwrap();
    let summary = create_reader(&info).read_summary(&path).unwrap();

    let hints = UserHints {
        experiment_type: Some("phosphorylation".to_string()),
        ..Default::default()
    };

    let recommender = ParamRecommender;
    let result = recommender.recommend(&summary, Some(&hints)).unwrap();

    // Should include Phospho modification
    assert!(result
        .decision
        .variable_modifications
        .iter()
        .any(|m| m.name == "Phospho"));

    // With hint, confidence should be ≥ 0.90
    assert!(
        result.confidence >= 0.90,
        "confidence with hint: {}",
        result.confidence
    );
}

#[test]
fn recommended_params_are_valid_with_database() {
    let path = fixtures_dir().join("small.mgf");
    let info = detect_format(&path).unwrap();
    let summary = create_reader(&info).read_summary(&path).unwrap();

    let recommender = ParamRecommender;
    let result = recommender.recommend(&summary, None).unwrap();

    // Set a real database path and validate
    let mut params = result.decision;
    params.database_path = "/data/human.fasta".to_string();
    assert!(params.validate().is_ok());
}
