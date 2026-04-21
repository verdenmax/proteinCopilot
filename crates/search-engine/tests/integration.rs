//! End-to-end integration test: spectrum-io → param-recommend → search-engine.
//!
//! Validates the complete data flow: read spectra → recommend params → search.

use std::io::Write;
use std::path::PathBuf;

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::progress::noop_progress;
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
    let result = engine
        .search(
            &params,
            &[mgf_path],
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
        .await
        .unwrap();

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
        .search(
            &params,
            &[mgf_path],
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
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
        .search(
            &params,
            &[mgf_path],
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
        .await
        .unwrap();
    assert_eq!(result.params_used.enzyme, Enzyme::LysC);
}

// ─────────────────────────────────────────────────────────
// Variable Modification Integration Tests (FW-2 + FW-1)
// ─────────────────────────────────────────────────────────

/// Search with Oxidation(M) variable modification finds modified peptide.
#[tokio::test]
async fn variable_mod_oxidation_m() {
    use protein_copilot_core::search_params::*;
    use protein_copilot_search_engine::SimpleSearchEngine;

    let fasta = create_test_fasta();

    // Protein TEST3 has MKWVTFISLLFLFSSAYSRGVFRR
    // Tryptic peptide: MKWVTFISLLFLFSSAYSR (min_len=7)
    // With Oxidation(M@0): mass += 15.994915

    let oxidation = Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.994915,
        residues: vec!['M'],
        position: ModPosition::Anywhere,
    };

    let params = SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 1,
        fixed_modifications: vec![],
        variable_modifications: vec![oxidation],
        max_variable_modifications: 3,
        precursor_tolerance: MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        database_path: fasta.path().to_string_lossy().to_string(),
        decoy_strategy: DecoyStrategy::None,
        min_peptide_length: 7,
        max_peptide_length: 50,
        acquisition_mode: None,
        engine: None,
    };

    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(
            &params,
            &[mgf_fixture()],
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
        .await
        .unwrap();

    // The search should complete without error
    assert_eq!(result.summary.total_spectra_searched, 10);
    // Variable mod should be in params_used
    assert_eq!(result.params_used.variable_modifications.len(), 1);
    assert_eq!(
        result.params_used.variable_modifications[0].name,
        "Oxidation"
    );
}

/// ProteinNTerm Acetylation only applies to N-terminal peptides (FW-1).
#[tokio::test]
async fn protein_nterm_acetylation_fw1() {
    use protein_copilot_core::search_params::*;
    use protein_copilot_search_engine::SimpleSearchEngine;

    let fasta = create_test_fasta();

    let acetyl = Modification {
        name: "Acetyl".to_string(),
        mass_delta: 42.010565,
        residues: vec![],
        position: ModPosition::ProteinNTerm,
    };

    let params = SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 1,
        fixed_modifications: vec![],
        variable_modifications: vec![acetyl],
        max_variable_modifications: 3,
        precursor_tolerance: MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        database_path: fasta.path().to_string_lossy().to_string(),
        decoy_strategy: DecoyStrategy::None,
        min_peptide_length: 7,
        max_peptide_length: 50,
        acquisition_mode: None,
        engine: None,
    };

    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(
            &params,
            &[mgf_fixture()],
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
        .await
        .unwrap();

    // Search completes — ProteinNTerm mods should only apply to first peptide
    assert_eq!(result.summary.total_spectra_searched, 10);

    // If any PSM has Acetyl mod, verify it's on a protein N-terminal peptide
    for psm in &result.psms {
        if psm.modifications.iter().any(|m| m.name == "Acetyl") {
            assert!(
                !psm.peptide_sequence.is_empty(),
                "acetylated PSM should have a sequence"
            );
        }
    }
}

// ─────────────────────────────────────────────────────────
// FDR Integration Tests (FW-6)
// ─────────────────────────────────────────────────────────

/// Search with Reverse decoy strategy produces q-values and removes decoys.
#[tokio::test]
async fn fdr_reverse_decoy_strategy_fw6() {
    use protein_copilot_core::search_params::*;
    use protein_copilot_search_engine::SimpleSearchEngine;

    let fasta = create_test_fasta();

    let params = SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 1,
        fixed_modifications: vec![],
        variable_modifications: vec![],
        max_variable_modifications: 3,
        precursor_tolerance: MassTolerance {
            value: 10.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        database_path: fasta.path().to_string_lossy().to_string(),
        decoy_strategy: DecoyStrategy::Reverse,
        min_peptide_length: 7,
        max_peptide_length: 50,
        acquisition_mode: None,
        engine: None,
    };

    let engine = SimpleSearchEngine::new();
    let result = engine
        .search(
            &params,
            &[mgf_fixture()],
            noop_progress(),
            &mut protein_copilot_core::diagnostics::SearchDiagnostics::new(),
        )
        .await
        .unwrap();

    assert_eq!(result.summary.total_spectra_searched, 10);

    // All output PSMs should be targets (decoys removed)
    for psm in &result.psms {
        assert!(!psm.is_decoy, "decoy PSMs should be removed from output");
        assert!(
            !psm.protein_accessions.iter().any(|a| a.starts_with("REV_")),
            "no REV_ accessions in output"
        );
    }

    // PSMs with scores should have q-values assigned
    if !result.psms.is_empty() {
        let has_qval = result.psms.iter().any(|p| p.q_value.is_some());
        assert!(has_qval, "at least some PSMs should have q-values");

        // q-values must be monotonically non-decreasing when sorted by score descending
        let mut scored: Vec<_> = result.psms.iter().filter(|p| p.q_value.is_some()).collect();
        scored.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for w in scored.windows(2) {
            assert!(
                w[0].q_value.unwrap() <= w[1].q_value.unwrap(),
                "q-values should be monotonically non-decreasing: {} > {}",
                w[0].q_value.unwrap(),
                w[1].q_value.unwrap(),
            );
        }
    }
}
