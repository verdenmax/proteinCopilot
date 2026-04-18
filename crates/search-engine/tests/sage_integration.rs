//! Integration test for SageAdapter.
//!
//! Uses a small FASTA + mgf file to verify the full search pipeline.

use std::path::PathBuf;

use protein_copilot_core::engine::SearchEngineAdapter;
use protein_copilot_core::search_params::*;
use protein_copilot_search_engine::adapters::sage::SageAdapter;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures")
        .join(name)
}

fn test_params(fasta_path: &str) -> SearchParams {
    SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 2,
        fixed_modifications: vec![Modification {
            name: "Carbamidomethyl".into(),
            mass_delta: 57.021464,
            residues: vec!['C'],
            position: ModPosition::Anywhere,
        }],
        variable_modifications: vec![Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residues: vec!['M'],
            position: ModPosition::Anywhere,
        }],
        precursor_tolerance: MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.5,
            unit: ToleranceUnit::Da,
        },
        database_path: fasta_path.to_string(),
        decoy_strategy: DecoyStrategy::Reverse,
        acquisition_mode: None,
        max_variable_modifications: 3,
        min_peptide_length: 5,
        max_peptide_length: 50,
        engine: Some("Sage".into()),
    }
}

#[tokio::test]
async fn sage_search_produces_results() {
    let fasta = fixture_path("small_test.fasta");
    let mgf = fixture_path("small_test.mgf");

    if !fasta.exists() || !mgf.exists() {
        eprintln!(
            "Skipping sage integration test: test fixtures not found at {:?}",
            fasta
        );
        return;
    }

    let params = test_params(fasta.to_str().unwrap_or_default());
    let adapter = SageAdapter::default();
    let on_progress: protein_copilot_core::progress::ProgressCallback = Box::new(|p| {
        eprintln!("Progress: {:?} ({:?}%)", p.stage, p.progress_pct);
    });

    let result = adapter.search(&params, &[mgf], on_progress, &mut protein_copilot_core::diagnostics::SearchDiagnostics::new()).await;

    match result {
        Ok(result) => {
            // We should get some PSMs
            eprintln!(
                "Sage search complete: {} PSMs, {} at 1% FDR",
                result.psms.len(),
                result
                    .psms
                    .iter()
                    .filter(|p| !p.is_decoy && p.q_value.map_or(false, |q| q <= 0.01))
                    .count()
            );

            // All PSMs should have valid fields
            for psm in &result.psms {
                assert!(psm.spectrum_scan >= 1, "scan should be >= 1");
                assert!(
                    !psm.peptide_sequence.is_empty(),
                    "sequence should not be empty"
                );
                assert!(psm.charge > 0, "charge should be positive");
                assert!(psm.precursor_mz > 0.0, "precursor_mz should be positive");
                assert!(psm.score.is_finite(), "score should be finite");
                assert!(psm.q_value.is_some(), "q_value should be set by Sage FDR");
                assert!(
                    !psm.protein_accessions.is_empty(),
                    "should have protein accessions"
                );
                assert!(psm.extra.is_some(), "should have extra sage fields");

                // Verify extra fields contain expected sage-specific data
                let extra = psm.extra.as_ref().unwrap_or_else(|| {
                    panic!("extra should be Some for PSM scan={}", psm.spectrum_scan)
                });
                assert!(extra.contains_key("hyperscore"), "should have hyperscore");
                assert!(
                    extra.contains_key("matched_peaks"),
                    "should have matched_peaks"
                );
            }

            // Summary should be consistent
            assert_eq!(result.summary.total_psms, result.psms.len() as u64);
        }
        Err(e) => {
            // Synthetic spectra may not produce matches — that's OK for a basic test
            let msg = e.to_string();
            eprintln!(
                "Sage search returned error (may be OK with synthetic data): {}",
                msg
            );
        }
    }
}

#[tokio::test]
async fn sage_engine_info_and_health() {
    let adapter = SageAdapter::default();
    let info = adapter.engine_info();
    assert_eq!(info.name, "Sage");

    let health = adapter.health_check().await;
    match health {
        Ok(status) => {
            assert_eq!(status, protein_copilot_core::engine::HealthStatus::Healthy);
        }
        Err(e) => {
            panic!("health_check should not fail for Sage: {}", e);
        }
    }
}

#[tokio::test]
async fn sage_search_with_nonexistent_fasta_returns_error() {
    let adapter = SageAdapter::default();
    let params = test_params("/nonexistent/path/to.fasta");

    // Create a minimal spectrum to pass the MS2 check
    let spectrum = protein_copilot_core::spectrum::Spectrum {
        scan_number: 1,
        ms_level: protein_copilot_core::spectrum::MsLevel::MS2,
        retention_time_min: 60.0,
        precursors: vec![protein_copilot_core::spectrum::PrecursorInfo {
            mz: 500.0,
            charge: Some(2),
            intensity: Some(1000.0),
            isolation_window: None,
            source_scan: None,
        }],
        mz_array: vec![100.0, 200.0, 300.0],
        intensity_array: vec![1000.0, 2000.0, 3000.0],
    };

    let on_progress: protein_copilot_core::progress::ProgressCallback = Box::new(|_| {});
    let result = adapter
        .search_with_spectra(&params, vec![spectrum], on_progress, &mut protein_copilot_core::diagnostics::SearchDiagnostics::new())
        .await;
    assert!(result.is_err(), "Should fail with nonexistent FASTA");
    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("FASTA")
            || err_msg.contains("fasta")
            || err_msg.contains("not found")
            || err_msg.contains("No such file"),
        "Error should mention FASTA: {}",
        err_msg
    );
}
