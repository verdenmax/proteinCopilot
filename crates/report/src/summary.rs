//! Summary generation — statistical aggregation of search results.

use std::collections::{HashMap, HashSet};

use protein_copilot_core::search_result::{SearchResult, SearchResultSummary};
use protein_copilot_core::util::compute_median;

/// Generates a statistical summary with optional FDR filtering.
pub(crate) fn generate_summary(result: &SearchResult) -> SearchResultSummary {
    let _span = tracing::info_span!("generate_summary").entered();

    let total_spectra = result.summary.total_spectra_searched;
    let total_psms = result.psms.len() as u64;

    // Filter PSMs at 1% FDR (q_value ≤ 0.01), or keep all if no q-values
    let has_qvalues = result.psms.iter().any(|p| p.q_value.is_some());
    let filtered_psms: Vec<_> = if has_qvalues {
        result
            .psms
            .iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .collect()
    } else {
        result.psms.iter().collect()
    };

    let psms_at_1pct_fdr = filtered_psms.len() as u64;

    // Unique peptides and proteins at 1% FDR
    let unique_peptides: HashSet<&str> = filtered_psms
        .iter()
        .map(|p| p.peptide_sequence.as_str())
        .collect();
    let unique_proteins: HashSet<&str> = filtered_psms
        .iter()
        .flat_map(|p| p.protein_accessions.iter().map(|a| a.as_str()))
        .collect();

    // Identification rate
    let identification_rate = if total_spectra > 0 {
        psms_at_1pct_fdr as f64 / total_spectra as f64
    } else {
        0.0
    };

    // Modification and charge distributions (from filtered PSMs)
    let mut mod_dist: HashMap<String, u64> = HashMap::new();
    let mut charge_dist: HashMap<i32, u64> = HashMap::new();
    let mut scores: Vec<f64> = Vec::new();
    let mut delta_ppms: Vec<f64> = Vec::new();

    for psm in &filtered_psms {
        for m in &psm.modifications {
            *mod_dist.entry(m.name.clone()).or_insert(0) += 1;
        }
        if psm.charge > 0 {
            *charge_dist.entry(psm.charge).or_insert(0) += 1;
        }
        if psm.score.is_finite() {
            scores.push(psm.score);
        }
        if psm.delta_mass_ppm.is_finite() {
            delta_ppms.push(psm.delta_mass_ppm);
        }
    }

    // Median calculations
    scores.sort_by(|a, b| a.total_cmp(b));
    delta_ppms.sort_by(|a, b| a.total_cmp(b));

    let median_score = compute_median(&scores);
    let median_delta = compute_median(&delta_ppms);

    tracing::info!(
        id_rate = format!("{:.1}%", identification_rate * 100.0),
        psms_1pct = psms_at_1pct_fdr,
        "summary generated"
    );

    if identification_rate < 0.10 {
        tracing::warn!(
            id_rate = format!("{:.1}%", identification_rate * 100.0),
            psms_1pct = psms_at_1pct_fdr,
            total_spectra = total_spectra,
            "low identification rate — check search parameters or database"
        );
    }

    SearchResultSummary {
        total_spectra_searched: total_spectra,
        total_psms,
        psms_at_1pct_fdr,
        unique_peptides_at_1pct_fdr: unique_peptides.len() as u64,
        protein_groups_at_1pct_fdr: unique_proteins.len() as u64,
        median_score,
        median_delta_mass_ppm: median_delta,
        identification_rate,
        modification_distribution: mod_dist,
        charge_distribution: charge_dist,
        search_duration_sec: result.summary.search_duration_sec,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::engine::EngineInfo;
    use protein_copilot_core::run_metadata::RunMetadata;
    use protein_copilot_core::search_params::*;
    use protein_copilot_core::search_result::*;
    use std::path::PathBuf;

    fn sample_params() -> SearchParams {
        SearchParams {
            database_path: "/data/db.fasta".to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 2,
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
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: None,
        }
    }

    fn sample_engine_info() -> EngineInfo {
        EngineInfo {
            name: "SimpleSearch".to_string(),
            version: "0.1.0".to_string(),
            supported_features: vec![],
        }
    }

    fn sample_result_with_qvalues() -> SearchResult {
        let meta = RunMetadata::new(
            sample_params(),
            sample_engine_info(),
            vec![PathBuf::from("/data/test.mgf")],
        );
        SearchResult {
            run_id: meta.run_id,
            engine_info: meta.engine_info.clone(),
            params_used: meta.params_used.clone(),
            psms: vec![
                Psm {
                    spectrum_scan: 1,
                    peptide_sequence: "PEPTIDER".to_string(),
                    modifications: vec![],
                    charge: 2,
                    precursor_mz: 471.25,
                    calculated_mz: 471.25,
                    delta_mass_ppm: 0.5,
                    score: 0.8,
                    q_value: Some(0.001),
                    protein_accessions: vec!["P001".to_string()],
                    is_decoy: false,
                    extra: None,
                },
                Psm {
                    spectrum_scan: 2,
                    peptide_sequence: "ANOTHERR".to_string(),
                    modifications: vec![],
                    charge: 3,
                    precursor_mz: 300.0,
                    calculated_mz: 300.0,
                    delta_mass_ppm: 1.0,
                    score: 0.6,
                    q_value: Some(0.005),
                    protein_accessions: vec!["P001".to_string()],
                    is_decoy: false,
                    extra: None,
                },
                Psm {
                    spectrum_scan: 3,
                    peptide_sequence: "DECOYSEQ".to_string(),
                    modifications: vec![],
                    charge: 2,
                    precursor_mz: 400.0,
                    calculated_mz: 400.0,
                    delta_mass_ppm: 2.0,
                    score: 0.3,
                    q_value: Some(0.05),
                    protein_accessions: vec!["REV_P002".to_string()],
                    is_decoy: true,
                    extra: None,
                },
            ],
            peptides: vec![],
            proteins: vec![],
            summary: SearchResultSummary {
                total_spectra_searched: 10,
                total_psms: 3,
                psms_at_1pct_fdr: 2,
                unique_peptides_at_1pct_fdr: 2,
                protein_groups_at_1pct_fdr: 1,
                median_score: 0.6,
                median_delta_mass_ppm: 1.0,
                identification_rate: 0.2,
                modification_distribution: HashMap::new(),
                charge_distribution: HashMap::new(),
                search_duration_sec: 1.0,
            },
            metadata: meta,
        }
    }

    #[test]
    fn summary_filters_by_fdr() {
        let result = sample_result_with_qvalues();
        let summary = generate_summary(&result);
        // PSM 3 has q_value=0.05 > 0.01, should be filtered out
        assert_eq!(summary.psms_at_1pct_fdr, 2);
        assert_eq!(summary.total_psms, 3);
    }

    #[test]
    fn summary_unique_peptides_at_fdr() {
        let result = sample_result_with_qvalues();
        let summary = generate_summary(&result);
        assert_eq!(summary.unique_peptides_at_1pct_fdr, 2); // PEPTIDER + ANOTHERR
    }

    #[test]
    fn summary_identification_rate() {
        let result = sample_result_with_qvalues();
        let summary = generate_summary(&result);
        // 2 PSMs at FDR / 10 spectra = 0.2
        assert!((summary.identification_rate - 0.2).abs() < 0.01);
    }

    #[test]
    fn summary_median_score() {
        let result = sample_result_with_qvalues();
        let summary = generate_summary(&result);
        // Filtered PSMs scores: [0.6, 0.8], proper median = (0.6+0.8)/2 = 0.7
        assert!((summary.median_score - 0.7).abs() < 0.01);
    }

    #[test]
    fn summary_charge_distribution() {
        let result = sample_result_with_qvalues();
        let summary = generate_summary(&result);
        // Filtered: charge 2 (PSM1) + charge 3 (PSM2)
        assert_eq!(*summary.charge_distribution.get(&2).unwrap_or(&0), 1);
        assert_eq!(*summary.charge_distribution.get(&3).unwrap_or(&0), 1);
    }

    #[test]
    fn summary_no_qvalues_keeps_all() {
        let mut result = sample_result_with_qvalues();
        for psm in &mut result.psms {
            psm.q_value = None;
        }
        let summary = generate_summary(&result);
        assert_eq!(summary.psms_at_1pct_fdr, 3); // all kept
    }
}
