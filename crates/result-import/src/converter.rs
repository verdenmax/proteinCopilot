//! Converts matched ImportedPsms to core::SearchResult.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use protein_copilot_core::engine::EngineInfo;
use protein_copilot_core::run_metadata::RunMetadata;
use protein_copilot_core::search_params::{
    DecoyStrategy, Enzyme, MassTolerance, SearchParams, ToleranceUnit,
};
use protein_copilot_core::search_result::{
    PeptideResult, ProteinResult, Psm, SearchResult, SearchResultSummary,
};
use protein_copilot_search_engine::chemistry::{peptide_mass, peptide_mz};
use uuid::Uuid;

use crate::{ImportResult, ImportedPsm, MatchReport};

/// Convert matched ImportedPsms into a standard SearchResult.
///
/// Only PSMs with `matched_scan.is_some()` are included.
/// Returns `(SearchResult, ImportResult)`.
pub fn build_search_result(
    psms: &[ImportedPsm],
    match_report: MatchReport,
    format_name: &str,
    input_files: Vec<PathBuf>,
) -> (SearchResult, ImportResult) {
    let run_id = Uuid::new_v4();

    let core_psms: Vec<Psm> = psms
        .iter()
        .filter(|p| p.matched_scan.is_some())
        .map(|p| to_core_psm(p))
        .collect();

    let peptides = aggregate_peptides(&core_psms);
    let proteins = aggregate_proteins(&core_psms);

    let unique_peptide_count = core_psms
        .iter()
        .map(|p| p.peptide_sequence.as_str())
        .collect::<HashSet<_>>()
        .len();
    let unique_protein_count = core_psms
        .iter()
        .flat_map(|p| p.protein_accessions.iter().map(|a| a.as_str()))
        .collect::<HashSet<_>>()
        .len();

    let scores: Vec<f64> = core_psms.iter().map(|p| p.score).collect();
    let delta_ppms: Vec<f64> = core_psms.iter().map(|p| p.delta_mass_ppm).collect();
    let median_score = median(&scores);
    let median_delta = median(&delta_ppms);

    let mut charge_dist: HashMap<i32, u64> = HashMap::new();
    let mut mod_dist: HashMap<String, u64> = HashMap::new();
    for psm in &core_psms {
        *charge_dist.entry(psm.charge).or_default() += 1;
        for m in &psm.modifications {
            *mod_dist.entry(m.name.clone()).or_default() += 1;
        }
    }

    // PSMs at 1% FDR: if we have q-values, filter; otherwise count all
    let has_qvalues = core_psms.iter().any(|p| p.q_value.is_some());
    let psms_at_fdr = if has_qvalues {
        core_psms
            .iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .count() as u64
    } else {
        core_psms.len() as u64
    };
    let peptides_at_fdr = if has_qvalues {
        core_psms
            .iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .map(|p| p.peptide_sequence.as_str())
            .collect::<HashSet<_>>()
            .len() as u64
    } else {
        unique_peptide_count as u64
    };

    let summary = SearchResultSummary {
        total_spectra_searched: core_psms.len() as u64,
        total_psms: core_psms.len() as u64,
        psms_at_1pct_fdr: psms_at_fdr,
        unique_peptides_at_1pct_fdr: peptides_at_fdr,
        protein_groups_at_1pct_fdr: unique_protein_count as u64,
        median_score,
        median_delta_mass_ppm: median_delta,
        identification_rate: 1.0, // all imported PSMs are "identified"
        modification_distribution: mod_dist,
        charge_distribution: charge_dist,
        search_duration_sec: 0.0,
    };

    let engine_info = EngineInfo {
        name: "imported".to_string(),
        version: format_name.to_string(),
        supported_features: vec![],
    };

    // Placeholder SearchParams — external results don't carry search params
    let params = SearchParams {
        enzyme: Enzyme::Trypsin,
        missed_cleavages: 2,
        fixed_modifications: vec![],
        variable_modifications: vec![],
        precursor_tolerance: MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        },
        fragment_tolerance: MassTolerance {
            value: 0.02,
            unit: ToleranceUnit::Da,
        },
        database_path: "imported".to_string(),
        decoy_strategy: DecoyStrategy::None,
        acquisition_mode: None,
        max_variable_modifications: 3,
        min_peptide_length: 7,
        max_peptide_length: 50,
    };

    let mut metadata = RunMetadata::new(params.clone(), engine_info.clone(), input_files);
    metadata.run_id = run_id;

    let raw_files: Vec<String> = psms
        .iter()
        .map(|p| p.raw_name.clone())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let search_result = SearchResult {
        run_id,
        engine_info,
        params_used: params,
        psms: core_psms,
        peptides,
        proteins,
        summary,
        metadata,
    };

    let import_result = ImportResult {
        run_id: run_id.to_string(),
        match_report,
        imported_psm_count: search_result.psms.len(),
        unique_peptides: unique_peptide_count,
        protein_count: unique_protein_count,
        raw_files,
    };

    (search_result, import_result)
}

fn to_core_psm(imported: &ImportedPsm) -> Psm {
    let scan = imported.matched_scan.unwrap_or(0);

    // Calculate theoretical m/z from sequence + modification masses
    let mod_mass: f64 = imported.modifications.iter().map(|m| m.mass_delta).sum();
    let calculated_mz = peptide_mass(&imported.sequence)
        .map(|neutral| peptide_mz(neutral + mod_mass, imported.charge))
        .unwrap_or(imported.precursor_mz); // fallback for non-standard AAs

    let delta_ppm = if calculated_mz > 0.0 {
        (imported.precursor_mz - calculated_mz) / calculated_mz * 1e6
    } else {
        0.0
    };

    Psm {
        spectrum_scan: scan,
        peptide_sequence: imported.sequence.clone(),
        modifications: imported.modifications.clone(),
        charge: imported.charge,
        precursor_mz: imported.precursor_mz,
        calculated_mz,
        delta_mass_ppm: delta_ppm,
        score: imported.score.unwrap_or(0.0),
        q_value: imported.q_value,
        protein_accessions: imported.protein_accessions.clone(),
        is_decoy: false,
    }
}

fn aggregate_peptides(psms: &[Psm]) -> Vec<PeptideResult> {
    let mut map: HashMap<&str, (Vec<&str>, f64, u64, Option<f64>)> = HashMap::new();
    for psm in psms {
        let entry = map
            .entry(psm.peptide_sequence.as_str())
            .or_insert_with(|| (Vec::new(), f64::MIN, 0, None));
        for acc in &psm.protein_accessions {
            if !entry.0.contains(&acc.as_str()) {
                entry.0.push(acc.as_str());
            }
        }
        if psm.score > entry.1 {
            entry.1 = psm.score;
        }
        entry.2 += 1;
        if let Some(q) = psm.q_value {
            entry.3 = Some(entry.3.map_or(q, |prev: f64| prev.min(q)));
        }
    }

    map.into_iter()
        .map(|(seq, (proteins, best_score, count, q))| PeptideResult {
            sequence: seq.to_string(),
            protein_accessions: proteins.into_iter().map(|s| s.to_string()).collect(),
            best_score: if best_score == f64::MIN { 0.0 } else { best_score },
            q_value: q,
            psm_count: count,
        })
        .collect()
}

fn aggregate_proteins(psms: &[Psm]) -> Vec<ProteinResult> {
    let mut map: HashMap<&str, HashSet<&str>> = HashMap::new();
    for psm in psms {
        for acc in &psm.protein_accessions {
            map.entry(acc.as_str())
                .or_default()
                .insert(psm.peptide_sequence.as_str());
        }
    }

    map.into_iter()
        .map(|(acc, peptides)| ProteinResult {
            accession: acc.to_string(),
            description: String::new(),
            coverage: 0.0, // cannot compute without protein sequence
            peptide_count: peptides.len() as u64,
            unique_peptide_count: peptides.len() as u64,
        })
        .collect()
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    sorted[sorted.len() / 2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{ModPosition, Modification};

    fn sample_psms() -> Vec<ImportedPsm> {
        vec![
            ImportedPsm {
                sequence: "PEPTIDE".to_string(),
                charge: 2,
                precursor_mz: 400.19,
                rt_sec: 600.0,
                modifications: vec![],
                score: Some(0.001),
                q_value: Some(0.001),
                protein_accessions: vec!["P12345".to_string()],
                raw_name: "test".to_string(),
                matched_scan: Some(100),
                rt_delta_sec: Some(1.5),
            },
            ImportedPsm {
                sequence: "PEPTMIDE".to_string(),
                charge: 3,
                precursor_mz: 310.5,
                rt_sec: 1200.0,
                modifications: vec![Modification {
                    name: "Oxidation".to_string(),
                    mass_delta: 15.994915,
                    residues: vec!['M'],
                    position: ModPosition::Anywhere,
                }],
                score: Some(0.005),
                q_value: Some(0.005),
                protein_accessions: vec!["P12345".to_string()],
                raw_name: "test".to_string(),
                matched_scan: Some(200),
                rt_delta_sec: Some(0.5),
            },
            ImportedPsm {
                sequence: "PEPTIDE".to_string(),
                charge: 2,
                precursor_mz: 400.19,
                rt_sec: 1800.0,
                modifications: vec![],
                score: None,
                q_value: None,
                protein_accessions: vec!["P67890".to_string()],
                raw_name: "test".to_string(),
                matched_scan: None, // unmatched — should be excluded
                rt_delta_sec: None,
            },
        ]
    }

    #[test]
    fn build_search_result_excludes_unmatched() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3,
            matched: 2,
            unmatched: 1,
            median_rt_delta_sec: 1.0,
            max_rt_delta_sec: 1.5,
            per_file: HashMap::new(),
        };
        let (result, import) = build_search_result(&psms, report, "custom_json", vec![]);
        assert_eq!(result.psms.len(), 2);
        assert_eq!(import.imported_psm_count, 2);
    }

    #[test]
    fn build_search_result_calculates_mz() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3,
            matched: 2,
            unmatched: 1,
            median_rt_delta_sec: 0.0,
            max_rt_delta_sec: 0.0,
            per_file: HashMap::new(),
        };
        let (result, _) = build_search_result(&psms, report, "custom_json", vec![]);
        for psm in &result.psms {
            assert!(psm.calculated_mz > 0.0, "calculated_mz should be positive");
            assert!(
                psm.delta_mass_ppm.is_finite(),
                "delta_mass_ppm should be finite"
            );
        }
    }

    #[test]
    fn build_search_result_aggregates_peptides() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3,
            matched: 2,
            unmatched: 1,
            median_rt_delta_sec: 0.0,
            max_rt_delta_sec: 0.0,
            per_file: HashMap::new(),
        };
        let (result, _) = build_search_result(&psms, report, "custom_json", vec![]);
        assert_eq!(result.peptides.len(), 2); // PEPTIDE and PEPTMIDE
    }

    #[test]
    fn build_search_result_engine_info() {
        let psms = sample_psms();
        let report = MatchReport {
            total_psms: 3,
            matched: 2,
            unmatched: 1,
            median_rt_delta_sec: 0.0,
            max_rt_delta_sec: 0.0,
            per_file: HashMap::new(),
        };
        let (result, _) = build_search_result(&psms, report, "diann_parquet", vec![]);
        assert_eq!(result.engine_info.name, "imported");
        assert_eq!(result.engine_info.version, "diann_parquet");
    }
}
