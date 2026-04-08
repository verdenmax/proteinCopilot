//! Simplified built-in search engine for MVP validation.
//!
//! Implements the full `SearchEngineAdapter` trait using:
//! - FASTA parsing → in-silico digestion → precursor matching → b/y scoring
//!
//! This is a placeholder engine that validates the architecture and data flow.
//! It should be replaced by pFind/MSFragger adapters for production use.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use protein_copilot_core::engine::{EngineInfo, HealthStatus, SearchEngineAdapter};
use protein_copilot_core::error::CoreError;
use protein_copilot_core::progress::{ProgressCallback, SearchProgress};
use protein_copilot_core::run_metadata::{RunMetadata, RunStatus};
use protein_copilot_core::search_params::{DecoyStrategy, SearchParams};
use protein_copilot_core::search_result::{
    PeptideResult, ProteinResult, Psm, SearchResult, SearchResultSummary,
};
use protein_copilot_core::spectrum::{MsLevel, Spectrum};
use uuid::Uuid;

use crate::digest::{digest_with_length, DigestedPeptide};
use crate::error::SearchEngineError;
use crate::fasta::parse_fasta;
use crate::matching::{match_spectrum, match_spectrum_all, PeptideMatch};

/// A simplified search engine that runs entirely in-process.
///
/// Performs: FASTA → digest → precursor match → b/y score → SearchResult.
/// No SSH, no external binaries, no statistical scoring model.
pub struct SimpleSearchEngine;

impl SimpleSearchEngine {
    /// Creates a new SimpleSearchEngine instance.
    pub fn new() -> Self {
        Self
    }

    /// Core search logic (synchronous, called from async trait method).
    fn run_search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: &dyn Fn(SearchProgress),
    ) -> Result<SearchResult, SearchEngineError> {
        let start = Instant::now();
        let run_id = Uuid::new_v4();

        // Helper to report progress at a given stage
        let report = |stage: &str, pct: f64| {
            on_progress(SearchProgress {
                run_id,
                status: "Running".to_string(),
                stage: Some(stage.to_string()),
                progress_pct: Some(pct),
                elapsed_sec: start.elapsed().as_secs_f64(),
                estimated_remaining_sec: None,
            });
        };

        // Step 1: Validate parameters
        params
            .validate()
            .map_err(|e| SearchEngineError::InvalidParams {
                detail: e.to_string(),
            })?;

        if input_files.is_empty() {
            return Err(SearchEngineError::NoInputSpectra);
        }

        // Step 2: Read FASTA database and digest
        report("Reading FASTA database", 0.02);
        let fasta_path = Path::new(&params.database_path);
        let proteins = parse_fasta(fasta_path)?;

        report("Digesting proteins", 0.08);
        let mut all_peptides: Vec<DigestedPeptide> = Vec::new();
        for protein in &proteins {
            let peptides = digest_with_length(
                &protein.sequence,
                &protein.accession,
                &params.enzyme,
                params.missed_cleavages,
                params.min_peptide_length,
                params.max_peptide_length,
            );
            all_peptides.extend(peptides);
        }

        if all_peptides.is_empty() {
            return Err(SearchEngineError::ExecutionError {
                detail: format!(
                    "no candidate peptides generated from {} proteins \
                     (all peptides may be shorter than {} or longer than {} residues)",
                    proteins.len(),
                    params.min_peptide_length,
                    params.max_peptide_length,
                ),
            });
        }

        // Generate and digest decoy database if strategy is active
        if params.decoy_strategy != DecoyStrategy::None {
            report("Generating decoy database", 0.10);
            let target_tuples: Vec<(String, String, String)> = proteins
                .iter()
                .map(|p| (p.accession.clone(), p.description.clone(), p.sequence.clone()))
                .collect();
            let decoy_proteins =
                protein_copilot_fdr::generate_decoys(&target_tuples, params.decoy_strategy);
            for decoy in &decoy_proteins {
                let peptides = digest_with_length(
                    &decoy.sequence,
                    &decoy.accession,
                    &params.enzyme,
                    params.missed_cleavages,
                    params.min_peptide_length,
                    params.max_peptide_length,
                );
                all_peptides.extend(peptides);
            }
        }

        // Step 3: Read all spectra from input files
        report("Reading spectra", 0.15);
        let mut all_spectra: Vec<Spectrum> = Vec::new();
        for file_path in input_files {
            let info = protein_copilot_spectrum_io::detect_format(file_path).map_err(|e| {
                SearchEngineError::IoError {
                    detail: e.to_string(),
                }
            })?;
            let reader = protein_copilot_spectrum_io::create_reader(&info);
            let spectra = reader
                .read_all(file_path)
                .map_err(|e| SearchEngineError::IoError {
                    detail: e.to_string(),
                })?;
            all_spectra.extend(spectra);
        }

        if all_spectra.is_empty() {
            return Err(SearchEngineError::NoInputSpectra);
        }

        self.run_search_on_spectra(
            params,
            all_spectra,
            all_peptides,
            &proteins,
            input_files,
            run_id,
            start,
            on_progress,
        )
    }

    /// Core search logic operating on pre-loaded spectra.
    ///
    /// Shared by `run_search` (file-based) and `search_with_spectra` (pre-loaded).
    #[allow(clippy::too_many_arguments)]
    fn run_search_on_spectra(
        &self,
        params: &SearchParams,
        all_spectra: Vec<Spectrum>,
        all_peptides: Vec<DigestedPeptide>,
        proteins: &[crate::fasta::FastaEntry],
        input_files: &[PathBuf],
        run_id: Uuid,
        start: Instant,
        on_progress: &dyn Fn(SearchProgress),
    ) -> Result<SearchResult, SearchEngineError> {
        let report = |stage: &str, pct: f64| {
            on_progress(SearchProgress {
                run_id,
                status: "Running".to_string(),
                stage: Some(stage.to_string()),
                progress_pct: Some(pct),
                elapsed_sec: start.elapsed().as_secs_f64(),
                estimated_remaining_sec: None,
            });
        };

        // Filter to MS2 only (MS1 survey scans have no precursors to match)
        let ms2_spectra: Vec<&Spectrum> = all_spectra
            .iter()
            .filter(|s| s.ms_level == MsLevel::MS2)
            .collect();
        if ms2_spectra.is_empty() {
            return Err(SearchEngineError::NoInputSpectra);
        }
        let total_spectra = ms2_spectra.len();
        let mut psms: Vec<Psm> = Vec::new();

        for (i, spectrum) in ms2_spectra.iter().enumerate() {
            if i % 50 == 0 || i + 1 == total_spectra {
                let pct = 0.15 + 0.75 * (i as f64 / total_spectra.max(1) as f64);
                report(
                    &format!("Matching spectra ({}/{})", i + 1, total_spectra),
                    pct,
                );
            }

            if spectrum.precursors.len() > 1 {
                // DIA mode: multiple precursors, collect all matches
                let matches = match_spectrum_all(
                    spectrum,
                    &all_peptides,
                    &params.precursor_tolerance,
                    &params.fragment_tolerance,
                    &params.fixed_modifications,
                    &params.variable_modifications,
                    params.max_variable_modifications,
                );
                for m in &matches {
                    psms.push(build_psm(spectrum, m, &params.fixed_modifications));
                }
            } else {
                // DDA mode: single precursor, use original function
                if let Some(m) = match_spectrum(
                    spectrum,
                    &all_peptides,
                    &params.precursor_tolerance,
                    &params.fragment_tolerance,
                    &params.fixed_modifications,
                    &params.variable_modifications,
                    params.max_variable_modifications,
                ) {
                    psms.push(build_psm(spectrum, &m, &params.fixed_modifications));
                }
            }
        }

        // Calculate FDR and assign q-values if decoy strategy is active
        if params.decoy_strategy != DecoyStrategy::None {
            report("Calculating FDR", 0.88);
            let scored: Vec<protein_copilot_fdr::calculation::ScoredPsm> = psms
                .iter()
                .enumerate()
                .map(|(i, p)| protein_copilot_fdr::calculation::ScoredPsm {
                    index: i,
                    score: p.score,
                    is_decoy: p.is_decoy,
                })
                .collect();

            if let Ok(qvalues) = protein_copilot_fdr::calculate_fdr(&scored) {
                for (idx, q) in qvalues {
                    if idx < psms.len() {
                        psms[idx].q_value = Some(q);
                    }
                }
            } else {
                tracing::warn!(
                    "FDR calculation failed (likely no decoy hits); q-values not assigned"
                );
            }

            // Remove decoy PSMs from final results
            psms.retain(|p| !p.is_decoy);
        }

        // Aggregate peptide and protein level results
        report("Aggregating results", 0.92);
        let peptides = aggregate_peptides(&psms);
        let protein_results = aggregate_proteins(&psms, proteins);

        // Build summary
        let duration = start.elapsed().as_secs_f64();
        let summary = build_summary(&psms, ms2_spectra.len() as u64, duration);

        // Build metadata
        let engine_info = self.engine_info();
        let mut metadata =
            RunMetadata::new(params.clone(), engine_info.clone(), input_files.to_vec());
        metadata.run_id = run_id;
        metadata.status = RunStatus::Completed;
        metadata.duration_sec = Some(duration);

        Ok(SearchResult {
            run_id,
            engine_info,
            params_used: params.clone(),
            psms,
            peptides,
            proteins: protein_results,
            summary,
            metadata,
        })
    }
}

impl Default for SimpleSearchEngine {
    fn default() -> Self {
        Self::new()
    }
}

#[async_trait::async_trait]
impl SearchEngineAdapter for SimpleSearchEngine {
    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        self.run_search(params, input_files, &*on_progress)
            .map_err(CoreError::from)
    }

    fn engine_info(&self) -> EngineInfo {
        EngineInfo {
            name: "SimpleSearch".to_string(),
            version: "0.1.0".to_string(),
            supported_features: vec!["basic_search".to_string(), "b_y_scoring".to_string()],
        }
    }

    async fn health_check(&self) -> Result<HealthStatus, CoreError> {
        Ok(HealthStatus::Healthy)
    }

    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        let start = Instant::now();
        let run_id = Uuid::new_v4();

        let progress_fn = move |p: SearchProgress| {
            on_progress(p);
        };

        let report = |stage: &str, pct: f64| {
            progress_fn(SearchProgress {
                run_id,
                status: "Running".to_string(),
                stage: Some(stage.to_string()),
                progress_pct: Some(pct),
                elapsed_sec: start.elapsed().as_secs_f64(),
                estimated_remaining_sec: None,
            });
        };

        // Validate params
        params.validate().map_err(|e| {
            CoreError::from(SearchEngineError::InvalidParams {
                detail: e.to_string(),
            })
        })?;

        if spectra.is_empty() {
            return Err(CoreError::from(SearchEngineError::NoInputSpectra));
        }

        // Digest database (same as run_search)
        report("Reading FASTA database", 0.02);
        let fasta_path = Path::new(&params.database_path);
        let proteins = parse_fasta(fasta_path).map_err(CoreError::from)?;

        report("Digesting proteins", 0.08);
        let mut all_peptides: Vec<DigestedPeptide> = Vec::new();
        for protein in &proteins {
            let peptides = digest_with_length(
                &protein.sequence,
                &protein.accession,
                &params.enzyme,
                params.missed_cleavages,
                params.min_peptide_length,
                params.max_peptide_length,
            );
            all_peptides.extend(peptides);
        }

        if all_peptides.is_empty() {
            return Err(CoreError::from(SearchEngineError::ExecutionError {
                detail: format!(
                    "no candidate peptides generated from {} proteins \
                     (all peptides may be shorter than {} or longer than {} residues)",
                    proteins.len(),
                    params.min_peptide_length,
                    params.max_peptide_length,
                ),
            }));
        }

        // Generate and digest decoy database if strategy is active
        if params.decoy_strategy != DecoyStrategy::None {
            report("Generating decoy database", 0.10);
            let target_tuples: Vec<(String, String, String)> = proteins
                .iter()
                .map(|p| (p.accession.clone(), p.description.clone(), p.sequence.clone()))
                .collect();
            let decoy_proteins =
                protein_copilot_fdr::generate_decoys(&target_tuples, params.decoy_strategy);
            for decoy in &decoy_proteins {
                let peptides = digest_with_length(
                    &decoy.sequence,
                    &decoy.accession,
                    &params.enzyme,
                    params.missed_cleavages,
                    params.min_peptide_length,
                    params.max_peptide_length,
                );
                all_peptides.extend(peptides);
            }
        }

        // Use the shared core logic (no input files for spectra-based search)
        self.run_search_on_spectra(
            params,
            spectra,
            all_peptides,
            &proteins,
            &[],
            run_id,
            start,
            &progress_fn,
        )
        .map_err(CoreError::from)
    }
}

// ---------------------------------------------------------------------------
// Result building helpers
// ---------------------------------------------------------------------------

fn build_psm(
    spectrum: &Spectrum,
    m: &PeptideMatch,
    fixed_mods: &[protein_copilot_core::search_params::Modification],
) -> Psm {
    // Collect fixed modifications that apply to this peptide
    let mut mods = Vec::new();
    for fm in fixed_mods {
        for ch in m.peptide.sequence.chars() {
            if fm.residues.contains(&ch) {
                mods.push(fm.clone());
                break; // one per mod type
            }
        }
    }

    // Add variable modifications from the match
    for (vm, _pos) in &m.applied_variable_mods {
        mods.push(vm.clone());
    }

    let is_decoy = m.peptide.protein_accession.starts_with("REV_")
        || m.peptide.protein_accession.starts_with("SHUF_");

    Psm {
        spectrum_scan: spectrum.scan_number,
        peptide_sequence: m.peptide.sequence.clone(),
        modifications: mods,
        charge: m.charge,
        precursor_mz: m.observed_mz,
        calculated_mz: m.theoretical_mz,
        delta_mass_ppm: m.delta_mass_ppm,
        score: m.score,
        q_value: None,
        protein_accessions: vec![m.peptide.protein_accession.clone()],
        is_decoy,
    }
}

fn aggregate_peptides(psms: &[Psm]) -> Vec<PeptideResult> {
    let mut peptide_map: HashMap<String, (f64, Vec<String>, u64)> = HashMap::new();

    for psm in psms {
        let entry = peptide_map
            .entry(psm.peptide_sequence.clone())
            .or_insert_with(|| (f64::MIN, Vec::new(), 0));
        if psm.score > entry.0 {
            entry.0 = psm.score;
        }
        for acc in &psm.protein_accessions {
            if !entry.1.contains(acc) {
                entry.1.push(acc.clone());
            }
        }
        entry.2 += 1;
    }

    peptide_map
        .into_iter()
        .map(|(seq, (best_score, accessions, count))| PeptideResult {
            sequence: seq,
            protein_accessions: accessions,
            best_score,
            q_value: None,
            psm_count: count,
        })
        .collect()
}

fn aggregate_proteins(psms: &[Psm], proteins: &[crate::fasta::FastaEntry]) -> Vec<ProteinResult> {
    // Count peptides per protein
    let mut protein_peptides: HashMap<String, Vec<String>> = HashMap::new();
    for psm in psms {
        for acc in &psm.protein_accessions {
            let peptides = protein_peptides.entry(acc.clone()).or_default();
            if !peptides.contains(&psm.peptide_sequence) {
                peptides.push(psm.peptide_sequence.clone());
            }
        }
    }

    protein_peptides
        .into_iter()
        .map(|(acc, peptides)| {
            let description = proteins
                .iter()
                .find(|p| p.accession == acc)
                .map(|p| p.description.clone())
                .unwrap_or_default();

            // Coverage: track which positions are covered by identified peptides
            let protein_seq = proteins
                .iter()
                .find(|p| p.accession == acc)
                .map(|p| p.sequence.as_str())
                .unwrap_or("");
            let protein_len = protein_seq.len().max(1);

            let mut covered = vec![false; protein_len];
            for pep in &peptides {
                // Find all occurrences of this peptide in the protein
                let mut start = 0;
                while let Some(pos) = protein_seq[start..].find(pep.as_str()) {
                    let abs_pos = start + pos;
                    for i in abs_pos..abs_pos + pep.len() {
                        if i < covered.len() {
                            covered[i] = true;
                        }
                    }
                    start = abs_pos + 1;
                }
            }
            let covered_count = covered.iter().filter(|&&c| c).count();
            let coverage = covered_count as f64 / protein_len as f64;

            let peptide_count = peptides.len() as u64;

            ProteinResult {
                accession: acc,
                description,
                coverage,
                peptide_count,
                unique_peptide_count: peptide_count, // simplified: all unique
            }
        })
        .collect()
}

fn build_summary(psms: &[Psm], total_spectra: u64, duration: f64) -> SearchResultSummary {
    let total_psms = psms.len() as u64;

    // FDR-filtered counts: if q-values are available, count only PSMs at ≤1% FDR
    let has_qvalues = psms.iter().any(|p| p.q_value.is_some());
    let fdr_filtered: Vec<&Psm> = if has_qvalues {
        psms.iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .collect()
    } else {
        psms.iter().collect()
    };

    let psms_at_fdr = fdr_filtered.len() as u64;
    let unique_peptides_fdr: std::collections::HashSet<&str> = fdr_filtered
        .iter()
        .map(|p| p.peptide_sequence.as_str())
        .collect();
    let unique_proteins_fdr: std::collections::HashSet<&str> = fdr_filtered
        .iter()
        .flat_map(|p| p.protein_accessions.iter().map(|a| a.as_str()))
        .collect();

    let mut mod_dist: HashMap<String, u64> = HashMap::new();
    let mut charge_dist: HashMap<i32, u64> = HashMap::new();
    let mut scores: Vec<f64> = Vec::new();
    let mut delta_ppms: Vec<f64> = Vec::new();

    for psm in psms {
        for m in &psm.modifications {
            *mod_dist.entry(m.name.clone()).or_insert(0) += 1;
        }
        *charge_dist.entry(psm.charge).or_insert(0) += 1;
        scores.push(psm.score);
        delta_ppms.push(psm.delta_mass_ppm);
    }

    // Filter out non-finite values before sorting (defensive)
    scores.retain(|v| v.is_finite());
    delta_ppms.retain(|v| v.is_finite());

    scores.sort_by(|a, b| a.total_cmp(b));
    delta_ppms.sort_by(|a, b| a.total_cmp(b));

    let median_score = protein_copilot_core::util::compute_median(&scores);
    let median_delta = protein_copilot_core::util::compute_median(&delta_ppms);
    let id_rate = if total_spectra > 0 {
        psms_at_fdr as f64 / total_spectra as f64
    } else {
        0.0
    };

    SearchResultSummary {
        total_spectra_searched: total_spectra,
        total_psms,
        psms_at_1pct_fdr: psms_at_fdr,
        unique_peptides_at_1pct_fdr: unique_peptides_fdr.len() as u64,
        protein_groups_at_1pct_fdr: unique_proteins_fdr.len() as u64,
        median_score,
        median_delta_mass_ppm: median_delta,
        identification_rate: id_rate,
        modification_distribution: mod_dist,
        charge_distribution: charge_dist,
        search_duration_sec: duration,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::progress::noop_progress;
    use protein_copilot_core::search_params::{
        DecoyStrategy, Enzyme, MassTolerance, ToleranceUnit,
    };
    use std::io::Write;

    fn create_test_fasta() -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(
            f,
            ">sp|P001|TEST1 Test protein 1\n\
             PEPTIDEKANSTHERRLASTR\n\
             >sp|P002|TEST2 Test protein 2\n\
             AVCDEFGKHIKLMNPQRST\n"
        )
        .unwrap();
        f
    }

    fn test_params(fasta_path: &str) -> SearchParams {
        SearchParams {
            database_path: fasta_path.to_string(),
            enzyme: Enzyme::Trypsin,
            missed_cleavages: 1,
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
        }
    }

    #[test]
    fn engine_info_correct() {
        let engine = SimpleSearchEngine::new();
        let info = engine.engine_info();
        assert_eq!(info.name, "SimpleSearch");
        assert!(!info.version.is_empty());
    }

    #[tokio::test]
    async fn health_check_returns_healthy() {
        let engine = SimpleSearchEngine::new();
        let status = engine.health_check().await.unwrap();
        assert_eq!(status, HealthStatus::Healthy);
    }

    #[test]
    fn search_with_no_files_errors() {
        let fasta = create_test_fasta();
        let params = test_params(&fasta.path().to_string_lossy());
        let engine = SimpleSearchEngine::new();
        let result = engine.run_search(&params, &[], &|_| {});
        assert!(result.is_err());
    }

    #[test]
    fn search_with_invalid_params_errors() {
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
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
        };
        let engine = SimpleSearchEngine::new();
        let result = engine.run_search(&params, &[PathBuf::from("test.mgf")], &|_| {});
        assert!(result.is_err());
    }

    #[test]
    fn search_produces_valid_result() {
        let fasta = create_test_fasta();
        let params = test_params(&fasta.path().to_string_lossy());

        // Use the MGF fixture from spectrum-io
        let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("spectrum-io")
            .join("tests")
            .join("fixtures")
            .join("small.mgf");

        let engine = SimpleSearchEngine::new();
        let result = engine.run_search(&params, &[fixture], &|_| {}).unwrap();

        // Basic structural checks
        assert_eq!(result.engine_info.name, "SimpleSearch");
        assert_eq!(result.params_used.enzyme, Enzyme::Trypsin);
        assert_eq!(result.summary.total_spectra_searched, 10);
        assert!(result.summary.search_duration_sec >= 0.0);
        assert!(result.summary.identification_rate >= 0.0);
        assert!(result.summary.identification_rate <= 1.0);
        assert_eq!(result.run_id, result.metadata.run_id);
        assert_eq!(result.engine_info, result.metadata.engine_info);

        // Metadata should be complete
        assert_eq!(result.metadata.status, RunStatus::Completed);
        assert!(result.metadata.duration_sec.is_some());
    }

    #[tokio::test]
    async fn search_via_trait() {
        let fasta = create_test_fasta();
        let params = test_params(&fasta.path().to_string_lossy());
        let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("spectrum-io")
            .join("tests")
            .join("fixtures")
            .join("small.mgf");

        let engine = SimpleSearchEngine::new();
        let result = engine
            .search(&params, &[fixture], noop_progress())
            .await
            .unwrap();
        assert_eq!(result.engine_info.name, "SimpleSearch");
    }

    #[test]
    fn registry_with_simple_engine() {
        let mut registry = crate::registry::EngineRegistry::new();
        registry.register(Box::new(SimpleSearchEngine::new()));

        assert_eq!(registry.len(), 1);
        assert!(registry.get("SimpleSearch").is_some());
        assert_eq!(registry.list_available()[0].name, "SimpleSearch");
    }

    #[test]
    fn search_with_short_proteins_errors() {
        // All proteins too short to produce peptides (< 6 aa after digestion)
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(f, ">tiny\nACK\n>tiny2\nGR\n").unwrap();

        let params = test_params(&f.path().to_string_lossy());
        let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("spectrum-io")
            .join("tests")
            .join("fixtures")
            .join("small.mgf");

        let engine = SimpleSearchEngine::new();
        let result = engine.run_search(&params, &[fixture], &|_| {});
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("no candidate peptides"));
    }

    #[tokio::test]
    async fn search_reports_progress_stages() {
        use std::sync::{Arc, Mutex};

        let stages: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let stages_clone = Arc::clone(&stages);

        let on_progress: ProgressCallback = Box::new(move |p: SearchProgress| {
            if let Some(stage) = p.stage {
                let mut s = stages_clone.lock().unwrap();
                // Only push if stage changed
                if s.last().map(|l| l != &stage).unwrap_or(true) {
                    s.push(stage);
                }
            }
        });

        let fasta = create_test_fasta();
        let params = test_params(&fasta.path().to_string_lossy());
        let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("spectrum-io")
            .join("tests")
            .join("fixtures")
            .join("small.mgf");

        let engine = SimpleSearchEngine::new();
        let _result = engine
            .search(&params, &[fixture], on_progress)
            .await
            .unwrap();

        let recorded = stages.lock().unwrap();
        assert!(
            recorded.len() >= 4,
            "Expected at least 4 stages, got: {recorded:?}"
        );
        assert!(
            recorded[0].contains("FASTA"),
            "First stage should be FASTA reading, got: {}",
            recorded[0]
        );
        assert!(recorded.iter().any(|s| s.contains("Matching")));
        assert!(recorded.last().unwrap().contains("Aggregating"));
    }
}
