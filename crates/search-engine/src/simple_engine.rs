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
use protein_copilot_core::progress::ProgressCallback;
use protein_copilot_core::run_metadata::{RunMetadata, RunStatus};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{
    PeptideResult, ProteinResult, Psm, SearchResult, SearchResultSummary,
};
use protein_copilot_core::spectrum::Spectrum;

use crate::digest::{digest, DigestedPeptide};
use crate::error::SearchEngineError;
use crate::fasta::parse_fasta;
use crate::matching::{match_spectrum, PeptideMatch};

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
    ) -> Result<SearchResult, SearchEngineError> {
        let start = Instant::now();

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
        let fasta_path = Path::new(&params.database_path);
        let proteins = parse_fasta(fasta_path)?;

        let mut all_peptides: Vec<DigestedPeptide> = Vec::new();
        for protein in &proteins {
            let peptides = digest(
                &protein.sequence,
                &protein.accession,
                &params.enzyme,
                params.missed_cleavages,
            );
            all_peptides.extend(peptides);
        }

        if all_peptides.is_empty() {
            return Err(SearchEngineError::ExecutionError {
                detail: format!(
                    "no candidate peptides generated from {} proteins \
                     (all peptides may be shorter than 6 or longer than 50 residues)",
                    proteins.len()
                ),
            });
        }

        // Step 3: Read all spectra from input files
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

        // Step 4: Match each spectrum against peptide candidates
        let mut psms: Vec<Psm> = Vec::new();

        for spectrum in &all_spectra {
            if let Some(m) = match_spectrum(
                spectrum,
                &all_peptides,
                &params.precursor_tolerance,
                &params.fragment_tolerance,
                &params.fixed_modifications,
            ) {
                psms.push(build_psm(spectrum, &m, &params.fixed_modifications));
            }
        }

        // Step 5: Aggregate peptide and protein level results
        let peptides = aggregate_peptides(&psms);
        let protein_results = aggregate_proteins(&psms, &proteins);

        // Step 6: Build summary
        let duration = start.elapsed().as_secs_f64();
        let summary = build_summary(&psms, all_spectra.len() as u64, duration);

        // Step 7: Build metadata
        let engine_info = self.engine_info();
        let mut metadata =
            RunMetadata::new(params.clone(), engine_info.clone(), input_files.to_vec());
        metadata.status = RunStatus::Completed;
        metadata.duration_sec = Some(duration);

        Ok(SearchResult {
            run_id: metadata.run_id,
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
        _on_progress: ProgressCallback,
    ) -> Result<SearchResult, CoreError> {
        self.run_search(params, input_files)
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
}

// ---------------------------------------------------------------------------
// Result building helpers
// ---------------------------------------------------------------------------

fn build_psm(
    spectrum: &Spectrum,
    m: &PeptideMatch,
    fixed_mods: &[protein_copilot_core::search_params::Modification],
) -> Psm {
    // Collect modifications that apply to this peptide
    let mut mods = Vec::new();
    for fm in fixed_mods {
        for ch in m.peptide.sequence.chars() {
            if fm.residues.contains(&ch) {
                mods.push(fm.clone());
                break; // one per mod type
            }
        }
    }

    Psm {
        spectrum_scan: spectrum.scan_number,
        peptide_sequence: m.peptide.sequence.clone(),
        modifications: mods,
        charge: m.charge,
        precursor_mz: m.observed_mz,
        calculated_mz: m.theoretical_mz,
        delta_mass_ppm: m.delta_mass_ppm,
        score: m.score,
        q_value: None, // simplified engine doesn't compute FDR
        protein_accessions: vec![m.peptide.protein_accession.clone()],
        is_decoy: false,
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

    // Without FDR, treat all PSMs as passing
    let unique_peptides: std::collections::HashSet<&str> =
        psms.iter().map(|p| p.peptide_sequence.as_str()).collect();
    let unique_proteins: std::collections::HashSet<&str> = psms
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
        total_psms as f64 / total_spectra as f64
    } else {
        0.0
    };

    SearchResultSummary {
        total_spectra_searched: total_spectra,
        total_psms,
        psms_at_1pct_fdr: total_psms, // no FDR filtering
        unique_peptides_at_1pct_fdr: unique_peptides.len() as u64,
        protein_groups_at_1pct_fdr: unique_proteins.len() as u64,
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
        let result = engine.run_search(&params, &[]);
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
        };
        let engine = SimpleSearchEngine::new();
        let result = engine.run_search(&params, &[PathBuf::from("test.mgf")]);
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
        let result = engine.run_search(&params, &[fixture]).unwrap();

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
        let result = engine.run_search(&params, &[fixture]);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("no candidate peptides"));
    }
}
