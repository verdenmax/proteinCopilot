//! Sage search engine adapter.
//!
//! Integrates sage-core as a library for high-performance proteomics search.
//! Sage runs entirely in-process using rayon for parallelism, bridged to
//! tokio via `spawn_blocking`.

pub mod config;
pub mod convert;

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use protein_copilot_core::diagnostics::SearchDiagnostics;
use protein_copilot_core::engine::{EngineInfo, HealthStatus, SearchEngineAdapter};
use protein_copilot_core::error::CoreError;
use protein_copilot_core::progress::{ProgressCallback, SearchProgress};
use protein_copilot_core::run_metadata::{RunMetadata, RunStatus};
use protein_copilot_core::search_params::SearchParams;
use protein_copilot_core::search_result::{
    PeptideResult, ProteinResult, Psm, SearchResult, SearchResultSummary,
};
use protein_copilot_core::spectrum::{MsLevel, Spectrum};

use rayon::prelude::*;
use sage_core::database::IndexedDatabase;
use sage_core::fasta::Fasta;
use sage_core::mass::PROTON;
use sage_core::scoring::{Feature, ScoreType, Scorer};
use sage_core::spectrum::SpectrumProcessor;
use uuid::Uuid;

use self::config::build_sage_parameters;
use self::convert::{mass_tolerance_to_sage, spectrum_to_raw};

/// Sage search engine adapter.
///
/// Wraps sage-core for in-process proteomics search with rayon parallelism.
/// Each `search()` call creates an independent `IndexedDatabase` and `Scorer`.
#[derive(Debug, Default)]
pub struct SageAdapter {
    /// Number of threads for rayon (0 = use rayon default = all cores).
    /// Reserved for future per-search thread pool configuration.
    #[allow(dead_code)]
    thread_count: usize,
}

impl SageAdapter {
    /// Create a new `SageAdapter`.
    ///
    /// `thread_count` of 0 means "use all available cores" (rayon default).
    pub fn new(thread_count: usize) -> Self {
        Self { thread_count }
    }
}

#[async_trait::async_trait]
impl SearchEngineAdapter for SageAdapter {
    fn engine_info(&self) -> EngineInfo {
        EngineInfo {
            name: "Sage".to_string(),
            version: "0.15.0".to_string(),
            supported_features: vec![
                "open_search".to_string(),
                "lfq".to_string(),
                "tmt".to_string(),
                "chimera".to_string(),
            ],
        }
    }

    async fn health_check(&self) -> Result<HealthStatus, CoreError> {
        // Sage runs in-process — always healthy if linked successfully.
        Ok(HealthStatus::Healthy)
    }

    async fn search(
        &self,
        params: &SearchParams,
        input_files: &[PathBuf],
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError> {
        let mut all_spectra: Vec<Spectrum> = Vec::new();
        for path in input_files {
            let file_info = protein_copilot_spectrum_io::detect_format(path).map_err(|e| {
                CoreError::SearchEngineError {
                    engine: "Sage".into(),
                    detail: format!("Failed to detect format for {}: {}", path.display(), e),
                    suggestion: "Check that the input file exists and is a valid mzML/mgf file"
                        .into(),
                }
            })?;
            let reader = protein_copilot_spectrum_io::create_reader(&file_info);
            let spectra = reader
                .read_all(path)
                .map_err(|e| CoreError::SearchEngineError {
                    engine: "Sage".into(),
                    detail: format!("Error reading spectra from {}: {}", path.display(), e),
                    suggestion: "Check input file format".into(),
                })?;
            all_spectra.extend(spectra);
        }

        let mut result = self
            .search_with_spectra(params, all_spectra, on_progress, diagnostics)
            .await?;

        // Populate input_files in metadata (search_with_spectra doesn't have them)
        result.metadata.input_files = input_files.to_vec();

        Ok(result)
    }

    async fn search_with_spectra(
        &self,
        params: &SearchParams,
        spectra: Vec<Spectrum>,
        on_progress: ProgressCallback,
        diagnostics: &mut SearchDiagnostics,
    ) -> Result<SearchResult, CoreError> {
        let _ = &diagnostics; // Will be used in Task 6
        let start = Instant::now();
        let run_id = Uuid::new_v4();

        // ── Phase 1: filter to MS2 ───────────────────────────────────────
        let ms2_spectra: Vec<&Spectrum> = spectra
            .iter()
            .filter(|s| s.ms_level == MsLevel::MS2)
            .collect();

        if ms2_spectra.is_empty() {
            return Err(CoreError::SearchEngineError {
                engine: "Sage".into(),
                detail: "No MS2 spectra found in input".into(),
                suggestion: "Check that input files contain MS2 spectra".into(),
            });
        }

        let total_spectra = ms2_spectra.len();

        // Convert Spectrum → sage RawSpectrum
        let raw_spectra: Vec<sage_core::spectrum::RawSpectrum> = ms2_spectra
            .iter()
            .enumerate()
            .map(|(i, s)| spectrum_to_raw(s, i))
            .collect();

        on_progress(SearchProgress {
            run_id,
            status: "Running".to_string(),
            stage: Some("Spectra loaded".to_string()),
            progress_pct: Some(5.0),
            elapsed_sec: start.elapsed().as_secs_f64(),
            estimated_remaining_sec: None,
            error_category: None,
            has_diagnostics: false,
        });

        // ── Phase 2: Build DB + Score (rayon via spawn_blocking) ─────────
        let sage_params = build_sage_parameters(params);
        let precursor_tol = mass_tolerance_to_sage(&params.precursor_tolerance);
        let fragment_tol = mass_tolerance_to_sage(&params.fragment_tolerance);

        // Read FASTA content
        let fasta_path = &params.database_path;
        let fasta_content = tokio::fs::read_to_string(fasta_path).await.map_err(|e| {
            CoreError::SearchEngineError {
                engine: "Sage".into(),
                detail: format!("Failed to read FASTA file {}: {}", fasta_path, e),
                suggestion: "Check that database_path points to a valid FASTA file".into(),
            }
        })?;

        // Progress counter shared between rayon workers and the tokio poll task.
        let progress_counter = Arc::new(AtomicUsize::new(0));
        let counter_for_poll = Arc::clone(&progress_counter);

        // Wrap the progress callback in Arc so we can share between tasks.
        // ProgressCallback is `Box<dyn Fn(SearchProgress) + Send + Sync>`.
        let on_progress_arc: Arc<dyn Fn(SearchProgress) + Send + Sync> = Arc::from(on_progress);
        let on_progress_for_poll = Arc::clone(&on_progress_arc);

        let progress_run_id = run_id;
        let progress_start = start;
        let total_for_poll = total_spectra;

        let progress_handle = tokio::spawn(async move {
            let mut interval = tokio::time::interval(Duration::from_millis(500));
            loop {
                interval.tick().await;
                let done = counter_for_poll.load(Ordering::Relaxed);
                let pct = (done as f64 / total_for_poll.max(1) as f64) * 100.0;
                on_progress_for_poll(SearchProgress {
                    run_id: progress_run_id,
                    status: "Running".to_string(),
                    stage: Some(format!("Searching ({}/{})", done, total_for_poll)),
                    progress_pct: Some(pct.min(100.0)),
                    elapsed_sec: progress_start.elapsed().as_secs_f64(),
                    estimated_remaining_sec: None,
                    error_category: None,
                    has_diagnostics: false,
                });
                if done >= total_for_poll {
                    break;
                }
            }
        });

        // Main search work in a blocking thread pool.
        let counter_for_search = Arc::clone(&progress_counter);
        let generate_decoys = sage_params.generate_decoys;
        let decoy_tag = sage_params.decoy_tag.clone();
        let search_result = tokio::task::spawn_blocking(move || {
            // Parse FASTA
            let fasta = Fasta::parse(fasta_content, decoy_tag.as_str(), generate_decoys);

            // Build indexed database
            let db = sage_params.build(fasta);

            // Create spectrum processor (top 150 peaks, deisotope, min_deisotope_mz=0)
            let sp = SpectrumProcessor::new(150, true, 0.0);

            // Create scorer
            let scorer = Scorer {
                db: &db,
                precursor_tol,
                fragment_tol,
                min_matched_peaks: 4,
                min_isotope_err: -1,
                max_isotope_err: 3,
                min_precursor_charge: 2,
                max_precursor_charge: 4,
                override_precursor_charge: false,
                max_fragment_charge: None,
                chimera: false,
                report_psms: 1,
                wide_window: false,
                annotate_matches: false,
                score_type: ScoreType::SageHyperScore,
            };

            // Score all spectra in parallel
            let mut features: Vec<Feature> = raw_spectra
                .into_par_iter()
                .flat_map(|raw| {
                    let processed = sp.process(raw);
                    let results = scorer.score(&processed);
                    counter_for_search.fetch_add(1, Ordering::Relaxed);
                    results
                })
                .collect();

            // ── Phase 3: FDR ─────────────────────────────────────────────
            // LDA rescoring
            if sage_core::ml::linear_discriminant::score_psms(&mut features, precursor_tol)
                .is_none()
            {
                // Fallback: heuristic discriminant score (guard against NaN from ln_1p)
                features.par_iter_mut().for_each(|feat| {
                    let p = (-feat.poisson as f32).max(-0.999);
                    feat.discriminant_score = p.ln_1p() + feat.longest_y_pct / 3.0;
                });
            }

            // Sort by discriminant score descending
            features
                .par_sort_unstable_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));

            // Spectrum-level q-value
            sage_core::ml::qvalue::spectrum_q_value(&mut features);

            // Picked peptide FDR
            sage_core::fdr::picked_peptide(&db, &mut features);

            // Picked protein FDR
            sage_core::fdr::picked_protein(&db, &mut features);

            Ok::<(Vec<Feature>, IndexedDatabase), CoreError>((features, db))
        })
        .await
        .map_err(|e| CoreError::SearchEngineError {
            engine: "Sage".into(),
            detail: format!("Search task panicked: {}", e),
            suggestion: "This is likely a bug in the Sage adapter".into(),
        })??;

        // Stop the progress polling task.
        progress_handle.abort();

        let (features, db) = search_result;
        let duration = start.elapsed().as_secs_f64();

        on_progress_arc(SearchProgress {
            run_id,
            status: "Running".to_string(),
            stage: Some("Converting results".to_string()),
            progress_pct: Some(95.0),
            elapsed_sec: duration,
            estimated_remaining_sec: None,
            error_category: None,
            has_diagnostics: false,
        });

        // ── Phase 4: Convert Feature → Psm ──────────────────────────────
        let psms: Vec<Psm> = features
            .iter()
            .map(|feat| feature_to_psm(feat, &db))
            .collect();

        // ── Build peptide-level results ──────────────────────────────────
        let mut peptide_map: HashMap<String, PeptideResult> = HashMap::new();
        for psm in &psms {
            if psm.is_decoy {
                continue;
            }
            let entry = peptide_map
                .entry(psm.peptide_sequence.clone())
                .or_insert_with(|| PeptideResult {
                    sequence: psm.peptide_sequence.clone(),
                    protein_accessions: psm.protein_accessions.clone(),
                    best_score: psm.score,
                    q_value: psm.q_value,
                    psm_count: 0,
                });
            entry.psm_count += 1;
            if psm.score > entry.best_score {
                entry.best_score = psm.score;
                entry.q_value = psm.q_value;
            }
        }
        let peptides: Vec<PeptideResult> = peptide_map.into_values().collect();

        // ── Build protein-level results ──────────────────────────────────
        let mut protein_map: HashMap<String, ProteinResult> = HashMap::new();
        // Track distinct peptides per protein.
        let mut protein_peptides: HashMap<String, HashSet<String>> = HashMap::new();
        // Track which proteins each peptide maps to (for unique_peptide_count).
        let mut peptide_proteins: HashMap<String, HashSet<String>> = HashMap::new();
        for psm in &psms {
            if psm.is_decoy {
                continue;
            }
            for acc in &psm.protein_accessions {
                protein_map
                    .entry(acc.clone())
                    .or_insert_with(|| ProteinResult {
                        accession: acc.clone(),
                        description: String::new(),
                        coverage: 0.0,
                        peptide_count: 0,
                        unique_peptide_count: 0,
                    });
                protein_peptides
                    .entry(acc.clone())
                    .or_default()
                    .insert(psm.peptide_sequence.clone());
                peptide_proteins
                    .entry(psm.peptide_sequence.clone())
                    .or_default()
                    .insert(acc.clone());
            }
        }
        for (acc, pep_set) in &protein_peptides {
            if let Some(entry) = protein_map.get_mut(acc) {
                entry.peptide_count = pep_set.len() as u64;
                // Unique peptides = those mapped to only this protein
                entry.unique_peptide_count = pep_set
                    .iter()
                    .filter(|pep| {
                        peptide_proteins
                            .get(*pep)
                            .is_none_or(|prots| prots.len() == 1)
                    })
                    .count() as u64;
            }
        }
        let proteins: Vec<ProteinResult> = protein_map.into_values().collect();

        // ── Build summary ────────────────────────────────────────────────
        let target_psms: Vec<&Psm> = psms.iter().filter(|p| !p.is_decoy).collect();
        let psms_at_1pct: u64 = target_psms
            .iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .count() as u64;
        let peptides_at_1pct: u64 = peptides
            .iter()
            .filter(|p| p.q_value.is_some_and(|q| q <= 0.01))
            .count() as u64;

        // Median score
        let median_score = {
            let mut scores: Vec<f64> = psms.iter().map(|p| p.score).collect();
            scores.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            if scores.is_empty() {
                0.0
            } else {
                scores[scores.len() / 2]
            }
        };

        // Median delta mass ppm
        let median_delta_mass_ppm = {
            let mut deltas: Vec<f64> = psms.iter().map(|p| p.delta_mass_ppm).collect();
            deltas.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            if deltas.is_empty() {
                0.0
            } else {
                deltas[deltas.len() / 2]
            }
        };

        let identification_rate = if total_spectra > 0 {
            psms_at_1pct as f64 / total_spectra as f64
        } else {
            0.0
        };

        // Modification distribution
        let mut modification_distribution: HashMap<String, u64> = HashMap::new();
        for psm in &psms {
            for m in &psm.modifications {
                *modification_distribution.entry(m.name.clone()).or_default() += 1;
            }
        }

        // Charge distribution
        let mut charge_distribution: HashMap<i32, u64> = HashMap::new();
        for psm in &psms {
            *charge_distribution.entry(psm.charge).or_default() += 1;
        }

        // Count proteins at 1% FDR using protein_q from sage features
        let proteins_at_1pct: u64 = {
            let mut protein_accs_at_1pct: HashSet<&str> = HashSet::new();
            for feat in &features {
                if feat.label != -1 && feat.protein_q <= 0.01 {
                    for acc in &db[feat.peptide_idx].proteins {
                        protein_accs_at_1pct.insert(acc.as_ref());
                    }
                }
            }
            protein_accs_at_1pct.len() as u64
        };

        let summary = SearchResultSummary {
            total_spectra_searched: total_spectra as u64,
            total_psms: psms.len() as u64,
            psms_at_1pct_fdr: psms_at_1pct,
            unique_peptides_at_1pct_fdr: peptides_at_1pct,
            protein_groups_at_1pct_fdr: proteins_at_1pct,
            median_score,
            median_delta_mass_ppm,
            identification_rate,
            modification_distribution,
            charge_distribution,
            search_duration_sec: duration,
        };

        // ── Build metadata ───────────────────────────────────────────────
        let engine_info = self.engine_info();
        let metadata = RunMetadata {
            run_id,
            engine_info: engine_info.clone(),
            params_used: params.clone(),
            input_files: vec![],
            status: RunStatus::Completed,
            created_at: chrono::Utc::now(),
            duration_sec: Some(duration),
        };

        on_progress_arc(SearchProgress {
            run_id,
            status: "Completed".to_string(),
            stage: Some("Done".to_string()),
            progress_pct: Some(100.0),
            elapsed_sec: duration,
            estimated_remaining_sec: Some(0.0),
            error_category: None,
            has_diagnostics: false,
        });

        Ok(SearchResult {
            run_id,
            engine_info,
            params_used: params.clone(),
            psms,
            peptides,
            proteins,
            summary,
            metadata,
        })
    }

    async fn cancel(&self, _run_id: Uuid) -> Result<(), CoreError> {
        // Sage runs in-process; cancellation is handled by the caller
        // aborting the tokio task.
        Ok(())
    }
}

/// Convert a sage `Feature` + `IndexedDatabase` lookup into our `Psm`.
fn feature_to_psm(feat: &Feature, db: &IndexedDatabase) -> Psm {
    let peptide = &db[feat.peptide_idx];

    let peptide_sequence = String::from_utf8_lossy(&peptide.sequence).to_string();
    let protein_accessions: Vec<String> = peptide.proteins.iter().map(|p| p.to_string()).collect();

    // Convert mass → m/z:  mz = (mass + charge * PROTON) / charge
    // PROTON is f32 in sage-core.
    let charge = f64::from(feat.charge);
    let proton = f64::from(PROTON);
    let precursor_mz = (f64::from(feat.expmass) + charge * proton) / charge;
    let calculated_mz = (f64::from(feat.calcmass) + charge * proton) / charge;

    let delta_mass_ppm = if calculated_mz > 0.0 {
        (precursor_mz - calculated_mz) / calculated_mz * 1e6
    } else {
        0.0
    };

    // spec_id is the scan number stored as a string (see convert::spectrum_to_raw)
    let spectrum_scan = feat.spec_id.parse::<u32>().unwrap_or_else(|_| {
        tracing::warn!(spec_id = %feat.spec_id, "Failed to parse spec_id as scan number, defaulting to 1");
        1
    });

    // Preserve sage-specific scoring details in extra fields.
    let mut extra = HashMap::new();
    extra.insert(
        "matched_peaks".into(),
        serde_json::json!(feat.matched_peaks),
    );
    extra.insert("longest_b".into(), serde_json::json!(feat.longest_b));
    extra.insert("longest_y".into(), serde_json::json!(feat.longest_y));
    extra.insert("delta_next".into(), serde_json::json!(feat.delta_next));
    extra.insert("delta_best".into(), serde_json::json!(feat.delta_best));
    extra.insert(
        "discriminant_score".into(),
        serde_json::json!(feat.discriminant_score),
    );
    extra.insert(
        "posterior_error".into(),
        serde_json::json!(feat.posterior_error),
    );
    extra.insert("hyperscore".into(), serde_json::json!(feat.hyperscore));
    extra.insert("peptide_q".into(), serde_json::json!(feat.peptide_q));
    extra.insert("protein_q".into(), serde_json::json!(feat.protein_q));

    let modifications = convert_sage_modifications(peptide);

    Psm {
        spectrum_scan,
        peptide_sequence,
        modifications,
        charge: feat.charge as i32,
        precursor_mz,
        calculated_mz,
        delta_mass_ppm,
        score: feat.hyperscore,
        q_value: Some(f64::from(feat.spectrum_q)),
        protein_accessions,
        is_decoy: feat.label == -1,
        extra: Some(extra),
    }
}

/// Convert sage `Peptide` modifications to our `Modification` structs.
fn convert_sage_modifications(
    peptide: &sage_core::peptide::Peptide,
) -> Vec<protein_copilot_core::search_params::Modification> {
    use protein_copilot_core::search_params::{ModPosition, Modification};

    let mut mods = Vec::new();

    // N-terminal modification
    if let Some(nterm_mass) = peptide.nterm {
        if nterm_mass.abs() > 1e-6 {
            mods.push(Modification {
                name: format!("NTerm({:.4})", nterm_mass),
                mass_delta: f64::from(nterm_mass),
                residues: vec![],
                position: ModPosition::AnyNTerm,
            });
        }
    }

    // Residue modifications
    for (i, &mass) in peptide.modifications.iter().enumerate() {
        if mass.abs() > 1e-6 {
            let residue = if i < peptide.sequence.len() {
                peptide.sequence[i] as char
            } else {
                '?'
            };
            mods.push(Modification {
                name: format!("{}{:.4}", residue, mass),
                mass_delta: f64::from(mass),
                residues: vec![residue],
                position: ModPosition::Anywhere,
            });
        }
    }

    // C-terminal modification
    if let Some(cterm_mass) = peptide.cterm {
        if cterm_mass.abs() > 1e-6 {
            mods.push(Modification {
                name: format!("CTerm({:.4})", cterm_mass),
                mass_delta: f64::from(cterm_mass),
                residues: vec![],
                position: ModPosition::AnyCTerm,
            });
        }
    }

    mods
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn engine_info() {
        let adapter = SageAdapter::default();
        let info = adapter.engine_info();
        assert_eq!(info.name, "Sage");
        assert!(!info.version.is_empty());
    }

    #[tokio::test]
    async fn health_check_always_healthy() {
        let adapter = SageAdapter::default();
        let status = adapter.health_check().await.unwrap();
        assert_eq!(status, HealthStatus::Healthy);
    }

    #[tokio::test]
    async fn search_with_empty_spectra_returns_error() {
        let adapter = SageAdapter::default();
        let params = protein_copilot_core::search_params::SearchParams {
            enzyme: protein_copilot_core::search_params::Enzyme::Trypsin,
            missed_cleavages: 2,
            fixed_modifications: vec![],
            variable_modifications: vec![],
            precursor_tolerance: protein_copilot_core::search_params::MassTolerance {
                value: 20.0,
                unit: protein_copilot_core::search_params::ToleranceUnit::Ppm,
            },
            fragment_tolerance: protein_copilot_core::search_params::MassTolerance {
                value: 0.5,
                unit: protein_copilot_core::search_params::ToleranceUnit::Da,
            },
            database_path: "/nonexistent.fasta".into(),
            decoy_strategy: protein_copilot_core::search_params::DecoyStrategy::Reverse,
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
            engine: Some("Sage".into()),
        };
        let on_progress: ProgressCallback = Box::new(|_| {});
        let result = adapter
            .search_with_spectra(&params, vec![], on_progress, &mut protein_copilot_core::diagnostics::SearchDiagnostics::new())
            .await;
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("No MS2 spectra"), "Error was: {}", err);
    }
}
