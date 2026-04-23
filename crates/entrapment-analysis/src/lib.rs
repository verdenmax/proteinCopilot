//! Entrapment analysis - classify trap-database PSM hits by homology to target proteome.
//!
//! Provides L0-L4 discriminability levels for each trap PSM, identifying
//! razor attribution errors, L/I isomers, near-identical homologs, and true trap hits.

pub mod config;
pub mod digest;
pub mod error;
pub mod levenshtein;
pub mod loader;
pub mod mirror_plot;
pub mod output;
pub mod mod_parser;
pub mod provenance;
pub mod coelution;
pub mod multi_provenance;
pub mod multi_report;
pub mod report;
pub mod similarity;
pub mod tagger;
pub mod types;

pub use error::EntrapmentError;
pub use types::{
    ClassifiedPsm, DiscriminabilityLevel, EntrapmentSummary, LevelCounts, PsmGroup,
    SubstitutionType, UnifiedPsm,
};

use std::path::Path;

use config::EntrapmentConfig;
use digest::TargetDigestIndex;
use similarity::classify_single;
use tagger::Tagger;

// ---------------------------------------------------------------------------
// EntrapmentAnalyzer – public API
// ---------------------------------------------------------------------------

/// High-level API that ties together configuration, tagging, digest index,
/// and similarity classification into a single entry point.
pub struct EntrapmentAnalyzer {
    config: EntrapmentConfig,
    tagger: Tagger,
    index: TargetDigestIndex,
}

impl EntrapmentAnalyzer {
    /// Build a new analyser from a config and a target FASTA database.
    ///
    /// Compiles all tagging rules and builds the in-silico tryptic digest
    /// index from the FASTA file.
    pub fn new(config: EntrapmentConfig, fasta_path: &Path) -> Result<Self, EntrapmentError> {
        let tagger = Tagger::new(&config)?;
        let index = TargetDigestIndex::from_fasta(
            fasta_path,
            config.similarity.max_missed_cleavages,
            config.similarity.max_mismatches,
        )?;
        Ok(Self {
            config,
            tagger,
            index,
        })
    }

    /// Classify a single PSM: tag its group and determine similarity level.
    pub fn classify(&self, psm: &UnifiedPsm) -> Result<ClassifiedPsm, EntrapmentError> {
        let group = self.tagger.tag(&psm.protein_ids)?;
        Ok(classify_single(
            psm,
            group,
            &self.index,
            &self.config.similarity,
        ))
    }

    /// Classify a batch of PSMs.
    pub fn classify_all(&self, psms: &[UnifiedPsm]) -> Result<Vec<ClassifiedPsm>, EntrapmentError> {
        psms.iter().map(|psm| self.classify(psm)).collect()
    }

    /// Compute summary statistics from a set of classified PSMs.
    pub fn summary(&self, classified: &[ClassifiedPsm]) -> EntrapmentSummary {
        let mut level_counts = LevelCounts::default();
        let mut target_psms = 0usize;
        let mut trap_psms = 0usize;
        let mut ambiguous_psms = 0usize;

        // Track L0 razor families
        let mut razor_families: std::collections::HashMap<String, (usize, String, String, String)> =
            std::collections::HashMap::new();

        for cp in classified {
            match cp.group {
                PsmGroup::Target => target_psms += 1,
                PsmGroup::Trap => {
                    trap_psms += 1;
                    level_counts.increment(cp.level);

                    if cp.level == DiscriminabilityLevel::L0 {
                        let family = extract_family_name(&cp.psm.protein_ids);
                        let entry = razor_families.entry(family.clone()).or_insert_with(|| {
                            (
                                0,
                                cp.psm.peptide.clone(),
                                cp.psm.protein_ids.clone(),
                                cp.best_target_protein.clone().unwrap_or_default(),
                            )
                        });
                        entry.0 += 1;
                    }
                }
                PsmGroup::Ambiguous => ambiguous_psms += 1,
            }
        }

        let mut top_razor_families: Vec<types::RazorFamily> = razor_families
            .into_iter()
            .map(|(family, (count, pep, trap, target))| types::RazorFamily {
                family,
                count,
                example_peptide: pep,
                example_trap_protein: trap,
                example_target_protein: target,
            })
            .collect();
        top_razor_families.sort_by(|a, b| b.count.cmp(&a.count));
        top_razor_families.truncate(10);

        EntrapmentSummary {
            total_psms: classified.len(),
            target_psms,
            trap_psms,
            ambiguous_psms,
            level_counts,
            top_razor_families,
        }
    }
}

// ---------------------------------------------------------------------------
// Batch provenance tracing
// ---------------------------------------------------------------------------

/// Run fragment ion provenance tracing on classified PSMs.
///
/// For each PSM whose level is in `config.provenance.levels_to_trace` and has
/// a `best_target_peptide`, reads the corresponding MS2 spectrum from the mzML
/// file and calls [`provenance::trace_provenance`].
///
/// # Arguments
/// * `classified` — mutable slice of classified PSMs (provenance field will be set)
/// * `mzml_dir` — directory containing mzML files (matched by `spectrum_file` name)
/// * `config` — entrapment config with provenance settings
///
/// # Returns
/// Number of PSMs successfully traced
pub fn trace_provenance_batch(
    classified: &mut [ClassifiedPsm],
    mzml_dir: &Path,
    config: &EntrapmentConfig,
) -> Result<u32, EntrapmentError> {
    use std::collections::{HashMap, HashSet};

    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    use crate::provenance::trace_provenance;

    let tolerance = MassTolerance {
        value: config.provenance.fragment_tolerance_ppm,
        unit: ToleranceUnit::Ppm,
    };

    let levels_to_trace: HashSet<&str> = config
        .provenance
        .levels_to_trace
        .iter()
        .map(|s| s.as_str())
        .collect();

    let mut traced_count = 0u32;

    // Group PSMs by spectrum_file for efficient file reading.
    // Build a map: spectrum_file → list of indices to trace.
    // Each entry carries either a scan_number or RT+precursor_mz for lookup.
    let mut file_groups: HashMap<String, Vec<usize>> = HashMap::new();

    for (idx, cpsm) in classified.iter().enumerate() {
        // Only trace trap PSMs (target PSMs get L4 from classify_single but
        // should not be provenance-traced).
        if cpsm.group != PsmGroup::Trap {
            continue;
        }
        if !levels_to_trace.contains(cpsm.level.as_str()) {
            continue;
        }
        if cpsm.provenance.is_some() {
            continue; // already traced
        }

        // Must have either scan_number or retention_time + precursor_mz for RT-based lookup.
        let has_scan = cpsm.psm.scan_number.is_some_and(|s| s > 0);
        let has_rt_mz = cpsm.psm.retention_time.is_some() && cpsm.psm.precursor_mz.is_some();
        if !has_scan && !has_rt_mz {
            continue;
        }

        let spectrum_file = match cpsm.psm.spectrum_file.as_deref() {
            Some(f) if !f.is_empty() => f.to_string(),
            _ => continue,
        };
        file_groups.entry(spectrum_file).or_default().push(idx);
    }

    // Process each file group.
    if file_groups.is_empty() {
        // Check if we had eligible PSMs but all lacked scan numbers
        let eligible_count = classified
            .iter()
            .filter(|c| levels_to_trace.contains(c.level.as_str()) && c.provenance.is_none())
            .count();
        if eligible_count > 0 {
            tracing::warn!(
                eligible_psms = eligible_count,
                "no PSMs have scan_number or spectrum_file — provenance tracing requires \
                 scan-level data. DIA-NN parquet results do not include scan numbers; \
                 consider using a search engine that provides scan-level identifiers."
            );
        }
    }

    for (spectrum_file, indices) in &file_groups {
        // Find the mzML file on disk — skip this file group if not found.
        let mzml_path = match find_mzml_file(mzml_dir, spectrum_file) {
            Ok(p) => p,
            Err(_) => {
                tracing::warn!(
                    spectrum_file = %spectrum_file,
                    skipped_psms = indices.len(),
                    "mzML file not found, skipping provenance for this run"
                );
                continue;
            }
        };

        // Create indexed reader for O(1) scan access.
        let reader = match protein_copilot_spectrum_io::create_indexed_reader(&mzml_path) {
            Ok(r) => r,
            Err(e) => {
                tracing::warn!(
                    file = %mzml_path.display(),
                    error = %e,
                    skipped_psms = indices.len(),
                    "failed to create spectrum reader, skipping this run"
                );
                continue;
            }
        };

        for &idx in indices {
            let cpsm = &classified[idx];

            // Resolve scan number: direct or via RT-based lookup.
            let scan_number = if let Some(s) = cpsm.psm.scan_number.filter(|&s| s > 0) {
                s
            } else if let (Some(rt), Some(mz)) =
                (cpsm.psm.retention_time, cpsm.psm.precursor_mz)
            {
                // Derive tolerance from RT.Start / RT.Stop, or use config default.
                let tol = match (cpsm.psm.rt_start, cpsm.psm.rt_stop) {
                    (Some(start), Some(stop)) if stop > start => (stop - start) / 2.0,
                    _ => config.provenance.rt_tolerance_min,
                };
                match reader.find_by_rt(&mzml_path, rt, mz, tol) {
                    Ok(Some((scan, delta))) => {
                        tracing::debug!(
                            peptide = %cpsm.psm.peptide,
                            rt_min = rt,
                            precursor_mz = mz,
                            resolved_scan = scan,
                            rt_delta_min = delta,
                            "resolved scan via RT-based lookup"
                        );
                        scan
                    }
                    Ok(None) => {
                        tracing::debug!(
                            peptide = %cpsm.psm.peptide,
                            rt_min = rt,
                            precursor_mz = mz,
                            "no matching MS2 scan found for RT-based lookup, skipping"
                        );
                        continue;
                    }
                    Err(e) => {
                        tracing::warn!(
                            peptide = %cpsm.psm.peptide,
                            error = %e,
                            "RT-based scan lookup failed, skipping"
                        );
                        continue;
                    }
                }
            } else {
                continue;
            };

            // Read the MS2 spectrum.
            let spectrum = match reader.read_spectrum(&mzml_path, scan_number) {
                Ok(s) => s,
                Err(e) => {
                    tracing::warn!(
                        scan = scan_number,
                        file = %mzml_path.display(),
                        error = %e,
                        "could not read scan, skipping provenance trace"
                    );
                    continue;
                }
            };

            // Skip spectra with too few peaks.
            if (spectrum.mz_array.len() as u32) < config.provenance.min_peaks_for_analysis {
                continue;
            }

            let trap_seq = &cpsm.psm.peptide;
            let target_seq = cpsm.best_target_peptide.as_deref().unwrap_or("");
            let modifications = &cpsm.psm.modifications;

            let mut prov = trace_provenance(
                &spectrum.mz_array,
                &spectrum.intensity_array,
                trap_seq,
                target_seq,
                modifications,
                &tolerance,
                config.provenance.max_fragment_charge,
            );

            // Set chimeric flag based on config threshold.
            prov.is_chimeric = prov.shared_ratio > config.provenance.chimera_threshold;

            classified[idx].provenance = Some(prov);
            traced_count += 1;
        }
    }

    Ok(traced_count)
}

/// Find the mzML file matching a raw/spectrum file name in the given directory.
///
/// Tries `{name}.mzML`, `{name}.mzml`, and the bare name (in case it already
/// has an extension).
fn find_mzml_file(dir: &Path, raw_file: &str) -> Result<std::path::PathBuf, EntrapmentError> {
    for ext in &["mzML", "mzml"] {
        let path = dir.join(format!("{}.{}", raw_file, ext));
        if path.exists() {
            return Ok(path);
        }
    }
    // Try if raw_file already has extension.
    let path = dir.join(raw_file);
    if path.exists() {
        return Ok(path);
    }

    Err(EntrapmentError::SpectrumError {
        path: dir.to_path_buf(),
        detail: format!("mzML file not found for spectrum_file '{}'", raw_file),
    })
}

// ---------------------------------------------------------------------------
// v4 Multi-Target Provenance Batch Pipeline
// ---------------------------------------------------------------------------

/// Extract DIA isolation windows from an mzML file via its MS2 scan metadata.
///
/// De-duplicates windows by rounding center/lower/upper to 2 decimal places.
/// Returns an empty list if the reader fails to list MS2 metadata.
pub fn extract_dia_windows(
    reader: &dyn protein_copilot_spectrum_io::SpectrumReader,
    mzml_path: &std::path::Path,
) -> Vec<coelution::DiaWindow> {
    let mut seen = std::collections::HashSet::new();
    let mut windows = Vec::new();

    if let Ok(meta_list) = reader.list_ms2_meta(mzml_path) {
        for meta in &meta_list {
            if let Some((center, lower, upper)) = meta.isolation_window {
                // Round to 2 decimal places to deduplicate.
                let key = (
                    (center * 100.0) as i64,
                    (lower * 100.0) as i64,
                    (upper * 100.0) as i64,
                );
                if seen.insert(key) {
                    windows.push(coelution::DiaWindow {
                        center,
                        low: center - lower,
                        high: center + upper,
                    });
                }
            }
        }
    }

    windows
}

/// Run multi-target provenance tracing on classified PSMs.
///
/// For each L2/L3 trap PSM:
/// 1. Find co-eluting targets from the full PSM list
/// 2. Read the MS2 spectrum from mzML
/// 3. Match observed peaks against all candidate theoretical ions
/// 4. Generate per-PSM HTML report
///
/// Returns `(traced_count, provenance_results)`.
pub fn trace_multi_target_provenance(
    classified: &[ClassifiedPsm],
    all_psms: &[UnifiedPsm],
    all_groups: &[PsmGroup],
    mzml_dir: &Path,
    config: &EntrapmentConfig,
    output_dir: &Path,
) -> Result<(u32, Vec<types::MultiTargetProvenance>), EntrapmentError> {
    use std::collections::HashSet;

    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    let tolerance = MassTolerance {
        value: config.provenance.fragment_tolerance_ppm,
        unit: ToleranceUnit::Ppm,
    };

    let levels_to_trace: HashSet<&str> = config
        .provenance
        .levels_to_trace
        .iter()
        .map(|s| s.as_str())
        .collect();

    // Collect eligible trap PSMs: group==Trap, level in levels_to_trace,
    // has scan_number or (retention_time + precursor_mz), has spectrum_file.
    struct Eligible<'a> {
        cpsm: &'a ClassifiedPsm,
        run: String,
    }

    let mut eligible: Vec<Eligible<'_>> = Vec::new();
    for cpsm in classified {
        if cpsm.group != PsmGroup::Trap {
            continue;
        }
        if !levels_to_trace.contains(cpsm.level.as_str()) {
            continue;
        }
        let has_scan = cpsm.psm.scan_number.is_some_and(|s| s > 0);
        let has_rt_mz = cpsm.psm.retention_time.is_some() && cpsm.psm.precursor_mz.is_some();
        if !has_scan && !has_rt_mz {
            continue;
        }
        let run = match cpsm.psm.spectrum_file.as_deref() {
            Some(f) if !f.is_empty() => f.to_string(),
            _ => continue,
        };
        eligible.push(Eligible { cpsm, run });
    }

    // Get unique run names (preserving order).
    let mut run_names: Vec<String> = Vec::new();
    {
        let mut seen_runs: HashSet<&str> = HashSet::new();
        for e in &eligible {
            if seen_runs.insert(&e.run) {
                run_names.push(e.run.clone());
            }
        }
    }

    // Create provenance output directory if per-PSM reports requested.
    if config.provenance.generate_per_psm_reports {
        let prov_dir = output_dir.join("provenance");
        std::fs::create_dir_all(&prov_dir).map_err(|e| EntrapmentError::IoError {
            path: prov_dir.clone(),
            detail: format!("failed to create provenance directory: {e}"),
        })?;
    }

    let mut traced_count = 0u32;
    let mut results: Vec<types::MultiTargetProvenance> = Vec::new();

    for run_name in &run_names {
        // Find mzML file.
        let mzml_path = match find_mzml_file(mzml_dir, run_name) {
            Ok(p) => p,
            Err(_) => {
                let skip_count = eligible.iter().filter(|e| &e.run == run_name).count();
                tracing::warn!(
                    spectrum_file = %run_name,
                    skipped_psms = skip_count,
                    "mzML file not found, skipping multi-target provenance for this run"
                );
                continue;
            }
        };

        // Create indexed reader for O(1) scan access.
        let reader = match protein_copilot_spectrum_io::create_indexed_reader(&mzml_path) {
            Ok(r) => r,
            Err(e) => {
                let skip_count = eligible.iter().filter(|e| &e.run == run_name).count();
                tracing::warn!(
                    file = %mzml_path.display(),
                    error = %e,
                    skipped_psms = skip_count,
                    "failed to create spectrum reader, skipping this run"
                );
                continue;
            }
        };

        // Extract DIA windows.
        let dia_windows = extract_dia_windows(reader.as_ref(), &mzml_path);

        // Build co-elution index.
        let silac_ref = config.provenance.silac.as_ref();
        let index = coelution::CoElutionIndex::build(
            all_psms,
            all_groups,
            &dia_windows,
            silac_ref,
            config.provenance.max_co_eluting_candidates,
        );

        // Process each eligible PSM in this run.
        for e in eligible.iter().filter(|e| &e.run == run_name) {
            let cpsm = e.cpsm;

            // Resolve scan number: direct or via RT-based lookup.
            let scan_number = if let Some(s) = cpsm.psm.scan_number.filter(|&s| s > 0) {
                s
            } else if let (Some(rt), Some(mz)) =
                (cpsm.psm.retention_time, cpsm.psm.precursor_mz)
            {
                let tol = match (cpsm.psm.rt_start, cpsm.psm.rt_stop) {
                    (Some(start), Some(stop)) if stop > start => (stop - start) / 2.0,
                    _ => config.provenance.rt_tolerance_min,
                };
                match reader.find_by_rt(&mzml_path, rt, mz, tol) {
                    Ok(Some((scan, _delta))) => scan,
                    Ok(None) => {
                        tracing::debug!(
                            peptide = %cpsm.psm.peptide,
                            rt_min = rt,
                            precursor_mz = mz,
                            "no matching MS2 scan found for RT-based lookup, skipping"
                        );
                        continue;
                    }
                    Err(err) => {
                        tracing::warn!(
                            peptide = %cpsm.psm.peptide,
                            error = %err,
                            "RT-based scan lookup failed, skipping"
                        );
                        continue;
                    }
                }
            } else {
                continue;
            };

            // Read the MS2 spectrum.
            let spectrum = match reader.read_spectrum(&mzml_path, scan_number) {
                Ok(s) => s,
                Err(err) => {
                    tracing::warn!(
                        scan = scan_number,
                        file = %mzml_path.display(),
                        error = %err,
                        "could not read scan, skipping multi-target provenance"
                    );
                    continue;
                }
            };

            // Skip spectra with too few peaks.
            if (spectrum.mz_array.len() as u32) < config.provenance.min_peaks_for_analysis {
                continue;
            }

            // Find co-eluting candidates.
            let candidates = index.find_co_eluting(&cpsm.psm, run_name);
            if candidates.is_empty() {
                continue;
            }

            // --- Light mirror ---
            let mut light_mirror = multi_provenance::trace_multi_target(
                &spectrum.mz_array,
                &spectrum.intensity_array,
                &cpsm.psm.peptide,
                &cpsm.psm.modifications,
                &candidates,
                &tolerance,
                config.provenance.max_fragment_charge,
            );
            light_mirror.scan_number = scan_number;

            // --- Heavy mirror ---
            let heavy_mirror = if let Some(silac) = &config.provenance.silac {
                let seq_chars: Vec<char> = cpsm.psm.peptide.chars().collect();
                let total_delta: f64 = seq_chars
                    .iter()
                    .map(|&c| match c {
                        'K' => silac.heavy_k_delta,
                        'R' => silac.heavy_r_delta,
                        _ => 0.0,
                    })
                    .sum();

                let has_heavy_candidates = candidates.iter().any(|c| !matches!(c.label_form, types::LabelForm::Light));

                if total_delta > 0.0 && has_heavy_candidates {
                    let charge = cpsm.psm.charge.unwrap_or(2) as f64;
                    let heavy_mz = cpsm.psm.precursor_mz.unwrap_or(0.0) + total_delta / charge;

                    let heavy_tol = match (cpsm.psm.rt_start, cpsm.psm.rt_stop) {
                        (Some(start), Some(stop)) if stop > start => (stop - start) / 2.0,
                        _ => config.provenance.rt_tolerance_min,
                    };
                    let heavy_rt = cpsm.psm.retention_time.unwrap_or(0.0);

                    match reader.find_by_rt(&mzml_path, heavy_rt, heavy_mz, heavy_tol) {
                        Ok(Some((heavy_scan, _delta))) => {
                            match reader.read_spectrum(&mzml_path, heavy_scan) {
                                Ok(heavy_spec) => {
                                    // Generate trap heavy theoretical ions.
                                    let residue_deltas: Vec<(usize, f64)> = seq_chars
                                        .iter()
                                        .enumerate()
                                        .filter_map(|(i, &c)| match c {
                                            'K' => Some((i, silac.heavy_k_delta)),
                                            'R' => Some((i, silac.heavy_r_delta)),
                                            _ => None,
                                        })
                                        .collect();

                                    let trap_light_ions =
                                        crate::provenance::generate_theoretical_ions(
                                            &cpsm.psm.peptide,
                                            &cpsm.psm.modifications,
                                            config.provenance.max_fragment_charge,
                                        );
                                    let trap_heavy_ions =
                                        multi_provenance::shift_ions_heavy(
                                            &cpsm.psm.peptide,
                                            &trap_light_ions,
                                            &residue_deltas,
                                            config.provenance.max_fragment_charge,
                                        );

                                    let mut heavy_data =
                                        multi_provenance::trace_mirror_with_trap_ions(
                                            &heavy_spec.mz_array,
                                            &heavy_spec.intensity_array,
                                            &trap_heavy_ions,
                                            &candidates,
                                            &tolerance,
                                            config.provenance.max_fragment_charge,
                                        );
                                    heavy_data.scan_number = heavy_scan;

                                    Some(heavy_data)
                                }
                                Err(err) => {
                                    tracing::debug!(
                                        scan = heavy_scan,
                                        error = %err,
                                        "could not read heavy scan, skipping heavy mirror"
                                    );
                                    None
                                }
                            }
                        }
                        Ok(None) => {
                            tracing::debug!(
                                peptide = %cpsm.psm.peptide,
                                heavy_mz = heavy_mz,
                                "no heavy scan found, skipping heavy mirror"
                            );
                            None
                        }
                        Err(err) => {
                            tracing::debug!(error = %err, "heavy scan lookup failed");
                            None
                        }
                    }
                } else {
                    None
                }
            } else {
                None
            };

            // Compute trap heavy precursor m/z for metadata.
            let trap_precursor_mz_heavy = config.provenance.silac.as_ref().and_then(|silac| {
                let total_delta: f64 = cpsm.psm.peptide.chars()
                    .map(|c| match c {
                        'K' => silac.heavy_k_delta,
                        'R' => silac.heavy_r_delta,
                        _ => 0.0,
                    })
                    .sum();
                if total_delta > 0.0 {
                    let charge = cpsm.psm.charge.unwrap_or(2) as f64;
                    Some(cpsm.psm.precursor_mz.unwrap_or(0.0) + total_delta / charge)
                } else {
                    None
                }
            });

            let prov = types::MultiTargetProvenance {
                trap_peptide: cpsm.psm.peptide.clone(),
                trap_precursor_mz: cpsm.psm.precursor_mz.unwrap_or(0.0),
                trap_precursor_mz_heavy,
                trap_charge: cpsm.psm.charge.unwrap_or(0),
                spectrum_file: run_name.clone(),
                candidates,
                light: light_mirror,
                heavy: heavy_mirror,
            };

            // Generate per-PSM HTML report if configured.
            if config.provenance.generate_per_psm_reports {
                let report_path = output_dir
                    .join("provenance")
                    .join(format!(
                        "{}_{}_scan{}.html",
                        run_name, cpsm.psm.peptide, prov.light.scan_number
                    ));
                if let Err(err) = multi_report::render_multi_provenance_report(&prov, &report_path)
                {
                    tracing::warn!(
                        path = %report_path.display(),
                        error = %err,
                        "failed to write per-PSM provenance report"
                    );
                }
            }

            results.push(prov);
            traced_count += 1;
        }
    }

    // Write summary report.
    let summary_path = output_dir.join("provenance_summary.html");
    if let Err(err) = multi_report::render_provenance_summary(&results, &summary_path) {
        tracing::warn!(
            path = %summary_path.display(),
            error = %err,
            "failed to write provenance summary report"
        );
    }

    Ok((traced_count, results))
}

/// Extract the protein-family name from a UniProt-format accession.
///
/// E.g. `"sp|P12345|EF1A1_HUMAN"` → `"EF1A1"`.
/// Falls back to returning the first semicolon-delimited accession when
/// the format is not recognised.
fn extract_family_name(accession: &str) -> String {
    let first = accession.split(';').next().unwrap_or(accession);
    if let Some(bar_part) = first.split('|').nth(2) {
        if let Some(name) = bar_part.split('_').next() {
            return name.to_string();
        }
    }
    first.to_string()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_family_name_uniprot() {
        assert_eq!(extract_family_name("sp|P12345|EF1A1_HUMAN"), "EF1A1");
    }

    #[test]
    fn test_extract_family_name_multi_accession() {
        assert_eq!(
            extract_family_name("sp|P12345|EF1A1_HUMAN;sp|Q99999|EF1A2_HUMAN"),
            "EF1A1"
        );
    }

    #[test]
    fn test_extract_family_name_non_uniprot() {
        assert_eq!(extract_family_name("SOME_RANDOM_ID"), "SOME_RANDOM_ID");
    }

    // -----------------------------------------------------------------------
    // find_mzml_file tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_find_mzml_file_not_found() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let result = find_mzml_file(dir.path(), "nonexistent");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("not found"),
            "error should mention 'not found', got: {msg}"
        );
    }

    #[test]
    fn test_find_mzml_file_found_mzml() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let path = dir.path().join("test.mzML");
        std::fs::write(&path, "dummy").expect("write file");
        let result = find_mzml_file(dir.path(), "test");
        assert!(result.is_ok());
        assert_eq!(result.expect("should find"), path);
    }

    #[test]
    fn test_find_mzml_file_found_lowercase() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let path = dir.path().join("sample.mzml");
        std::fs::write(&path, "dummy").expect("write file");
        let result = find_mzml_file(dir.path(), "sample");
        assert!(result.is_ok());
        assert_eq!(result.expect("should find"), path);
    }

    #[test]
    fn test_find_mzml_file_bare_name() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let path = dir.path().join("data.raw.mzML");
        std::fs::write(&path, "dummy").expect("write file");
        // The bare name already includes extension
        let result = find_mzml_file(dir.path(), "data.raw.mzML");
        assert!(result.is_ok());
        assert_eq!(result.expect("should find"), path);
    }

    // -----------------------------------------------------------------------
    // trace_provenance_batch tests
    // -----------------------------------------------------------------------

    /// Helper to build a minimal config that doesn't require YAML validation.
    fn test_config() -> EntrapmentConfig {
        use config::*;
        EntrapmentConfig {
            version: 1,
            target: GroupConfig {
                rules: vec![Rule::AccessionContains {
                    any_of: vec!["_HUMAN".to_string()],
                }],
                fasta: None,
                accession_list: None,
            },
            trap: GroupConfig {
                rules: vec![Rule::AccessionContains {
                    any_of: vec!["_YEAST".to_string()],
                }],
                fasta: None,
                accession_list: None,
            },
            conflict_resolution: ConflictResolution::default(),
            unmatched: UnmatchedPolicy::default(),
            similarity: SimilarityConfig::default(),
            provenance: ProvenanceConfig::default(),
        }
    }

    #[test]
    fn test_trace_provenance_batch_empty_slice() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let config = test_config();
        let mut classified: Vec<ClassifiedPsm> = vec![];
        let result = trace_provenance_batch(&mut classified, dir.path(), &config);
        assert!(result.is_ok());
        assert_eq!(result.expect("success"), 0);
    }

    #[test]
    fn test_trace_provenance_batch_no_eligible_levels() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let mut config = test_config();
        // Only trace L4; but our PSM is L0
        config.provenance.levels_to_trace = vec!["L4".to_string()];

        let psm = UnifiedPsm {
            peptide: "PEPTIDE".into(),
            charge: Some(2),
            precursor_mz: Some(400.0),
            retention_time: None,
            rt_start: None,
            rt_stop: None,
            scan_number: Some(1),
            spectrum_file: Some("test".into()),
            protein_ids: "P1".into(),
            q_value: Some(0.01),
            modifications: vec![],
        };
        let mut classified = vec![ClassifiedPsm {
            psm,
            group: PsmGroup::Trap,
            level: DiscriminabilityLevel::L0,
            best_target_peptide: Some("PEPTIDE".into()),
            best_target_protein: Some("P2".into()),
            mismatches: Some(0),
            delta_mass_da: Some(0.0),
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: Some(0),
            alignment_detail: None,
            provenance: None,
        }];

        let result = trace_provenance_batch(&mut classified, dir.path(), &config);
        assert!(result.is_ok());
        assert_eq!(result.expect("success"), 0);
        assert!(classified[0].provenance.is_none());
    }

    #[test]
    fn test_trace_provenance_batch_skip_no_scan_number() {
        let dir = tempfile::tempdir().expect("create tempdir");
        let config = test_config();

        let psm = UnifiedPsm {
            peptide: "PEPTIDE".into(),
            charge: Some(2),
            precursor_mz: Some(400.0),
            retention_time: None,
            rt_start: None,
            rt_stop: None,
            scan_number: None, // no scan number
            spectrum_file: Some("test".into()),
            protein_ids: "P1".into(),
            q_value: Some(0.01),
            modifications: vec![],
        };
        let mut classified = vec![ClassifiedPsm {
            psm,
            group: PsmGroup::Trap,
            level: DiscriminabilityLevel::L3,
            best_target_peptide: Some("PEPTDDE".into()),
            best_target_protein: Some("P2".into()),
            mismatches: Some(1),
            delta_mass_da: Some(0.5),
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: Some(1),
            alignment_detail: None,
            provenance: None,
        }];

        let result = trace_provenance_batch(&mut classified, dir.path(), &config);
        assert!(result.is_ok());
        assert_eq!(result.expect("success"), 0);
        assert!(classified[0].provenance.is_none());
    }

    #[test]
    fn test_trace_provenance_batch_skip_already_traced() {
        use crate::provenance::FragmentProvenance;

        let dir = tempfile::tempdir().expect("create tempdir");
        let config = test_config();

        let psm = UnifiedPsm {
            peptide: "PEPTIDE".into(),
            charge: Some(2),
            precursor_mz: Some(400.0),
            retention_time: None,
            rt_start: None,
            rt_stop: None,
            scan_number: Some(1),
            spectrum_file: Some("test".into()),
            protein_ids: "P1".into(),
            q_value: Some(0.01),
            modifications: vec![],
        };
        let existing_prov = FragmentProvenance {
            trap_sequence: "PEPTIDE".into(),
            target_sequence: "PEPTIDE".into(),
            annotated_peaks: vec![],
            trap_matched_count: 0,
            target_matched_count: 0,
            shared_count: 0,
            unassigned_count: 0,
            shared_ratio: 0.0,
            is_chimeric: false,
        };
        let mut classified = vec![ClassifiedPsm {
            psm,
            group: PsmGroup::Trap,
            level: DiscriminabilityLevel::L3,
            best_target_peptide: Some("PEPTDDE".into()),
            best_target_protein: Some("P2".into()),
            mismatches: Some(1),
            delta_mass_da: Some(0.5),
            diff_positions: None,
            substitution_type: SubstitutionType::None,
            edit_distance: Some(1),
            alignment_detail: None,
            provenance: Some(existing_prov),
        }];

        // Already traced — should not try to read file or fail
        let result = trace_provenance_batch(&mut classified, dir.path(), &config);
        assert!(result.is_ok());
        assert_eq!(result.expect("success"), 0);
    }
}
