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
    let mut file_groups: HashMap<String, Vec<usize>> = HashMap::new();

    for (idx, cpsm) in classified.iter().enumerate() {
        if !levels_to_trace.contains(cpsm.level.as_str()) {
            continue;
        }
        if cpsm.provenance.is_some() {
            continue; // already traced
        }
        // Need a scan number and spectrum file to read the spectrum.
        let scan_number = match cpsm.psm.scan_number {
            Some(s) => s,
            None => continue,
        };
        if scan_number == 0 {
            continue;
        }
        let spectrum_file = match cpsm.psm.spectrum_file.as_deref() {
            Some(f) if !f.is_empty() => f.to_string(),
            _ => continue,
        };
        file_groups.entry(spectrum_file).or_default().push(idx);
    }

    // Process each file group.
    for (spectrum_file, indices) in &file_groups {
        // Find the mzML file on disk.
        let mzml_path = find_mzml_file(mzml_dir, spectrum_file)?;

        // Create indexed reader for O(1) scan access.
        let reader = protein_copilot_spectrum_io::create_indexed_reader(&mzml_path)
            .map_err(|e| EntrapmentError::SpectrumError {
                path: mzml_path.clone(),
                detail: format!("failed to create reader: {}", e),
            })?;

        for &idx in indices {
            let cpsm = &classified[idx];
            let scan_number = match cpsm.psm.scan_number {
                Some(s) => s,
                None => continue,
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
