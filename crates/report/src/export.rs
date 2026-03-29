//! Result export — TSV, JSON, and metadata output.

use std::fs;
use std::io::Write;
use std::path::Path;

use protein_copilot_core::run_metadata::RunMetadata;
use protein_copilot_core::search_result::SearchResult;

use crate::error::ReportError;

/// Sanitizes a string for TSV output by replacing tabs and newlines.
fn sanitize_tsv(s: &str) -> String {
    s.replace(['\t', '\n'], " ").replace('\r', "")
}

/// Exports search results as 3 TSV files in the given directory.
///
/// Creates: `psm.tsv`, `peptide.tsv`, `protein.tsv`
pub(crate) fn export_tsv(result: &SearchResult, output_dir: &Path) -> Result<(), ReportError> {
    fs::create_dir_all(output_dir).map_err(|e| ReportError::IoError {
        path: output_dir.to_path_buf(),
        detail: e.to_string(),
    })?;

    // PSM TSV
    let psm_path = output_dir.join("psm.tsv");
    let mut psm_file = create_file(&psm_path)?;
    writeln!(
        psm_file,
        "scan\tsequence\tcharge\tprecursor_mz\tcalculated_mz\tdelta_ppm\tscore\tq_value\tproteins\tis_decoy\tmodifications"
    )
    .map_err(|e| io_err(&psm_path, e))?;

    for psm in &result.psms {
        let mods: Vec<String> = psm
            .modifications
            .iter()
            .map(|m| sanitize_tsv(&m.name))
            .collect();
        writeln!(
            psm_file,
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.2}\t{:.6}\t{}\t{}\t{}\t{}",
            psm.spectrum_scan,
            sanitize_tsv(&psm.peptide_sequence),
            psm.charge,
            psm.precursor_mz,
            psm.calculated_mz,
            psm.delta_mass_ppm,
            psm.score,
            psm.q_value.map_or("NA".to_string(), |q| format!("{q:.6}")),
            psm.protein_accessions
                .iter()
                .map(|a| sanitize_tsv(a))
                .collect::<Vec<_>>()
                .join(";"),
            psm.is_decoy,
            mods.join(";"),
        )
        .map_err(|e| io_err(&psm_path, e))?;
    }

    // Peptide TSV
    let pep_path = output_dir.join("peptide.tsv");
    let mut pep_file = create_file(&pep_path)?;
    writeln!(
        pep_file,
        "sequence\tproteins\tbest_score\tq_value\tpsm_count"
    )
    .map_err(|e| io_err(&pep_path, e))?;

    for pep in &result.peptides {
        writeln!(
            pep_file,
            "{}\t{}\t{:.6}\t{}\t{}",
            sanitize_tsv(&pep.sequence),
            pep.protein_accessions
                .iter()
                .map(|a| sanitize_tsv(a))
                .collect::<Vec<_>>()
                .join(";"),
            pep.best_score,
            pep.q_value.map_or("NA".to_string(), |q| format!("{q:.6}")),
            pep.psm_count,
        )
        .map_err(|e| io_err(&pep_path, e))?;
    }

    // Protein TSV
    let prot_path = output_dir.join("protein.tsv");
    let mut prot_file = create_file(&prot_path)?;
    writeln!(
        prot_file,
        "accession\tdescription\tcoverage\tpeptide_count\tunique_peptide_count"
    )
    .map_err(|e| io_err(&prot_path, e))?;

    for prot in &result.proteins {
        writeln!(
            prot_file,
            "{}\t{}\t{:.4}\t{}\t{}",
            sanitize_tsv(&prot.accession),
            sanitize_tsv(&prot.description),
            prot.coverage,
            prot.peptide_count,
            prot.unique_peptide_count,
        )
        .map_err(|e| io_err(&prot_path, e))?;
    }

    Ok(())
}

/// Exports the complete search result as JSON.
pub(crate) fn export_json(result: &SearchResult, output_path: &Path) -> Result<(), ReportError> {
    let json = serde_json::to_string_pretty(result)
        .map_err(|e| ReportError::SerializationError(format!("SearchResult JSON: {e}")))?;
    fs::write(output_path, json).map_err(|e| ReportError::IoError {
        path: output_path.to_path_buf(),
        detail: e.to_string(),
    })
}

/// Exports run metadata as JSON.
pub(crate) fn export_metadata(
    metadata: &RunMetadata,
    output_path: &Path,
) -> Result<(), ReportError> {
    let json = serde_json::to_string_pretty(metadata)
        .map_err(|e| ReportError::SerializationError(format!("RunMetadata JSON: {e}")))?;
    fs::write(output_path, json).map_err(|e| ReportError::IoError {
        path: output_path.to_path_buf(),
        detail: e.to_string(),
    })
}

fn create_file(path: &Path) -> Result<fs::File, ReportError> {
    fs::File::create(path).map_err(|e| ReportError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    })
}

fn io_err(path: &Path, e: std::io::Error) -> ReportError {
    ReportError::IoError {
        path: path.to_path_buf(),
        detail: e.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::engine::EngineInfo;
    use protein_copilot_core::run_metadata::RunMetadata;
    use protein_copilot_core::search_params::*;
    use protein_copilot_core::search_result::*;
    use std::collections::HashMap;
    use std::path::PathBuf;

    fn sample_result() -> SearchResult {
        let params = SearchParams {
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
        };
        let engine_info = EngineInfo {
            name: "Test".to_string(),
            version: "1.0".to_string(),
            supported_features: vec![],
        };
        let meta = RunMetadata::new(
            params.clone(),
            engine_info.clone(),
            vec![PathBuf::from("/data/test.mgf")],
        );
        SearchResult {
            run_id: meta.run_id,
            engine_info,
            params_used: params,
            psms: vec![Psm {
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
            }],
            peptides: vec![PeptideResult {
                sequence: "PEPTIDER".to_string(),
                protein_accessions: vec!["P001".to_string()],
                best_score: 0.8,
                q_value: Some(0.001),
                psm_count: 1,
            }],
            proteins: vec![ProteinResult {
                accession: "P001".to_string(),
                description: "Test protein".to_string(),
                coverage: 0.4,
                peptide_count: 1,
                unique_peptide_count: 1,
            }],
            summary: SearchResultSummary {
                total_spectra_searched: 10,
                total_psms: 1,
                psms_at_1pct_fdr: 1,
                unique_peptides_at_1pct_fdr: 1,
                protein_groups_at_1pct_fdr: 1,
                median_score: 0.8,
                median_delta_mass_ppm: 0.5,
                identification_rate: 0.1,
                modification_distribution: HashMap::new(),
                charge_distribution: HashMap::new(),
                search_duration_sec: 1.0,
            },
            metadata: meta,
        }
    }

    #[test]
    fn export_tsv_creates_three_files() {
        let dir = tempfile::tempdir().unwrap();
        let result = sample_result();
        export_tsv(&result, dir.path()).unwrap();

        assert!(dir.path().join("psm.tsv").exists());
        assert!(dir.path().join("peptide.tsv").exists());
        assert!(dir.path().join("protein.tsv").exists());
    }

    #[test]
    fn psm_tsv_has_header_and_data() {
        let dir = tempfile::tempdir().unwrap();
        let result = sample_result();
        export_tsv(&result, dir.path()).unwrap();

        let content = fs::read_to_string(dir.path().join("psm.tsv")).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2); // header + 1 PSM
        assert!(lines[0].contains("scan"));
        assert!(lines[0].contains("sequence"));
        assert!(lines[1].contains("PEPTIDER"));
    }

    #[test]
    fn export_json_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let result = sample_result();
        let json_path = dir.path().join("result.json");
        export_json(&result, &json_path).unwrap();

        let content = fs::read_to_string(&json_path).unwrap();
        let back: SearchResult = serde_json::from_str(&content).unwrap();
        assert_eq!(result.run_id, back.run_id);
        assert_eq!(result.psms.len(), back.psms.len());
    }

    #[test]
    fn export_metadata_has_run_id() {
        let dir = tempfile::tempdir().unwrap();
        let result = sample_result();
        let meta_path = dir.path().join("run_metadata.json");
        export_metadata(&result.metadata, &meta_path).unwrap();

        let content = fs::read_to_string(&meta_path).unwrap();
        let back: RunMetadata = serde_json::from_str(&content).unwrap();
        assert_eq!(result.metadata.run_id, back.run_id);
    }
}
