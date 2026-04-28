//! FASTA database file parsing.
//!
//! Reads protein sequences from standard FASTA format files.
//! Each entry has a header line starting with `>` followed by
//! one or more sequence lines.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::SearchEngineError;

/// A single protein entry from a FASTA database.
#[derive(Debug, Clone)]
pub struct FastaEntry {
    /// Protein accession (first word after `>`).
    pub accession: String,
    /// Full description line (everything after `>`).
    pub description: String,
    /// Amino acid sequence (uppercase, no whitespace).
    pub sequence: String,
}

/// Parses a FASTA file and returns all protein entries.
pub fn parse_fasta(path: &Path) -> Result<Vec<FastaEntry>, SearchEngineError> {
    let _span = tracing::info_span!("parse_fasta", path = %path.display()).entered();

    let file = File::open(path).map_err(|e| SearchEngineError::FastaError {
        path: path.to_path_buf(),
        detail: format!("cannot open file: {e}"),
    })?;
    let reader = BufReader::new(file);

    let mut entries = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();
    let fasta_progress_interval = 5000;
    let fasta_loop_start = std::time::Instant::now();

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| SearchEngineError::FastaError {
            path: path.to_path_buf(),
            detail: format!("read error: {e}"),
        })?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if let Some(header) = trimmed.strip_prefix('>') {
            // Save previous entry
            if let Some(ref h) = current_header {
                if !current_seq.is_empty() {
                    entries.push(build_entry(h, &current_seq));

                    let count = entries.len();
                    if count % fasta_progress_interval == 0 {
                        let elapsed = fasta_loop_start.elapsed().as_secs_f64();
                        let rate = if elapsed > 0.0 { count as f64 / elapsed } else { 0.0 };
                        tracing::info!(
                            progress = count,
                            rate = format!("{:.0}/s", rate),
                            "parsing proteins"
                        );
                    }
                }
            }
            current_header = Some(header.to_string());
            current_seq.clear();
        } else {
            // Sequence line
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    // Save last entry
    if let Some(ref h) = current_header {
        if !current_seq.is_empty() {
            entries.push(build_entry(h, &current_seq));
        }
    }

    if entries.is_empty() {
        return Err(SearchEngineError::FastaError {
            path: path.to_path_buf(),
            detail: "no protein entries found".to_string(),
        });
    }

    tracing::info!(proteins = entries.len(), "FASTA parsed");

    Ok(entries)
}

fn build_entry(header: &str, sequence: &str) -> FastaEntry {
    let accession = header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string();
    FastaEntry {
        accession,
        description: header.to_string(),
        sequence: sequence.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_fasta(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        write!(f, "{content}").unwrap();
        f
    }

    #[test]
    fn parse_simple_fasta() {
        let f = write_temp_fasta(
            ">sp|P12345|ALBU_HUMAN Serum albumin\n\
             MKWVTFISLLFLFSSAYS\n\
             RGVFRRDAHKSEVAHRFK\n\
             >sp|Q67890|TEST_HUMAN Test protein\n\
             PEPTIDEK\n",
        );
        let entries = parse_fasta(f.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].accession, "sp|P12345|ALBU_HUMAN");
        assert_eq!(entries[0].sequence, "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFK");
        assert_eq!(entries[1].accession, "sp|Q67890|TEST_HUMAN");
        assert_eq!(entries[1].sequence, "PEPTIDEK");
    }

    #[test]
    fn parse_empty_fasta_errors() {
        let f = write_temp_fasta("");
        assert!(parse_fasta(f.path()).is_err());
    }

    #[test]
    fn parse_fasta_file_not_found() {
        let err = parse_fasta(Path::new("/nonexistent/db.fasta")).unwrap_err();
        assert!(err.to_string().contains("cannot open"));
    }

    #[test]
    fn parse_fasta_handles_blank_lines() {
        let f = write_temp_fasta(
            ">protein1\n\
             ACDE\n\
             \n\
             FGHIK\n\
             >protein2\n\
             LMNPQR\n",
        );
        let entries = parse_fasta(f.path()).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sequence, "ACDEFGHIK");
    }
}
