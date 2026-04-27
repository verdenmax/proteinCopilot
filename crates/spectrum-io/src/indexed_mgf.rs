//! Indexed MGF reader with O(1) scan lookup via [`ScanIndex`].
//!
//! On [`IndexedMgfReader::open`], the file is scanned once to record
//! byte offsets of each `BEGIN IONS` block. Subsequent
//! [`SpectrumReader::read_spectrum`] calls seek directly to the target
//! block, avoiding a full file scan.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;
use crate::index::{IndexSource, ScanIndex};
use crate::mgf::MgfReader;
use crate::reader::SpectrumReader;

/// MGF reader backed by a [`ScanIndex`] for O(1) scan lookup.
///
/// Bulk operations (`read_all`, `read_summary`, `for_each_spectrum`) delegate
/// to the standard streaming [`MgfReader`]; only `read_spectrum` uses the
/// seek-based path.
pub struct IndexedMgfReader {
    index: ScanIndex,
    path: PathBuf,
}

impl IndexedMgfReader {
    /// Opens an MGF file and builds a scan index.
    ///
    /// Two-pass approach:
    /// 1. Finds all `BEGIN IONS` byte positions
    /// 2. For each position, reads header lines to extract `SCANS=N`
    pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
        let index = build_mgf_index(path)?;
        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }

    /// Returns a reference to the underlying scan index.
    pub fn index(&self) -> &ScanIndex {
        &self.index
    }
}

impl SpectrumReader for IndexedMgfReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        MgfReader.read_all(path)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        MgfReader.read_summary(path)
    }

    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        if path != self.path {
            let canonical_self =
                std::fs::canonicalize(&self.path).unwrap_or_else(|_| self.path.clone());
            let canonical_path = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
            if canonical_self != canonical_path {
                tracing::warn!(
                    "IndexedMgfReader opened for {:?} but read_spectrum called with {:?}; using indexed file",
                    self.path, path
                );
            }
        }
        let offset = self
            .index
            .get_offset(scan)
            .ok_or_else(|| SpectrumIoError::ScanNotFound {
                path: self.path.clone(),
                scan,
            })?;
        parse_single_mgf_block(&self.path, offset, scan)
    }

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        MgfReader.for_each_spectrum(path, handler)
    }

    fn list_scan_meta(
        &self,
        _path: &Path,
    ) -> Result<Vec<crate::reader::ScanMetaInfo>, SpectrumIoError> {
        Ok(self
            .index
            .iter_meta()
            .map(|(&scan, meta)| crate::reader::ScanMetaInfo {
                scan_number: scan,
                ms_level: meta.ms_level,
                rt_min: meta.rt_seconds / 60.0,
                isolation_window: meta.isolation_window,
            })
            .collect())
    }
}

// ---------------------------------------------------------------------------
// Index builder
// ---------------------------------------------------------------------------

/// Builds a [`ScanIndex`] for an MGF file.
///
/// Pass 1: scan the file line-by-line recording the byte offset of each
/// `BEGIN IONS` line.
///
/// Pass 2: for each recorded offset, seek there and read header lines
/// to extract the `SCANS=N` value. If no `SCANS` header is present,
/// the block is numbered sequentially.
fn build_mgf_index(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    // --- Pass 1: collect BEGIN IONS byte offsets ---
    let mut reader = crate::util::open_buffered(path)?;
    let mut begin_offsets: Vec<u64> = Vec::new();
    let mut byte_pos: u64 = 0;
    let mut line = String::new();

    loop {
        let line_start = byte_pos;
        line.clear();
        let bytes_read = reader
            .read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;
        if bytes_read == 0 {
            break;
        }
        byte_pos += bytes_read as u64;

        if line.trim() == "BEGIN IONS" {
            begin_offsets.push(line_start);
        }
    }

    // --- Pass 2: read SCANS= from each block header ---
    let file = std::fs::File::open(path).map_err(|e| SpectrumIoError::IoError {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut buf_reader = BufReader::new(file);
    let mut offsets = HashMap::new();
    let mut fallback_scan: u32 = 0;

    for &offset in &begin_offsets {
        fallback_scan += 1;
        buf_reader
            .seek(SeekFrom::Start(offset))
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;

        let scan = read_scan_from_header(&mut buf_reader, path)?;
        let scan_num = scan.unwrap_or(fallback_scan);
        if offsets.contains_key(&scan_num) {
            tracing::warn!(
                "duplicate scan number {} in MGF file {:?}; keeping first occurrence, skipping later",
                scan_num,
                path,
            );
            continue;
        }
        offsets.insert(scan_num, offset);
    }

    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
}

/// Reads lines starting from the current reader position (expected to be
/// at `BEGIN IONS`) until a `SCANS=` header is found or the header section
/// ends (peak data or `END IONS`). Returns `Some(scan)` if found.
fn read_scan_from_header<R: Read + BufRead>(
    reader: &mut R,
    path: &Path,
) -> Result<Option<u32>, SpectrumIoError> {
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed == "BEGIN IONS" || trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if trimmed == "END IONS" {
            break;
        }
        // If the line does not contain '=', it's a peak data line — header section is over.
        if let Some((key, value)) = trimmed.split_once('=') {
            if key.trim().eq_ignore_ascii_case("SCANS") {
                return Ok(value.trim().parse::<u32>().ok());
            }
        } else {
            // Peak data line; no more headers to find.
            break;
        }
    }
    Ok(None)
}

// ---------------------------------------------------------------------------
// Single-block parser
// ---------------------------------------------------------------------------

/// Parses one MGF block (`BEGIN IONS` … `END IONS`) starting at `offset`.
///
/// Handles:
/// - `PEPMASS=<mz> [<intensity>]`
/// - `CHARGE=<n>[+|-]`
/// - `RTINSECONDS=<seconds>`
/// - `SCANS=<scan_number>`
/// - Two-column peak lines (`mz intensity`)
///
/// Returns a [`Spectrum`] with [`MsLevel::MS2`].
fn parse_single_mgf_block(
    path: &Path,
    offset: u64,
    expected_scan: u32,
) -> Result<Spectrum, SpectrumIoError> {
    let file = std::fs::File::open(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            SpectrumIoError::FileNotFound {
                path: path.to_path_buf(),
            }
        } else {
            SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            }
        }
    })?;

    let mut reader = BufReader::new(file);
    reader
        .seek(SeekFrom::Start(offset))
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;

    let mut scan: Option<u32> = None;
    let mut pepmass_mz: Option<f64> = None;
    let mut pepmass_intensity: Option<f64> = None;
    let mut charge: Option<i32> = None;
    let mut rt_min: Option<f64> = None;
    let mut mz_values: Vec<f64> = Vec::new();
    let mut intensity_values: Vec<f64> = Vec::new();
    let mut in_block = false;

    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;
        if n == 0 {
            break; // EOF — treat as implicit END IONS
        }
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "BEGIN IONS" {
            in_block = true;
            continue;
        }

        if trimmed == "END IONS" {
            break;
        }

        if !in_block {
            continue;
        }

        if let Some((key, value)) = trimmed.split_once('=') {
            match key.trim().to_uppercase().as_str() {
                "PEPMASS" => {
                    let (mz, int) = parse_pepmass(value);
                    pepmass_mz = mz;
                    pepmass_intensity = int;
                }
                "CHARGE" => {
                    charge = parse_charge(value);
                }
                "RTINSECONDS" => {
                    rt_min = value.trim().parse::<f64>().ok().map(|v| v / 60.0);
                }
                "SCANS" => {
                    scan = value.trim().parse::<u32>().ok();
                }
                _ => {} // TITLE etc.: skip
            }
        } else {
            // Peak data line
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(mz), Ok(int)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    mz_values.push(mz);
                    intensity_values.push(int);
                }
            }
        }
    }

    let scan_number = scan.unwrap_or(expected_scan);
    let mz = pepmass_mz.ok_or_else(|| SpectrumIoError::ParseError {
        path: path.to_path_buf(),
        line: 0,
        detail: format!("spectrum scan={scan_number} missing PEPMASS"),
    })?;

    let precursors = vec![PrecursorInfo {
        mz,
        charge,
        intensity: pepmass_intensity,
        isolation_window: None,
        source_scan: None,
    }];

    crate::util::sort_peaks_by_mz(&mut mz_values, &mut intensity_values);

    Spectrum::new(
        scan_number,
        MsLevel::MS2,
        rt_min.unwrap_or(0.0),
        precursors,
        mz_values,
        intensity_values,
    )
    .map_err(|e| SpectrumIoError::ValidationError {
        scan: scan_number,
        detail: e.to_string(),
    })
}

// ---------------------------------------------------------------------------
// Header-value parsers (duplicated from mgf.rs to keep module self-contained)
// ---------------------------------------------------------------------------

/// Parses a CHARGE field value like "2+", "3-", or "2".
fn parse_charge(s: &str) -> Option<i32> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    let value = if let Some(num_str) = s.strip_suffix('+') {
        num_str.trim().parse::<i32>().ok()
    } else if let Some(num_str) = s.strip_suffix('-') {
        num_str.trim().parse::<i32>().ok().map(|v| -v)
    } else {
        s.parse::<i32>().ok()
    };
    value.filter(|&v| v != 0)
}

/// Parses a PEPMASS field value like "471.2561" or "471.2561 1500000.0".
fn parse_pepmass(s: &str) -> (Option<f64>, Option<f64>) {
    let parts: Vec<&str> = s.split_whitespace().collect();
    let mz = parts.first().and_then(|v| v.parse::<f64>().ok());
    let intensity = parts.get(1).and_then(|v| v.parse::<f64>().ok());
    (mz, intensity)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mgf")
    }

    #[test]
    fn open_builds_index() {
        let reader = IndexedMgfReader::open(&fixture_path()).unwrap();
        assert_eq!(reader.index().len(), 10);
    }

    #[test]
    fn read_spectrum_by_index() {
        let reader = IndexedMgfReader::open(&fixture_path()).unwrap();
        let spec = reader.read_spectrum(&fixture_path(), 1).unwrap();
        assert_eq!(spec.scan_number, 1);
    }

    #[test]
    fn indexed_vs_standard_all_scans_match() {
        let indexed = IndexedMgfReader::open(&fixture_path()).unwrap();
        for scan in 1..=10 {
            let idx_spec = indexed.read_spectrum(&fixture_path(), scan).unwrap();
            let std_spec = MgfReader.read_spectrum(&fixture_path(), scan).unwrap();
            assert_eq!(idx_spec.scan_number, std_spec.scan_number);
            assert_eq!(idx_spec.mz_array.len(), std_spec.mz_array.len());
            assert!(
                (idx_spec.retention_time_min - std_spec.retention_time_min).abs() < 0.001,
                "RT mismatch for scan {scan}"
            );
        }
    }

    #[test]
    fn duplicate_scan_keeps_first_occurrence() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let mgf_path = dir.path().join("dup.mgf");
        let mut f = std::fs::File::create(&mgf_path).unwrap();
        // Two spectra with same SCANS=5
        write!(
            f,
            "BEGIN IONS\nSCANS=5\nPEPMASS=500.0\n100.0 1000\nEND IONS\n"
        )
        .unwrap();
        write!(
            f,
            "BEGIN IONS\nSCANS=5\nPEPMASS=600.0\n200.0 2000\nEND IONS\n"
        )
        .unwrap();
        drop(f);

        let index = build_mgf_index(&mgf_path).unwrap();
        // Should have only 1 entry for scan 5
        assert_eq!(index.len(), 1, "duplicate should be deduplicated");

        // Read the spectrum — should be the FIRST one (PEPMASS=500)
        let reader = IndexedMgfReader::open(&mgf_path).unwrap();
        let spec = reader.read_spectrum(&mgf_path, 5).unwrap();
        let prec_mz = spec.precursors[0].mz;
        assert!(
            prec_mz > 499.0 && prec_mz < 501.0,
            "should keep first occurrence with PEPMASS=500, got {}",
            prec_mz
        );
    }

    #[test]
    fn read_spectrum_not_found() {
        let reader = IndexedMgfReader::open(&fixture_path()).unwrap();
        let err = reader.read_spectrum(&fixture_path(), 999).unwrap_err();
        assert!(err.to_string().contains("999"));
    }
}
