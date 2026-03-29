//! MGF (Mascot Generic Format) spectrum reader.
//!
//! MGF is a text-based format where each spectrum is delimited by
//! `BEGIN IONS` / `END IONS` blocks containing header key-value pairs
//! and m/z + intensity peak lines.
//!
//! Supported header fields:
//! - `PEPMASS=<mz> [<intensity>]`
//! - `CHARGE=<n>[+|-]` or `<n>`
//! - `RTINSECONDS=<seconds>`
//! - `SCANS=<scan_number>`
//! - `TITLE=<text>` (parsed but not stored in Spectrum)

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use protein_copilot_core::spectrum::{
    MsLevel, PrecursorInfo, Spectrum, SpectrumFormat, SpectrumSummary,
};

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Reader for MGF (Mascot Generic Format) spectrum files.
pub struct MgfReader;

/// Intermediate state while parsing a single MGF spectrum block.
#[derive(Default)]
struct MgfBlock {
    scan: Option<u32>,
    pepmass_mz: Option<f64>,
    pepmass_intensity: Option<f64>,
    charge: Option<i32>,
    rt_sec: Option<f64>,
    mz_values: Vec<f64>,
    intensity_values: Vec<f64>,
}

impl MgfBlock {
    /// Converts the parsed block into a validated `Spectrum`.
    fn into_spectrum(
        self,
        fallback_scan: u32,
        path: &Path,
        line_start: usize,
    ) -> Result<Spectrum, SpectrumIoError> {
        let scan = self.scan.unwrap_or(fallback_scan);
        let mz = self.pepmass_mz.ok_or_else(|| SpectrumIoError::ParseError {
            path: path.to_path_buf(),
            line: line_start,
            detail: format!("spectrum scan={scan} missing PEPMASS"),
        })?;

        let precursors = vec![PrecursorInfo {
            mz,
            charge: self.charge,
            intensity: self.pepmass_intensity,
            isolation_window: None, // mgf does not carry isolation window
        }];

        // Sort m/z + intensity arrays together by m/z ascending.
        // Real-world MGF files often have unsorted peaks.
        let mut mz_values = self.mz_values;
        let mut intensity_values = self.intensity_values;
        if !mz_values.windows(2).all(|w| w[0] <= w[1]) {
            let mut indices: Vec<usize> = (0..mz_values.len()).collect();
            indices.sort_by(|&a, &b| mz_values[a].partial_cmp(&mz_values[b]).unwrap());
            mz_values = indices.iter().map(|&i| mz_values[i]).collect();
            intensity_values = indices.iter().map(|&i| intensity_values[i]).collect();
        }

        Spectrum::new(
            scan,
            MsLevel::MS2, // MGF spectra are always MS2
            self.rt_sec.unwrap_or(0.0),
            precursors,
            mz_values,
            intensity_values,
        )
        .map_err(|e| SpectrumIoError::ValidationError {
            scan,
            detail: e.to_string(),
        })
    }
}

/// Parses a CHARGE field value like "2+", "3-", or "2".
fn parse_charge(s: &str) -> Option<i32> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    if let Some(num_str) = s.strip_suffix('+') {
        num_str.trim().parse::<i32>().ok()
    } else if let Some(num_str) = s.strip_suffix('-') {
        num_str.trim().parse::<i32>().ok().map(|v| -v)
    } else {
        s.parse::<i32>().ok()
    }
}

/// Parses a PEPMASS field value like "471.2561" or "471.2561 1500000.0".
fn parse_pepmass(s: &str) -> (Option<f64>, Option<f64>) {
    let parts: Vec<&str> = s.split_whitespace().collect();
    let mz = parts.first().and_then(|v| v.parse::<f64>().ok());
    let intensity = parts.get(1).and_then(|v| v.parse::<f64>().ok());
    (mz, intensity)
}

/// Opens a buffered reader for the given path.
fn open_file(path: &Path) -> Result<BufReader<File>, SpectrumIoError> {
    let file = File::open(path).map_err(|e| {
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
    Ok(BufReader::new(file))
}

/// Streaming parser that calls `handler` for each parsed spectrum.
/// Returns the total number of spectra processed.
/// Streaming parser that calls `handler` for each parsed spectrum.
/// Handler returns `true` to continue parsing or `false` to stop early.
/// Returns the total number of spectra processed.
fn parse_mgf_streaming<F>(path: &Path, mut handler: F) -> Result<u32, SpectrumIoError>
where
    F: FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
{
    let reader = open_file(path)?;
    let mut block: Option<MgfBlock> = None;
    let mut block_start_line: usize = 0;
    let mut fallback_scan: u32 = 0;
    let mut count: u32 = 0;

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "BEGIN IONS" {
            block = Some(MgfBlock::default());
            block_start_line = line_idx + 1; // 1-based
            fallback_scan += 1;
            continue;
        }

        if trimmed == "END IONS" {
            if let Some(b) = block.take() {
                let spectrum = b.into_spectrum(fallback_scan, path, block_start_line)?;
                let keep_going = handler(spectrum)?;
                count += 1;
                if !keep_going {
                    return Ok(count);
                }
            }
            continue;
        }

        // Inside a block
        if let Some(ref mut b) = block {
            if let Some((key, value)) = trimmed.split_once('=') {
                match key.trim().to_uppercase().as_str() {
                    "PEPMASS" => {
                        let (mz, int) = parse_pepmass(value);
                        b.pepmass_mz = mz;
                        b.pepmass_intensity = int;
                    }
                    "CHARGE" => {
                        b.charge = parse_charge(value);
                    }
                    "RTINSECONDS" => {
                        b.rt_sec = value.trim().parse::<f64>().ok();
                    }
                    "SCANS" => {
                        b.scan = value.trim().parse::<u32>().ok();
                    }
                    _ => {} // TITLE and other fields: skip
                }
            } else {
                // Peak line: "mz intensity"
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let (Ok(mz), Ok(int)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                        b.mz_values.push(mz);
                        b.intensity_values.push(int);
                    }
                }
            }
        }
    }

    Ok(count)
}

impl SpectrumReader for MgfReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        let mut spectra = Vec::new();
        parse_mgf_streaming(path, |s| {
            spectra.push(s);
            Ok(true)
        })?;
        Ok(spectra)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        let mut total: u64 = 0;
        let mut mz_min = f64::MAX;
        let mut mz_max = f64::MIN;
        let mut rt_min = f64::MAX;
        let mut rt_max = f64::MIN;
        let mut charge_dist: HashMap<i32, u64> = HashMap::new();
        let mut peak_counts: Vec<u32> = Vec::new();

        parse_mgf_streaming(path, |s| {
            total += 1;

            if let Some(first) = s.mz_array.first() {
                if *first < mz_min {
                    mz_min = *first;
                }
            }
            if let Some(last) = s.mz_array.last() {
                if *last > mz_max {
                    mz_max = *last;
                }
            }

            if s.retention_time_sec < rt_min {
                rt_min = s.retention_time_sec;
            }
            if s.retention_time_sec > rt_max {
                rt_max = s.retention_time_sec;
            }

            for p in &s.precursors {
                if let Some(c) = p.charge {
                    *charge_dist.entry(c).or_insert(0) += 1;
                }
            }

            peak_counts.push(s.num_peaks() as u32);
            Ok(true)
        })?;

        // Handle empty file
        if total == 0 {
            mz_min = 0.0;
            mz_max = 0.0;
            rt_min = 0.0;
            rt_max = 0.0;
        }

        // Compute median peak count
        peak_counts.sort_unstable();
        let median_peaks = if peak_counts.is_empty() {
            0
        } else {
            peak_counts[peak_counts.len() / 2]
        };

        let summary = SpectrumSummary {
            file_path: path.to_string_lossy().to_string(),
            format: SpectrumFormat::Mgf,
            total_spectra: total,
            ms1_count: 0, // MGF only contains MS2
            ms2_count: total,
            mz_range: (mz_min, mz_max),
            rt_range_sec: (rt_min, rt_max),
            precursor_charge_distribution: charge_dist,
            median_peaks_per_spectrum: median_peaks,
        };
        summary
            .validate()
            .map_err(|e| SpectrumIoError::ValidationError {
                scan: 0,
                detail: format!("summary: {e}"),
            })?;
        Ok(summary)
    }

    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        let mut found: Option<Spectrum> = None;
        parse_mgf_streaming(path, |s| {
            if s.scan_number == scan {
                found = Some(s);
                Ok(false) // stop early
            } else {
                Ok(true)
            }
        })?;
        found.ok_or_else(|| SpectrumIoError::ScanNotFound {
            path: path.to_path_buf(),
            scan,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("small.mgf")
    }

    // -- parse helpers --------------------------------------------------

    #[test]
    fn parse_charge_formats() {
        assert_eq!(parse_charge("2+"), Some(2));
        assert_eq!(parse_charge("3+"), Some(3));
        assert_eq!(parse_charge("1+"), Some(1));
        assert_eq!(parse_charge("2-"), Some(-2));
        assert_eq!(parse_charge("2"), Some(2));
        assert_eq!(parse_charge(""), None);
        assert_eq!(parse_charge("4+"), Some(4));
    }

    #[test]
    fn parse_pepmass_formats() {
        let (mz, int) = parse_pepmass("471.2561 1500000.0");
        assert!((mz.unwrap() - 471.2561).abs() < 1e-4);
        assert!((int.unwrap() - 1500000.0).abs() < 1.0);

        let (mz, int) = parse_pepmass("523.7832");
        assert!((mz.unwrap() - 523.7832).abs() < 1e-4);
        assert!(int.is_none());
    }

    // -- read_all -------------------------------------------------------

    #[test]
    fn read_all_parses_10_spectra() {
        let reader = MgfReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        assert_eq!(spectra.len(), 10);
    }

    #[test]
    fn read_all_first_spectrum_correct() {
        let reader = MgfReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let s = &spectra[0];

        assert_eq!(s.scan_number, 1);
        assert_eq!(s.ms_level, MsLevel::MS2);
        assert!((s.retention_time_sec - 120.5).abs() < 0.01);
        assert_eq!(s.num_peaks(), 5);
        assert_eq!(s.precursors.len(), 1);
        assert!((s.precursors[0].mz - 471.2561).abs() < 1e-4);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].intensity.unwrap() - 1500000.0).abs() < 1.0);
    }

    #[test]
    fn read_all_spectrum_without_pepmass_intensity() {
        let reader = MgfReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let s = &spectra[1]; // scan 2: no PEPMASS intensity

        assert_eq!(s.scan_number, 2);
        assert!(s.precursors[0].intensity.is_none());
        assert_eq!(s.precursors[0].charge, Some(3));
    }

    #[test]
    fn read_all_various_charge_states() {
        let reader = MgfReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let charges: Vec<Option<i32>> = spectra.iter().map(|s| s.precursors[0].charge).collect();
        // scans 1-10 charges: 2, 3, 2, 2, 3, 2, 3, 1, 4, 2
        assert_eq!(
            charges,
            vec![
                Some(2),
                Some(3),
                Some(2),
                Some(2),
                Some(3),
                Some(2),
                Some(3),
                Some(1),
                Some(4),
                Some(2),
            ]
        );
    }

    #[test]
    fn read_all_peak_counts() {
        let reader = MgfReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let peaks: Vec<usize> = spectra.iter().map(|s| s.num_peaks()).collect();
        assert_eq!(peaks, vec![5, 4, 6, 7, 3, 5, 6, 3, 8, 4]);
    }

    #[test]
    fn read_all_mz_arrays_sorted() {
        let reader = MgfReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        for s in &spectra {
            for w in s.mz_array.windows(2) {
                assert!(w[0] <= w[1], "m/z not sorted in scan {}", s.scan_number);
            }
        }
    }

    // -- read_summary ---------------------------------------------------

    #[test]
    fn read_summary_correct() {
        let reader = MgfReader;
        let summary = reader.read_summary(&fixture_path()).unwrap();

        assert_eq!(summary.total_spectra, 10);
        assert_eq!(summary.ms1_count, 0);
        assert_eq!(summary.ms2_count, 10);
        assert_eq!(summary.format, SpectrumFormat::Mgf);

        // m/z range: smallest first peak to largest last peak
        assert!(summary.mz_range.0 > 0.0);
        assert!(summary.mz_range.1 > summary.mz_range.0);

        // RT range: 120.5 to 240.0
        assert!((summary.rt_range_sec.0 - 120.5).abs() < 0.1);
        assert!((summary.rt_range_sec.1 - 240.0).abs() < 0.1);

        // Charge distribution should have entries
        assert!(!summary.precursor_charge_distribution.is_empty());
        let charge_2_count = summary.precursor_charge_distribution.get(&2).unwrap_or(&0);
        assert_eq!(*charge_2_count, 5); // scans 1,3,4,6,10

        assert!(summary.median_peaks_per_spectrum > 0);
    }

    #[test]
    fn read_summary_validates() {
        let reader = MgfReader;
        let summary = reader.read_summary(&fixture_path()).unwrap();
        assert!(summary.validate().is_ok());
    }

    // -- read_spectrum --------------------------------------------------

    #[test]
    fn read_spectrum_by_scan() {
        let reader = MgfReader;
        let s = reader.read_spectrum(&fixture_path(), 5).unwrap();
        assert_eq!(s.scan_number, 5);
        assert_eq!(s.num_peaks(), 3);
        assert_eq!(s.precursors[0].charge, Some(3));
    }

    #[test]
    fn read_spectrum_not_found() {
        let reader = MgfReader;
        let err = reader.read_spectrum(&fixture_path(), 999).unwrap_err();
        assert!(err.to_string().contains("999"));
        assert!(err.to_string().contains("not found"));
    }

    // -- error cases ----------------------------------------------------

    #[test]
    fn read_all_file_not_found() {
        let reader = MgfReader;
        let err = reader
            .read_all(Path::new("/nonexistent/file.mgf"))
            .unwrap_err();
        assert!(err.to_string().contains("not found"));
    }
}
