//! Mass spectrometry spectrum data structures.
//!
//! This module defines the core types for representing mass spectrometry data:
//! - [`MsLevel`] — MS acquisition level (MS1, MS2, etc.)
//! - [`PrecursorInfo`] — Precursor ion information
//! - [`Spectrum`] — A single mass spectrum with peak data
//! - [`SpectrumSummary`] — Statistical summary of a spectrum file
//! - [`SpectrumFormat`] — Supported spectrum file formats
//! - [`SpectrumFileInfo`] — Metadata about a spectrum file on disk

use std::collections::HashMap;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use thiserror::Error;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors related to spectrum data validation.
#[derive(Debug, Error)]
pub enum SpectrumError {
    /// mz_array and intensity_array have different lengths.
    #[error("mz_array length ({mz_len}) does not match intensity_array length ({intensity_len})")]
    ArrayLengthMismatch {
        /// Length of mz_array.
        mz_len: usize,
        /// Length of intensity_array.
        intensity_len: usize,
    },
}

// ---------------------------------------------------------------------------
// MsLevel
// ---------------------------------------------------------------------------

/// MS acquisition level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
pub enum MsLevel {
    /// MS1 — survey / full scan.
    MS1,
    /// MS2 — tandem MS / fragmentation scan.
    MS2,
    /// Higher-order MSn or vendor-specific level.
    Other(u8),
}

impl std::fmt::Display for MsLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MsLevel::MS1 => write!(f, "MS1"),
            MsLevel::MS2 => write!(f, "MS2"),
            MsLevel::Other(n) => write!(f, "MS{n}"),
        }
    }
}

// ---------------------------------------------------------------------------
// PrecursorInfo
// ---------------------------------------------------------------------------

/// Precursor ion information for MS2+ spectra.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct PrecursorInfo {
    /// Precursor m/z value (Da).
    pub mz: f64,
    /// Charge state (`None` if undetermined).
    pub charge: Option<i32>,
    /// Intensity (`None` if not available).
    pub intensity: Option<f64>,
}

// ---------------------------------------------------------------------------
// Spectrum
// ---------------------------------------------------------------------------

/// A single mass spectrum with peak data.
///
/// `mz_array` and `intensity_array` must always have the same length.
/// `mz_array` is expected to be sorted in ascending order.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct Spectrum {
    /// Scan number (1-based, following mass spectrometry convention).
    pub scan_number: u32,
    /// MS level of this spectrum.
    pub ms_level: MsLevel,
    /// Retention time in seconds.
    pub retention_time_sec: f64,
    /// Precursor information (present for MS2+ spectra).
    pub precursor: Option<PrecursorInfo>,
    /// Array of m/z values (Da), sorted ascending.
    pub mz_array: Vec<f64>,
    /// Array of intensity values, same length as `mz_array`.
    pub intensity_array: Vec<f64>,
}

impl Spectrum {
    /// Creates a new `Spectrum` with validation.
    ///
    /// Returns `SpectrumError::ArrayLengthMismatch` if `mz_array` and
    /// `intensity_array` have different lengths.
    pub fn new(
        scan_number: u32,
        ms_level: MsLevel,
        retention_time_sec: f64,
        precursor: Option<PrecursorInfo>,
        mz_array: Vec<f64>,
        intensity_array: Vec<f64>,
    ) -> Result<Self, SpectrumError> {
        let spectrum = Self {
            scan_number,
            ms_level,
            retention_time_sec,
            precursor,
            mz_array,
            intensity_array,
        };
        spectrum.validate()?;
        Ok(spectrum)
    }

    /// Validates internal invariants.
    ///
    /// Call this after deserialization to ensure data consistency.
    pub fn validate(&self) -> Result<(), SpectrumError> {
        if self.mz_array.len() != self.intensity_array.len() {
            return Err(SpectrumError::ArrayLengthMismatch {
                mz_len: self.mz_array.len(),
                intensity_len: self.intensity_array.len(),
            });
        }
        Ok(())
    }

    /// Returns the number of peaks in this spectrum.
    pub fn num_peaks(&self) -> usize {
        self.mz_array.len()
    }
}

// ---------------------------------------------------------------------------
// SpectrumFormat & SpectrumFileInfo
// ---------------------------------------------------------------------------

/// Supported spectrum file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
pub enum SpectrumFormat {
    /// mzML format (PSI standard).
    MzML,
    /// Mascot Generic Format.
    Mgf,
}

impl std::fmt::Display for SpectrumFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpectrumFormat::MzML => write!(f, "mzML"),
            SpectrumFormat::Mgf => write!(f, "mgf"),
        }
    }
}

/// Metadata about a spectrum file on disk.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumFileInfo {
    /// File path.
    pub path: String,
    /// Detected file format.
    pub format: SpectrumFormat,
    /// File size in bytes.
    pub file_size_bytes: u64,
}

// ---------------------------------------------------------------------------
// SpectrumSummary
// ---------------------------------------------------------------------------

/// Statistical summary of a spectrum file.
///
/// This is the primary input for AI-driven parameter recommendation
/// and data quality assessment. The LLM reads this summary (via MCP tool)
/// to understand data characteristics before making recommendations.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumSummary {
    /// Source file path.
    pub file_path: String,
    /// File format.
    pub format: SpectrumFormat,
    /// Total number of spectra in the file.
    pub total_spectra: u64,
    /// Number of MS1 spectra.
    pub ms1_count: u64,
    /// Number of MS2 spectra.
    pub ms2_count: u64,
    /// m/z range as (min, max) in Da.
    pub mz_range: (f64, f64),
    /// Retention time range as (min, max) in seconds.
    pub rt_range_sec: (f64, f64),
    /// Distribution of precursor charge states (charge → count).
    pub precursor_charge_distribution: HashMap<i32, u64>,
    /// Median number of peaks per spectrum.
    pub median_peaks_per_spectrum: u32,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_precursor() -> PrecursorInfo {
        PrecursorInfo {
            mz: 500.2534,
            charge: Some(2),
            intensity: Some(1.5e6),
        }
    }

    fn sample_ms2_spectrum() -> Spectrum {
        Spectrum::new(
            42,
            MsLevel::MS2,
            120.5,
            Some(sample_precursor()),
            vec![100.0, 200.0, 300.0, 400.0],
            vec![1000.0, 2000.0, 500.0, 750.0],
        )
        .expect("test data should be valid")
    }

    fn sample_summary() -> SpectrumSummary {
        let mut charge_dist = HashMap::new();
        charge_dist.insert(2, 5000);
        charge_dist.insert(3, 3000);
        charge_dist.insert(4, 500);

        SpectrumSummary {
            file_path: "/data/sample.mzML".to_string(),
            format: SpectrumFormat::MzML,
            total_spectra: 12345,
            ms1_count: 2345,
            ms2_count: 10000,
            mz_range: (100.0, 2000.0),
            rt_range_sec: (0.0, 3600.0),
            precursor_charge_distribution: charge_dist,
            median_peaks_per_spectrum: 150,
        }
    }

    // -- MsLevel --------------------------------------------------------

    #[test]
    fn ms_level_serde_roundtrip() {
        for level in [MsLevel::MS1, MsLevel::MS2, MsLevel::Other(3)] {
            let json = serde_json::to_string(&level).unwrap();
            let back: MsLevel = serde_json::from_str(&json).unwrap();
            assert_eq!(level, back);
        }
    }

    #[test]
    fn ms_level_display() {
        assert_eq!(MsLevel::MS1.to_string(), "MS1");
        assert_eq!(MsLevel::MS2.to_string(), "MS2");
        assert_eq!(MsLevel::Other(5).to_string(), "MS5");
    }

    // -- PrecursorInfo --------------------------------------------------

    #[test]
    fn precursor_info_serde_roundtrip() {
        let info = sample_precursor();
        let json = serde_json::to_string_pretty(&info).unwrap();
        let back: PrecursorInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(info, back);
    }

    #[test]
    fn precursor_info_optional_fields_none() {
        let info = PrecursorInfo {
            mz: 400.0,
            charge: None,
            intensity: None,
        };
        let json = serde_json::to_string(&info).unwrap();
        let back: PrecursorInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(info, back);
    }

    // -- Spectrum -------------------------------------------------------

    #[test]
    fn spectrum_ms2_serde_roundtrip() {
        let spectrum = sample_ms2_spectrum();
        let json = serde_json::to_string_pretty(&spectrum).unwrap();
        let back: Spectrum = serde_json::from_str(&json).unwrap();
        assert_eq!(spectrum.scan_number, back.scan_number);
        assert_eq!(spectrum.ms_level, back.ms_level);
        assert_eq!(spectrum.mz_array, back.mz_array);
        assert_eq!(spectrum.intensity_array, back.intensity_array);
    }

    #[test]
    fn spectrum_ms1_no_precursor() {
        let spectrum = Spectrum::new(
            1,
            MsLevel::MS1,
            60.0,
            None,
            vec![200.0, 400.0],
            vec![1e5, 2e5],
        )
        .unwrap();
        let json = serde_json::to_string(&spectrum).unwrap();
        let back: Spectrum = serde_json::from_str(&json).unwrap();
        assert!(back.precursor.is_none());
        assert_eq!(back.ms_level, MsLevel::MS1);
    }

    #[test]
    fn spectrum_num_peaks() {
        let s = sample_ms2_spectrum();
        assert_eq!(s.num_peaks(), 4);
    }

    // -- SpectrumFormat -------------------------------------------------

    #[test]
    fn spectrum_format_serde_roundtrip() {
        for fmt in [SpectrumFormat::MzML, SpectrumFormat::Mgf] {
            let json = serde_json::to_string(&fmt).unwrap();
            let back: SpectrumFormat = serde_json::from_str(&json).unwrap();
            assert_eq!(fmt, back);
        }
    }

    #[test]
    fn spectrum_format_display() {
        assert_eq!(SpectrumFormat::MzML.to_string(), "mzML");
        assert_eq!(SpectrumFormat::Mgf.to_string(), "mgf");
    }

    // -- SpectrumFileInfo -----------------------------------------------

    #[test]
    fn spectrum_file_info_serde_roundtrip() {
        let info = SpectrumFileInfo {
            path: "/data/sample.mgf".to_string(),
            format: SpectrumFormat::Mgf,
            file_size_bytes: 1_048_576,
        };
        let json = serde_json::to_string_pretty(&info).unwrap();
        let back: SpectrumFileInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(info, back);
    }

    // -- SpectrumSummary ------------------------------------------------

    #[test]
    fn spectrum_summary_serde_roundtrip() {
        let summary = sample_summary();
        let json = serde_json::to_string_pretty(&summary).unwrap();
        let back: SpectrumSummary = serde_json::from_str(&json).unwrap();
        assert_eq!(summary.total_spectra, back.total_spectra);
        assert_eq!(summary.ms1_count, back.ms1_count);
        assert_eq!(summary.ms2_count, back.ms2_count);
        assert_eq!(summary.mz_range, back.mz_range);
        assert_eq!(summary.rt_range_sec, back.rt_range_sec);
        assert_eq!(
            summary.median_peaks_per_spectrum,
            back.median_peaks_per_spectrum
        );
        assert_eq!(
            summary.precursor_charge_distribution.len(),
            back.precursor_charge_distribution.len()
        );
        assert_eq!(
            summary.precursor_charge_distribution[&2],
            back.precursor_charge_distribution[&2]
        );
    }

    #[test]
    fn spectrum_summary_empty_charge_distribution() {
        let summary = SpectrumSummary {
            file_path: "/data/empty.mgf".to_string(),
            format: SpectrumFormat::Mgf,
            total_spectra: 0,
            ms1_count: 0,
            ms2_count: 0,
            mz_range: (0.0, 0.0),
            rt_range_sec: (0.0, 0.0),
            precursor_charge_distribution: HashMap::new(),
            median_peaks_per_spectrum: 0,
        };
        let json = serde_json::to_string(&summary).unwrap();
        let back: SpectrumSummary = serde_json::from_str(&json).unwrap();
        assert!(back.precursor_charge_distribution.is_empty());
    }

    // -- Validation -----------------------------------------------------

    #[test]
    fn spectrum_new_rejects_mismatched_arrays() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            None,
            vec![100.0, 200.0, 300.0],
            vec![1000.0, 2000.0], // one fewer
        );
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(
            err.to_string().contains("3"),
            "error should mention mz_len=3"
        );
        assert!(
            err.to_string().contains("2"),
            "error should mention intensity_len=2"
        );
    }

    #[test]
    fn spectrum_validate_catches_deserialized_bad_data() {
        // Simulate bad data arriving via JSON deserialization (bypasses new())
        let bad_json = r#"{
            "scan_number": 1,
            "ms_level": "MS2",
            "retention_time_sec": 10.0,
            "precursor": null,
            "mz_array": [100.0, 200.0],
            "intensity_array": [1000.0]
        }"#;
        let spectrum: Spectrum = serde_json::from_str(bad_json).unwrap();
        assert!(spectrum.validate().is_err());
    }

    #[test]
    fn spectrum_validate_passes_for_valid_data() {
        let s = sample_ms2_spectrum();
        assert!(s.validate().is_ok());
    }

    #[test]
    fn spectrum_new_accepts_empty_arrays() {
        let result = Spectrum::new(1, MsLevel::MS1, 0.0, None, vec![], vec![]);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().num_peaks(), 0);
    }
}
