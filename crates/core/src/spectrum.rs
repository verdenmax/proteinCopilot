//! Mass spectrometry spectrum data structures.
//!
//! This module defines the core types for representing mass spectrometry data:
//! - [`MsLevel`] — MS acquisition level (MS1, MS2, etc.)
//! - [`AcquisitionMode`] — Data acquisition mode (DDA, DIA, Unknown)
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

    /// A numeric field contains NaN or Infinity.
    #[error("{field} contains non-finite value")]
    NonFiniteValue {
        /// Name of the field with the invalid value.
        field: &'static str,
    },

    /// mz_array is not sorted in ascending order.
    #[error("mz_array is not sorted in ascending order")]
    MzArrayNotSorted,

    /// A range field has min > max.
    #[error("{field} has min ({min}) > max ({max})")]
    InvalidRange {
        /// Name of the range field.
        field: &'static str,
        /// The min value.
        min: f64,
        /// The max value.
        max: f64,
    },

    /// A field that must be strictly positive contains a zero or negative value.
    #[error("{field} contains non-positive value (must be > 0)")]
    NonPositiveValue {
        /// Name of the field with the invalid value.
        field: &'static str,
    },

    /// A field that must be non-negative contains a negative value.
    #[error("{field} contains negative value (must be ≥ 0)")]
    NegativeValue {
        /// Name of the field with the invalid value.
        field: &'static str,
    },

    /// Scan number is zero, violating 1-based indexing convention.
    #[error("scan_number must be ≥ 1 (1-based indexing)")]
    ZeroScanNumber,
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
// AcquisitionMode
// ---------------------------------------------------------------------------

/// Data acquisition mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[allow(clippy::upper_case_acronyms)]
pub enum AcquisitionMode {
    /// Data-Dependent Acquisition — narrow isolation window, single precursor.
    DDA,
    /// Data-Independent Acquisition — wide isolation window, multiple co-fragmented precursors.
    DIA,
    /// Acquisition mode could not be determined.
    Unknown,
}

impl std::fmt::Display for AcquisitionMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AcquisitionMode::DDA => write!(f, "DDA"),
            AcquisitionMode::DIA => write!(f, "DIA"),
            AcquisitionMode::Unknown => write!(f, "Unknown"),
        }
    }
}

// ---------------------------------------------------------------------------
// IsolationWindow
// ---------------------------------------------------------------------------

/// Isolation window used during precursor selection.
///
/// In DDA, the window is typically narrow (1–3 Da).
/// In DIA, the window is wide (e.g., 25 Da), covering many co-eluting ions.
/// Aligns with the mzML `<isolationWindow>` element.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct IsolationWindow {
    /// Center m/z of the isolation window.
    pub target_mz: f64,
    /// Offset below `target_mz` (m/z units, ≥ 0).
    pub lower_offset: f64,
    /// Offset above `target_mz` (m/z units, ≥ 0).
    pub upper_offset: f64,
}

// ---------------------------------------------------------------------------
// PrecursorInfo
// ---------------------------------------------------------------------------

/// Precursor ion information for MS2+ spectra.
///
/// In DDA, a spectrum typically has one precursor with a specific m/z and charge.
/// In DIA, the precursor describes the isolation window; m/z may be the window
/// center and charge may be unknown.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct PrecursorInfo {
    /// Precursor m/z value (mass-to-charge ratio).
    pub mz: f64,
    /// Charge state (`None` if undetermined, common in DIA).
    pub charge: Option<i32>,
    /// Intensity in detector counts (`None` if not available).
    pub intensity: Option<f64>,
    /// Isolation window (`None` if not recorded).
    pub isolation_window: Option<IsolationWindow>,
    /// MS1 scan number referenced by this precursor (from mzML `spectrumRef`).
    #[serde(default)]
    pub source_scan: Option<u32>,
}

// ---------------------------------------------------------------------------
// Spectrum
// ---------------------------------------------------------------------------

/// A single mass spectrum with peak data.
///
/// `mz_array` and `intensity_array` must always have the same length.
/// `mz_array` is expected to be sorted in ascending order.
///
/// The `precursors` field supports both DDA (typically 1 precursor)
/// and DIA (0 precursors, or 1 with a wide isolation window).
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct Spectrum {
    /// Scan number (1-based, following mass spectrometry convention).
    pub scan_number: u32,
    /// MS level of this spectrum.
    pub ms_level: MsLevel,
    /// Retention time in seconds.
    pub retention_time_min: f64,
    /// Precursor information. DDA: typically 1 entry. DIA: 0 or 1 with wide
    /// isolation window. MS1: empty.
    pub precursors: Vec<PrecursorInfo>,
    /// Array of m/z values (mass-to-charge ratio), sorted ascending.
    pub mz_array: Vec<f64>,
    /// Array of intensity values (detector counts), same length as `mz_array`.
    pub intensity_array: Vec<f64>,
}

impl Spectrum {
    /// Creates a new `Spectrum` with validation.
    ///
    /// Returns `SpectrumError` if any invariant is violated.
    pub fn new(
        scan_number: u32,
        ms_level: MsLevel,
        retention_time_min: f64,
        precursors: Vec<PrecursorInfo>,
        mz_array: Vec<f64>,
        intensity_array: Vec<f64>,
    ) -> Result<Self, SpectrumError> {
        let spectrum = Self {
            scan_number,
            ms_level,
            retention_time_min,
            precursors,
            mz_array,
            intensity_array,
        };
        spectrum.validate()?;
        Ok(spectrum)
    }

    /// Validates internal invariants.
    ///
    /// Call this after deserialization to ensure data consistency.
    /// Checks:
    /// - `scan_number` ≥ 1 (1-based indexing per mass spectrometry convention)
    /// - `mz_array` and `intensity_array` have the same length
    /// - All numeric fields are finite (no NaN or Infinity)
    /// - All precursor fields are finite; isolation window offsets ≥ 0
    /// - `mz_array` is sorted in ascending order
    pub fn validate(&self) -> Result<(), SpectrumError> {
        if self.scan_number == 0 {
            return Err(SpectrumError::ZeroScanNumber);
        }
        if self.mz_array.len() != self.intensity_array.len() {
            return Err(SpectrumError::ArrayLengthMismatch {
                mz_len: self.mz_array.len(),
                intensity_len: self.intensity_array.len(),
            });
        }
        if !self.retention_time_min.is_finite() {
            return Err(SpectrumError::NonFiniteValue {
                field: "retention_time_min",
            });
        }
        if self.retention_time_min < 0.0 {
            return Err(SpectrumError::NegativeValue {
                field: "retention_time_min",
            });
        }
        for (i, p) in self.precursors.iter().enumerate() {
            self.validate_precursor(i, p)?;
        }
        if self.mz_array.iter().any(|v| !v.is_finite()) {
            return Err(SpectrumError::NonFiniteValue { field: "mz_array" });
        }
        if self.mz_array.iter().any(|v| *v <= 0.0) {
            return Err(SpectrumError::NonPositiveValue { field: "mz_array" });
        }
        if self.intensity_array.iter().any(|v| !v.is_finite()) {
            return Err(SpectrumError::NonFiniteValue {
                field: "intensity_array",
            });
        }
        if self.intensity_array.iter().any(|v| *v < 0.0) {
            return Err(SpectrumError::NegativeValue {
                field: "intensity_array",
            });
        }
        if !self.mz_array.windows(2).all(|w| w[0] <= w[1]) {
            return Err(SpectrumError::MzArrayNotSorted);
        }
        Ok(())
    }

    fn validate_precursor(&self, idx: usize, p: &PrecursorInfo) -> Result<(), SpectrumError> {
        let mz_field = if idx == 0 {
            "precursors[0].mz"
        } else {
            "precursors[n].mz"
        };
        let intensity_field = if idx == 0 {
            "precursors[0].intensity"
        } else {
            "precursors[n].intensity"
        };

        if !p.mz.is_finite() {
            return Err(SpectrumError::NonFiniteValue { field: mz_field });
        }
        if p.mz <= 0.0 {
            return Err(SpectrumError::NonPositiveValue { field: mz_field });
        }
        if let Some(intensity) = p.intensity {
            if !intensity.is_finite() {
                return Err(SpectrumError::NonFiniteValue {
                    field: intensity_field,
                });
            }
            if intensity < 0.0 {
                return Err(SpectrumError::NegativeValue {
                    field: intensity_field,
                });
            }
        }
        if let Some(ref w) = p.isolation_window {
            if !w.target_mz.is_finite()
                || !w.lower_offset.is_finite()
                || !w.upper_offset.is_finite()
            {
                return Err(SpectrumError::NonFiniteValue {
                    field: "isolation_window",
                });
            }
            if w.target_mz <= 0.0 {
                return Err(SpectrumError::NonPositiveValue {
                    field: "isolation_window.target_mz",
                });
            }
            if w.lower_offset < 0.0 {
                return Err(SpectrumError::NegativeValue {
                    field: "isolation_window.lower_offset",
                });
            }
            if w.upper_offset < 0.0 {
                return Err(SpectrumError::NegativeValue {
                    field: "isolation_window.upper_offset",
                });
            }
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
    /// m/z range: \[min, max\].
    pub mz_range: [f64; 2],
    /// Retention time range: \[min, max\] in seconds.
    pub rt_range_min: [f64; 2],
    /// Distribution of precursor charge states (charge → count).
    pub precursor_charge_distribution: HashMap<i32, u64>,
    /// Median number of peaks per spectrum.
    pub median_peaks_per_spectrum: u32,
    /// Median isolation window width in Da (`None` if no isolation windows found).
    /// Useful for DIA detection: DDA windows are typically < 3 Da, DIA > 5 Da.
    #[serde(default)]
    pub median_isolation_window_da: Option<f64>,
}

impl SpectrumSummary {
    /// Returns `true` if the summary represents an empty file (no spectra).
    ///
    /// When `is_empty()` is true, `mz_range` and `rt_range_min` are set to
    /// `(0.0, 0.0)` as sentinel values and should not be interpreted as
    /// real data ranges.
    pub fn is_empty(&self) -> bool {
        self.total_spectra == 0
    }

    /// Validates that all numeric fields are finite and ranges are consistent.
    ///
    /// Checks:
    /// - `mz_range` values are finite and min ≤ max
    /// - `rt_range_min` values are finite and min ≤ max
    pub fn validate(&self) -> Result<(), SpectrumError> {
        if !self.mz_range[0].is_finite() || !self.mz_range[1].is_finite() {
            return Err(SpectrumError::NonFiniteValue { field: "mz_range" });
        }
        if !self.rt_range_min[0].is_finite() || !self.rt_range_min[1].is_finite() {
            return Err(SpectrumError::NonFiniteValue {
                field: "rt_range_min",
            });
        }
        if self.mz_range[0] > self.mz_range[1] {
            return Err(SpectrumError::InvalidRange {
                field: "mz_range",
                min: self.mz_range[0],
                max: self.mz_range[1],
            });
        }
        if self.rt_range_min[0] > self.rt_range_min[1] {
            return Err(SpectrumError::InvalidRange {
                field: "rt_range_min",
                min: self.rt_range_min[0],
                max: self.rt_range_min[1],
            });
        }
        if let Some(w) = self.median_isolation_window_da {
            if !w.is_finite() {
                return Err(SpectrumError::NonFiniteValue {
                    field: "median_isolation_window_da",
                });
            }
        }
        Ok(())
    }
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
            isolation_window: None,
            source_scan: None,
        }
    }

    fn sample_ms2_spectrum() -> Spectrum {
        Spectrum::new(
            42,
            MsLevel::MS2,
            120.5,
            vec![sample_precursor()],
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
            mz_range: [100.0, 2000.0],
            rt_range_min: [0.0, 60.0],
            precursor_charge_distribution: charge_dist,
            median_peaks_per_spectrum: 150,
            median_isolation_window_da: None,
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
            isolation_window: None,
            source_scan: None,
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
            vec![],
            vec![200.0, 400.0],
            vec![1e5, 2e5],
        )
        .unwrap();
        let json = serde_json::to_string(&spectrum).unwrap();
        let back: Spectrum = serde_json::from_str(&json).unwrap();
        assert!(back.precursors.is_empty());
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
        assert_eq!(summary.rt_range_min, back.rt_range_min);
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
            mz_range: [0.0, 0.0],
            rt_range_min: [0.0, 0.0],
            precursor_charge_distribution: HashMap::new(),
            median_peaks_per_spectrum: 0,
            median_isolation_window_da: None,
        };
        let json = serde_json::to_string(&summary).unwrap();
        let back: SpectrumSummary = serde_json::from_str(&json).unwrap();
        assert!(back.precursor_charge_distribution.is_empty());
    }

    #[test]
    fn spectrum_summary_validate_passes_for_valid_data() {
        let s = sample_summary();
        assert!(s.validate().is_ok());
    }

    #[test]
    fn spectrum_summary_validate_rejects_nan_mz_range() {
        let mut s = sample_summary();
        s.mz_range = [f64::NAN, 2000.0];
        assert!(s.validate().is_err());
        assert!(s.validate().unwrap_err().to_string().contains("mz_range"));
    }

    #[test]
    fn spectrum_summary_validate_rejects_infinity_rt_range() {
        let mut s = sample_summary();
        s.rt_range_min = [0.0, f64::INFINITY];
        assert!(s.validate().is_err());
        assert!(s
            .validate()
            .unwrap_err()
            .to_string()
            .contains("rt_range_min"));
    }

    #[test]
    fn spectrum_summary_validate_rejects_inverted_mz_range() {
        let mut s = sample_summary();
        s.mz_range = [2000.0, 100.0]; // min > max
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("mz_range"));
        assert!(err.to_string().contains("2000"));
    }

    #[test]
    fn spectrum_summary_validate_rejects_inverted_rt_range() {
        let mut s = sample_summary();
        s.rt_range_min = [60.0, 0.0]; // min > max
        let err = s.validate().unwrap_err();
        assert!(err.to_string().contains("rt_range_min"));
    }

    #[test]
    fn spectrum_summary_validate_accepts_equal_range() {
        let mut s = sample_summary();
        s.mz_range = [500.0, 500.0]; // min == max is OK (single value)
        s.rt_range_min = [100.0, 100.0];
        assert!(s.validate().is_ok());
    }

    #[test]
    fn spectrum_summary_is_empty() {
        let mut s = sample_summary();
        assert!(!s.is_empty());

        s.total_spectra = 0;
        assert!(s.is_empty());
    }

    // -- Validation -----------------------------------------------------

    #[test]
    fn spectrum_new_rejects_mismatched_arrays() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
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
            "retention_time_min": 10.0,
            "precursors": [],
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
        let result = Spectrum::new(1, MsLevel::MS1, 0.0, vec![], vec![], vec![]);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().num_peaks(), 0);
    }

    #[test]
    fn validate_rejects_zero_scan_number() {
        let result = Spectrum::new(0, MsLevel::MS2, 10.0, vec![], vec![100.0], vec![1000.0]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("scan_number"));
    }

    // -- NaN/Infinity validation ----------------------------------------

    #[test]
    fn validate_rejects_nan_retention_time() {
        let result = Spectrum::new(1, MsLevel::MS2, f64::NAN, vec![], vec![100.0], vec![1000.0]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("retention_time"));
    }

    #[test]
    fn validate_rejects_infinity_in_mz_array() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![100.0, f64::INFINITY],
            vec![1000.0, 2000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("mz_array"));
    }

    #[test]
    fn validate_rejects_nan_in_intensity_array() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![100.0, 200.0],
            vec![1000.0, f64::NAN],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("intensity_array"));
    }

    #[test]
    fn validate_rejects_nan_precursor_mz() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: f64::NAN,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("precursors"));
    }

    #[test]
    fn validate_rejects_infinity_precursor_intensity() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: Some(f64::INFINITY),
                isolation_window: None,
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("precursors"));
    }

    // -- Sortedness validation ------------------------------------------

    #[test]
    fn validate_rejects_unsorted_mz_array() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![300.0, 100.0, 200.0], // not sorted
            vec![1000.0, 2000.0, 500.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not sorted"));
    }

    #[test]
    fn validate_accepts_equal_adjacent_mz_values() {
        // Equal adjacent values (plateaus) are allowed
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![100.0, 100.0, 200.0],
            vec![1000.0, 2000.0, 500.0],
        );
        assert!(result.is_ok());
    }

    // -- DIA / IsolationWindow ------------------------------------------

    #[test]
    fn isolation_window_serde_roundtrip() {
        let w = IsolationWindow {
            target_mz: 500.0,
            lower_offset: 12.5,
            upper_offset: 12.5,
        };
        let json = serde_json::to_string_pretty(&w).unwrap();
        let back: IsolationWindow = serde_json::from_str(&json).unwrap();
        assert_eq!(w, back);
    }

    #[test]
    fn dia_spectrum_with_wide_isolation_window() {
        let spectrum = Spectrum::new(
            100,
            MsLevel::MS2,
            300.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: None, // DIA: charge often unknown
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: 500.0,
                    lower_offset: 12.5,
                    upper_offset: 12.5, // 25 Da window
                }),
                source_scan: None,
            }],
            vec![100.0, 200.0, 300.0],
            vec![500.0, 1000.0, 750.0],
        )
        .unwrap();
        assert_eq!(spectrum.precursors.len(), 1);
        assert!(spectrum.precursors[0].charge.is_none());
        assert!(spectrum.precursors[0].isolation_window.is_some());
    }

    #[test]
    fn precursor_with_isolation_window_serde_roundtrip() {
        let p = PrecursorInfo {
            mz: 750.0,
            charge: Some(3),
            intensity: Some(2e6),
            isolation_window: Some(IsolationWindow {
                target_mz: 750.0,
                lower_offset: 1.0,
                upper_offset: 1.0,
            }),
            source_scan: None,
        };
        let json = serde_json::to_string_pretty(&p).unwrap();
        let back: PrecursorInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(p, back);
    }

    #[test]
    fn validate_rejects_negative_isolation_window_offset() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: None,
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: 500.0,
                    lower_offset: -1.0, // invalid
                    upper_offset: 12.5,
                }),
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
    }

    #[test]
    fn validate_rejects_nan_isolation_window_target() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: None,
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: f64::NAN,
                    lower_offset: 12.5,
                    upper_offset: 12.5,
                }),
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("isolation_window"));
    }

    #[test]
    fn validate_rejects_zero_mz_in_array() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![0.0, 200.0],
            vec![1000.0, 2000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("non-positive"));
    }

    #[test]
    fn validate_rejects_negative_mz_in_array() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![-100.0, 200.0],
            vec![1000.0, 2000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("non-positive"));
    }

    #[test]
    fn validate_rejects_negative_intensity_in_array() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![100.0, 200.0],
            vec![1000.0, -500.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("negative"));
    }

    #[test]
    fn validate_accepts_zero_intensity() {
        // Zero intensity is valid (no signal detected at that m/z)
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![],
            vec![100.0, 200.0],
            vec![0.0, 2000.0],
        );
        assert!(result.is_ok());
    }

    #[test]
    fn validate_rejects_zero_precursor_mz() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: 0.0,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("non-positive"));
    }

    #[test]
    fn validate_rejects_negative_precursor_intensity() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: Some(-1.0),
                isolation_window: None,
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("negative"));
    }

    #[test]
    fn validate_rejects_zero_isolation_window_target_mz() {
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: None,
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: 0.0,
                    lower_offset: 12.5,
                    upper_offset: 12.5,
                }),
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        );
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("non-positive"));
    }

    #[test]
    fn multiple_precursors_valid() {
        // Chimeric DDA spectrum with 2 co-fragmented precursors
        let result = Spectrum::new(
            1,
            MsLevel::MS2,
            10.0,
            vec![
                PrecursorInfo {
                    mz: 400.0,
                    charge: Some(2),
                    intensity: Some(1e6),
                    isolation_window: None,
                    source_scan: None,
                },
                PrecursorInfo {
                    mz: 401.5,
                    charge: Some(3),
                    intensity: Some(5e5),
                    isolation_window: None,
                    source_scan: None,
                },
            ],
            vec![100.0, 200.0],
            vec![1000.0, 2000.0],
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap().precursors.len(), 2);
    }
}
