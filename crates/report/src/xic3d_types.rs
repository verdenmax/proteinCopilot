//! Data structures for the 3D MS2 annotation view (`extract_xic view=3d`).
//!
//! [`Xic3dData`] bundles, for one identified peptide, the per-scan b/y
//! annotation of every MS2 spectrum in the target isolation window's
//! ±n_cycles RT range. Built by `build_xic3d_data`
//! and rendered by `render_xic_3d`.

use protein_copilot_core::search_params::MassTolerance;
use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use serde::{Deserialize, Serialize};

/// Complete data for the 3D MS2 annotation view.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Xic3dData {
    /// Identified peptide sequence.
    pub peptide_sequence: String,
    /// Precursor charge state.
    pub charge: i32,
    /// Real precursor m/z (PSM-derived, not DIA window center).
    pub precursor_mz: f64,
    /// Scan number that identified the peptide.
    pub target_scan: u32,
    /// Source spectrum file name.
    pub source_file: String,
    /// Fragment m/z tolerance used for b/y matching.
    pub mz_tolerance: MassTolerance,
    /// Per-scan annotations (ascending by scan number).
    pub scans: Vec<Ms2ScanAnnotation>,
}

/// One MS2 spectrum's annotation entry within the 3D view.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ms2ScanAnnotation {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in minutes.
    pub retention_time_min: f64,
    /// Total number of peaks in this spectrum (= full `mz_array` length).
    pub total_peaks: usize,
    /// Whether this is the scan that identified the peptide.
    pub is_target: bool,
    /// Reused single-spectrum annotation: peaks (optionally b/y annotated),
    /// `b_ions`/`y_ions`, `matched_ions`/`total_ions`, etc.
    pub annotation: SpectrumAnnotation,
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

    #[test]
    fn xic3d_data_serde_round_trip() {
        let data = Xic3dData {
            peptide_sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            target_scan: 1852,
            source_file: "sample.mzML".to_string(),
            mz_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            scans: vec![],
        };
        let json = serde_json::to_string(&data).unwrap();
        let back: Xic3dData = serde_json::from_str(&json).unwrap();
        assert_eq!(back.peptide_sequence, "PEPTIDEK");
        assert_eq!(back.target_scan, 1852);
        assert_eq!(back.charge, 2);
        assert!(back.scans.is_empty());
    }
}
