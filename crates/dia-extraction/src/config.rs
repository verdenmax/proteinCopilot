//! Configuration and result types for DIA precursor extraction.

use protein_copilot_core::spectrum::AcquisitionMode;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Configuration for DIA precursor extraction.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct DiaExtractionConfig {
    /// Acquisition mode override. `None` = auto-detect from isolation window widths.
    pub acquisition_mode: Option<AcquisitionMode>,
    /// Isolation window width threshold (Da) for DIA detection. Default: 5.0.
    pub dia_threshold_da: f64,
}

impl Default for DiaExtractionConfig {
    fn default() -> Self {
        Self {
            acquisition_mode: None,
            dia_threshold_da: 5.0,
        }
    }
}

/// Statistics from a DIA precursor extraction run.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ExtractionStats {
    /// Number of MS1 spectra used for extraction.
    pub ms1_count: u32,
    /// Number of MS2 spectra processed.
    pub ms2_count: u32,
    /// Total precursors extracted across all MS2 spectra.
    pub total_precursors_extracted: u32,
    /// Average number of precursors extracted per MS2 spectrum.
    pub avg_precursors_per_ms2: f64,
    /// Distribution of assigned charge states.
    pub charge_distribution: HashMap<i32, u32>,
}

/// Result of extracting precursors from a single MS2 spectrum.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SingleSpectrumExtractionResult {
    /// The MS2 scan number that was analyzed.
    pub ms2_scan: u32,
    /// The MS1 scan used for isotope pattern extraction.
    pub ms1_scan_used: u32,
    /// How the MS1 was selected: "source_scan", "scan_order", or "rt_nearest".
    pub correlation_method: String,
    /// Isolation window of the MS2 spectrum (if available).
    pub isolation_window: Option<protein_copilot_core::spectrum::IsolationWindow>,
    /// Extracted precursor candidates from isotope pattern analysis.
    pub precursors: Vec<protein_copilot_core::spectrum::PrecursorInfo>,
}

/// Result of DIA precursor extraction.
#[derive(Debug, Clone)]
pub struct DiaExtractionResult {
    /// Detected (or overridden) acquisition mode.
    pub detected_mode: AcquisitionMode,
    /// MS2 spectra with extracted precursors populated.
    pub enhanced_spectra: Vec<protein_copilot_core::spectrum::Spectrum>,
    /// Extraction statistics.
    pub stats: ExtractionStats,
}

impl DiaExtractionResult {
    /// Expand each multi-precursor spectrum into multiple single-precursor pseudo-spectra.
    /// Useful for compatibility with search engines that expect one precursor per spectrum.
    pub fn expand_to_pseudo_spectra(&self) -> Vec<protein_copilot_core::spectrum::Spectrum> {
        let mut pseudo = Vec::new();
        for spectrum in &self.enhanced_spectra {
            if spectrum.precursors.len() <= 1 {
                pseudo.push(spectrum.clone());
            } else {
                for precursor in &spectrum.precursors {
                    let mut ps = spectrum.clone();
                    ps.precursors = vec![precursor.clone()];
                    pseudo.push(ps);
                }
            }
        }
        pseudo
    }
}
