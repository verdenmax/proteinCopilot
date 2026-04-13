//! Data types for the unified annotation + XIC HTML view.

use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use protein_copilot_xic::XicData;
use serde::Serialize;

// Re-export types that live in xic crate
pub use protein_copilot_xic::{IonMetadataEntry, RawScan, RawScanData};

/// Combined data for the unified HTML template.
#[derive(Debug, Clone, Serialize)]
pub struct UnifiedViewData {
    /// Source spectrum file name (e.g. "sample.mzML").
    pub source_file: String,
    /// Spectrum annotation (peaks, coverage, metadata).
    pub annotation: SpectrumAnnotation,
    /// Pre-computed XIC data (light + optional heavy traces).
    pub xic: Option<XicData>,
    /// Raw scan peak arrays for client-side SILAC recomputation.
    pub raw_scans: Option<RawScanData>,
    /// Fragment ion metadata with K/R counts for SILAC calculation.
    pub ion_metadata: Vec<IonMetadataEntry>,
    /// Peptide-level info for SILAC computation.
    pub peptide_info: PeptideInfo,
}

/// Peptide-level SILAC info.
#[derive(Debug, Clone, Serialize)]
pub struct PeptideInfo {
    /// Peptide amino acid sequence.
    pub sequence: String,
    /// Precursor charge state.
    pub charge: i32,
    /// Light precursor m/z.
    pub precursor_mz: f64,
    /// Total K count in the full peptide.
    pub total_k: u32,
    /// Total R count in the full peptide.
    pub total_r: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_info_serializes() {
        let info = PeptideInfo {
            sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            total_k: 1,
            total_r: 0,
        };
        let json = serde_json::to_string(&info).unwrap();
        assert!(json.contains("\"total_k\":1"));
        assert!(json.contains("\"total_r\":0"));
        assert!(json.contains("PEPTIDEK"));
    }

    #[test]
    fn raw_scan_data_serializes() {
        let data = RawScanData {
            ms1_scans: vec![RawScan {
                scan_number: 1,
                retention_time_min: 10.0,
                mz_array: vec![450.0],
                intensity_array: vec![1000.0],
            }],
            ms2_scans: vec![],
            ms2_heavy_scans: vec![],
        };
        let json = serde_json::to_string(&data).unwrap();
        assert!(json.contains("\"ms1_scans\""));
        assert!(json.contains("\"ms2_scans\":[]"));
    }
}
