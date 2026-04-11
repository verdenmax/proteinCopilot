//! Data types for the unified annotation + XIC HTML view.

use protein_copilot_search_engine::annotate::SpectrumAnnotation;
use protein_copilot_xic::{IonType, XicData};
use serde::Serialize;

/// Combined data for the unified HTML template.
#[derive(Debug, Clone, Serialize)]
pub struct UnifiedViewData {
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

/// Raw peak data from scans in the XIC RT window.
#[derive(Debug, Clone, Serialize)]
pub struct RawScanData {
    /// MS1 scans (trimmed to narrow m/z window around precursor).
    pub ms1_scans: Vec<RawScan>,
    /// MS2 scans (full peak lists).
    pub ms2_scans: Vec<RawScan>,
}

/// A single raw scan's peak list.
#[derive(Debug, Clone, Serialize)]
pub struct RawScan {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in seconds.
    pub retention_time_sec: f64,
    /// m/z values (sorted ascending).
    pub mz_array: Vec<f64>,
    /// Intensity values (parallel to mz_array).
    pub intensity_array: Vec<f64>,
}

/// Metadata for one fragment ion enabling client-side SILAC calculation.
#[derive(Debug, Clone, Serialize)]
pub struct IonMetadataEntry {
    /// Human-readable ion label (e.g. "y5¹⁺").
    pub label: String,
    /// Ion type.
    pub ion_type: IonType,
    /// Ion number (e.g. 5 for y5).
    pub ion_number: u32,
    /// Charge state.
    pub charge: u32,
    /// Theoretical m/z of the light (unlabeled) ion.
    pub light_mz: f64,
    /// Count of K (Lysine) residues in this fragment.
    pub k_count: u32,
    /// Count of R (Arginine) residues in this fragment.
    pub r_count: u32,
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
    fn ion_metadata_serializes() {
        let entry = IonMetadataEntry {
            label: "y5¹⁺".to_string(),
            ion_type: IonType::Y,
            ion_number: 5,
            charge: 1,
            light_mz: 574.28,
            k_count: 1,
            r_count: 0,
        };
        let json = serde_json::to_string(&entry).unwrap();
        assert!(json.contains("\"k_count\":1"));
        assert!(json.contains("\"light_mz\":574.28"));
    }

    #[test]
    fn raw_scan_serializes() {
        let scan = RawScan {
            scan_number: 42,
            retention_time_sec: 120.5,
            mz_array: vec![100.0, 200.0, 300.0],
            intensity_array: vec![1000.0, 5000.0, 2000.0],
        };
        let json = serde_json::to_string(&scan).unwrap();
        assert!(json.contains("\"scan_number\":42"));
        assert!(json.contains("[100.0,200.0,300.0]"));
    }

    #[test]
    fn raw_scan_data_serializes() {
        let data = RawScanData {
            ms1_scans: vec![RawScan {
                scan_number: 1,
                retention_time_sec: 10.0,
                mz_array: vec![450.0],
                intensity_array: vec![1000.0],
            }],
            ms2_scans: vec![],
        };
        let json = serde_json::to_string(&data).unwrap();
        assert!(json.contains("\"ms1_scans\""));
        assert!(json.contains("\"ms2_scans\":[]"));
    }
}
