//! # ProteinCopilot XIC
//!
//! Extracted Ion Chromatogram (XIC) computation for proteomics data.
//!
//! This crate handles:
//! - Fragment-ion-level XIC extraction from mzML files (DIA + DDA)
//! - MS1 precursor XIC extraction
//! - SILAC heavy-label m/z calculation for paired XIC traces
//!
//! It is a pure computation crate — no MCP, network, or LLM dependencies.

pub mod error;
pub mod extract;
pub mod heavy;

pub use error::XicError;

use protein_copilot_core::search_params::MassTolerance;
use serde::{Deserialize, Serialize};

/// Ion type for XIC trace identification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IonType {
    /// b-ion (N-terminal fragment).
    B,
    /// y-ion (C-terminal fragment).
    Y,
    /// Precursor ion (MS1 level).
    Precursor,
}

/// A single data point on an XIC trace.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicDataPoint {
    /// Retention time in seconds.
    pub retention_time_min: f64,
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Extracted intensity (0.0 if no peak found in tolerance window).
    pub intensity: f64,
    /// Observed m/z of the matched peak (`None` when no peak found).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub observed_mz: Option<f64>,
}

/// A single XIC trace for one ion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicTrace {
    /// Human-readable ion label (e.g. "y5¹⁺", "b3²⁺", "precursor").
    pub ion_label: String,
    /// Ion type.
    pub ion_type: IonType,
    /// Ion number (e.g. 5 for y5; 0 for precursor).
    pub ion_number: u32,
    /// Charge state.
    pub charge: u32,
    /// Theoretical m/z used for extraction.
    pub theoretical_mz: f64,
    /// Extracted data points.
    pub data_points: Vec<XicDataPoint>,
    /// Whether this is a heavy-label trace.
    pub is_heavy: bool,
}

/// Complete XIC extraction result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct XicData {
    /// Target peptide sequence.
    pub peptide_sequence: String,
    /// Target scan retention time (seconds).
    pub target_rt_min: f64,
    /// Target scan number.
    pub target_scan: u32,
    /// Precursor charge state.
    pub charge: i32,
    /// Real precursor m/z (from PSM or user input, not DIA window center).
    pub precursor_mz: f64,
    /// MS1 precursor XIC (light). `None` if MS1 data unavailable.
    pub ms1_precursor_xic: Option<XicTrace>,
    /// MS1 precursor XIC (heavy). `None` if no label or MS1 unavailable.
    pub ms1_heavy_precursor_xic: Option<XicTrace>,
    /// MS2 fragment ion XIC traces (light).
    pub fragment_xic_traces: Vec<XicTrace>,
    /// MS2 fragment ion XIC traces (heavy).
    pub heavy_fragment_xic_traces: Vec<XicTrace>,
    /// Parameters used for extraction.
    pub extraction_params: ExtractionParams,
}

/// Parameters controlling XIC extraction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtractionParams {
    /// Mass tolerance for peak matching.
    pub mz_tolerance: MassTolerance,
    /// Number of DIA cycles before/after target scan.
    pub n_cycles: u32,
    /// Number of top ions to display by default.
    pub top_n_ions: usize,
    /// Heavy-label type (None = no label).
    pub label_type: Option<LabelType>,
    /// How to extract intensity from peaks within tolerance.
    pub intensity_rule: IntensityRule,
}

/// Intensity extraction strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema)]
pub enum IntensityRule {
    /// Highest peak within tolerance window (default).
    #[default]
    MaxInWindow,
    /// Sum of all peaks within tolerance window.
    SumInWindow,
    /// Nearest peak to theoretical m/z.
    NearestPeak,
}

/// Heavy-label type for SILAC comparison.
#[derive(Debug, Clone, Serialize, Deserialize, schemars::JsonSchema)]
pub enum LabelType {
    /// SILAC heavy amino acids.
    Silac {
        /// Mass shift for heavy Lysine (default: 8.014199 Da for ¹³C₆¹⁵N₂-Lys).
        heavy_k_delta: f64,
        /// Mass shift for heavy Arginine (default: 10.008269 Da for ¹³C₆¹⁵N₄-Arg).
        heavy_r_delta: f64,
    },
    /// Custom residue mass shifts.
    Custom {
        /// (residue, mass_delta) pairs.
        residue_deltas: Vec<(char, f64)>,
    },
}

impl LabelType {
    /// Standard SILAC heavy labels (K+8, R+10).
    pub fn standard_silac() -> Self {
        LabelType::Silac {
            heavy_k_delta: 8.014199,
            heavy_r_delta: 10.008269,
        }
    }
}

/// Plotly.js loading mode for HTML output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, schemars::JsonSchema)]
pub enum PlotlyMode {
    /// Load from CDN (default, smaller file).
    #[default]
    Cdn,
    /// Embed plotly-basic.min.js inline (larger file, works offline).
    Embedded,
}

/// Metadata for one fragment ion enabling client-side SILAC calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
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

/// Raw peak data from scans in the XIC RT window.
/// Used for client-side SILAC recomputation in the unified HTML.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawScanData {
    /// MS1 scans (trimmed to narrow m/z window around precursor).
    pub ms1_scans: Vec<RawScan>,
    /// MS2 scans (full peak lists from matching isolation windows).
    pub ms2_scans: Vec<RawScan>,
}

/// A single raw scan's peak list.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawScan {
    /// Scan number (1-based).
    pub scan_number: u32,
    /// Retention time in seconds.
    pub retention_time_min: f64,
    /// m/z values (sorted ascending).
    pub mz_array: Vec<f64>,
    /// Intensity values (parallel to mz_array).
    pub intensity_array: Vec<f64>,
}
