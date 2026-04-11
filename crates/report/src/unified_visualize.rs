//! Unified annotation + XIC visualization — renders combined HTML
//! with interactive SILAC controls.

use std::fs;
use std::path::Path;

use protein_copilot_xic::PlotlyMode;

use crate::error::ReportError;
use crate::unified_types::UnifiedViewData;

const PLOTLY_CDN: &str = "https://cdn.plot.ly/plotly-2.35.2.min.js";
const TEMPLATE: &str = include_str!("../templates/unified.html");

/// Renders a [`UnifiedViewData`] into a standalone HTML file.
///
/// The template has two placeholders:
/// - `/*__UNIFIED_JSON__*/` → replaced with `window.__UNIFIED_DATA__ = {...};`
/// - `__PLOTLY_SRC__` → replaced with the Plotly.js CDN URL
pub fn render_unified_html(
    data: &UnifiedViewData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(data)
        .map_err(|e| ReportError::SerializationError(e.to_string()))?;
    let json = crate::escape_json_for_html(&json);

    let plotly_src = match plotly_mode {
        PlotlyMode::Cdn | PlotlyMode::Embedded => PLOTLY_CDN.to_string(),
    };

    let html = TEMPLATE
        .replace(
            "/*__UNIFIED_JSON__*/",
            &format!("window.__UNIFIED_DATA__ = {};", json),
        )
        .replace("__PLOTLY_SRC__", &plotly_src);

    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent).map_err(|e| ReportError::IoError {
            path: parent.to_path_buf(),
            detail: e.to_string(),
        })?;
    }

    fs::write(output_path, &html).map_err(|e| ReportError::IoError {
        path: output_path.to_path_buf(),
        detail: e.to_string(),
    })?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unified_types::{PeptideInfo, UnifiedViewData};
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
    use protein_copilot_search_engine::annotate::{
        AnnotatedPeak, IonAnnotation, IonType as AnnotIonType, SpectrumAnnotation, TheoreticalIon,
    };
    use protein_copilot_xic::{
        ExtractionParams, IntensityRule, IonMetadataEntry, IonType as XicIonType, PlotlyMode,
        RawScan, RawScanData, XicData, XicDataPoint, XicTrace,
    };

    fn sample_annotation() -> SpectrumAnnotation {
        SpectrumAnnotation {
            scan_number: 42,
            retention_time_min: 120.0,
            peptide_sequence: "PEPTIDEK".to_string(),
            charge: 2,
            precursor_mz: 450.25,
            theoretical_mz: 450.2510,
            delta_mass_ppm: -2.2,
            score: 0.75,
            matched_ions: 5,
            total_ions: 14,
            protein_accessions: vec!["P12345".to_string()],
            peaks: vec![
                AnnotatedPeak {
                    mz: 200.0,
                    intensity: 1000.0,
                    annotation: Some(IonAnnotation {
                        ion_type: AnnotIonType::B,
                        ion_number: 3,
                        charge: 1,
                        theoretical_mz: 200.001,
                        delta_mz: -0.001,
                        delta_ppm: -5.0,
                    }),
                },
                AnnotatedPeak {
                    mz: 400.0,
                    intensity: 500.0,
                    annotation: None,
                },
            ],
            b_ions: vec![TheoreticalIon {
                ion_type: AnnotIonType::B,
                number: 3,
                charge: 1,
                theoretical_mz: 200.001,
                matched: true,
                matched_mz: Some(200.0),
                delta_ppm: Some(-5.0),
            }],
            y_ions: vec![TheoreticalIon {
                ion_type: AnnotIonType::Y,
                number: 1,
                charge: 1,
                theoretical_mz: 148.06,
                matched: false,
                matched_mz: None,
                delta_ppm: None,
            }],
            modifications: vec![],
        }
    }

    fn sample_xic() -> XicData {
        XicData {
            peptide_sequence: "PEPTIDEK".to_string(),
            target_rt_min: 120.0,
            target_scan: 100,
            charge: 2,
            precursor_mz: 450.25,
            ms1_precursor_xic: Some(XicTrace {
                ion_label: "precursor".to_string(),
                ion_type: XicIonType::Precursor,
                ion_number: 0,
                charge: 2,
                theoretical_mz: 450.25,
                data_points: vec![
                    XicDataPoint {
                        retention_time_min: 115.0,
                        scan_number: 90,
                        intensity: 1000.0,
                    },
                    XicDataPoint {
                        retention_time_min: 120.0,
                        scan_number: 100,
                        intensity: 5000.0,
                    },
                ],
                is_heavy: false,
            }),
            ms1_heavy_precursor_xic: None,
            fragment_xic_traces: vec![],
            heavy_fragment_xic_traces: Vec::new(),
            extraction_params: ExtractionParams {
                mz_tolerance: MassTolerance {
                    value: 20.0,
                    unit: ToleranceUnit::Ppm,
                },
                n_cycles: 5,
                top_n_ions: 6,
                label_type: None,
                intensity_rule: IntensityRule::MaxInWindow,
            },
        }
    }

    fn sample_unified_data() -> UnifiedViewData {
        UnifiedViewData {
            source_file: "test_sample.mzML".to_string(),
            annotation: sample_annotation(),
            xic: Some(sample_xic()),
            raw_scans: Some(RawScanData {
                ms1_scans: vec![RawScan {
                    scan_number: 90,
                    retention_time_min: 115.0,
                    mz_array: vec![450.0, 450.25, 451.0],
                    intensity_array: vec![100.0, 5000.0, 200.0],
                }],
                ms2_scans: vec![],
            }),
            ion_metadata: vec![IonMetadataEntry {
                label: "b3".to_string(),
                ion_type: XicIonType::B,
                ion_number: 3,
                charge: 1,
                light_mz: 200.001,
                k_count: 0,
                r_count: 0,
            }],
            peptide_info: PeptideInfo {
                sequence: "PEPTIDEK".to_string(),
                charge: 2,
                precursor_mz: 450.25,
                total_k: 1,
                total_r: 0,
            },
        }
    }

    #[test]
    fn render_unified_creates_html_with_data() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("unified.html");
        let data = sample_unified_data();
        render_unified_html(&data, &out, PlotlyMode::Cdn).unwrap();
        let html = fs::read_to_string(&out).unwrap();
        assert!(
            html.contains("__UNIFIED_DATA__"),
            "missing unified data injection"
        );
        assert!(html.contains("plotly-2.35.2"), "missing Plotly CDN");
        assert!(html.contains("PEPTIDEK"), "missing peptide sequence");
        assert!(html.contains("P12345"), "missing protein accession");
    }

    #[test]
    fn render_unified_without_xic() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("no_xic.html");
        let mut data = sample_unified_data();
        data.xic = None;
        data.raw_scans = None;
        render_unified_html(&data, &out, PlotlyMode::Cdn).unwrap();
        let html = fs::read_to_string(&out).unwrap();
        assert!(html.contains("__UNIFIED_DATA__"));
        // XIC null → JS hides sections
        assert!(html.contains("\"xic\":null"));
    }

    #[test]
    fn render_unified_creates_parent_dirs() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("sub").join("deep").join("unified.html");
        render_unified_html(&sample_unified_data(), &out, PlotlyMode::Cdn).unwrap();
        assert!(out.exists());
    }
}
