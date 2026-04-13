//! XIC visualization — renders [`XicData`] to a self-contained HTML file
//! with interactive Plotly.js charts.

use std::fs;
use std::path::Path;

use protein_copilot_xic::{PlotlyMode, XicData};

use crate::error::ReportError;

const PLOTLY_CDN: &str = "https://cdn.plot.ly/plotly-2.35.2.min.js";

const TEMPLATE: &str = include_str!("../templates/xic.html");

/// Renders XIC data into a standalone HTML file with Plotly.js charts.
pub fn render_xic_html(
    xic_data: &XicData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(xic_data)
        .map_err(|e| ReportError::SerializationError(e.to_string()))?;
    let json = crate::escape_json_for_html(&json);

    let plotly_src = match plotly_mode {
        PlotlyMode::Cdn => PLOTLY_CDN.to_string(),
        PlotlyMode::Embedded => {
            // For embedded mode, fall back to CDN for MVP.
            PLOTLY_CDN.to_string()
        }
    };

    let html = TEMPLATE
        .replace("/*__XIC_JSON__*/", &json)
        .replace("__PLOTLY_SRC__", &plotly_src)
        .replace("__PEPTIDE_PLACEHOLDER__", &xic_data.peptide_sequence);

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
    use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};
    use protein_copilot_xic::*;

    fn sample_xic_data() -> XicData {
        XicData {
            peptide_sequence: "PEPTIDEK".to_string(),
            target_rt_min: 120.0,
            target_scan: 100,
            charge: 2,
            precursor_mz: 450.25,
            ms1_precursor_xic: Some(XicTrace {
                ion_label: "precursor".to_string(),
                ion_type: IonType::Precursor,
                ion_number: 0,
                charge: 2,
                theoretical_mz: 450.25,
                data_points: vec![
                    XicDataPoint { retention_time_min: 115.0, scan_number: 90, intensity: 1000.0, observed_mz: None },
                    XicDataPoint { retention_time_min: 120.0, scan_number: 100, intensity: 5000.0, observed_mz: None },
                    XicDataPoint { retention_time_min: 125.0, scan_number: 110, intensity: 2000.0, observed_mz: None },
                ],
                is_heavy: false,
            }),
            ms1_heavy_precursor_xic: None,
            fragment_xic_traces: vec![
                XicTrace {
                    ion_label: "y5\u{00b9}\u{207a}".to_string(),
                    ion_type: IonType::Y,
                    ion_number: 5,
                    charge: 1,
                    theoretical_mz: 574.28,
                    data_points: vec![
                        XicDataPoint { retention_time_min: 115.0, scan_number: 91, intensity: 800.0, observed_mz: None },
                        XicDataPoint { retention_time_min: 120.0, scan_number: 101, intensity: 3000.0, observed_mz: None },
                        XicDataPoint { retention_time_min: 125.0, scan_number: 111, intensity: 1200.0, observed_mz: None },
                    ],
                    is_heavy: false,
                },
            ],
            heavy_fragment_xic_traces: Vec::new(),
            extraction_params: ExtractionParams {
                mz_tolerance: MassTolerance { value: 20.0, unit: ToleranceUnit::Ppm },
                n_cycles: 5,
                top_n_ions: 6,
                label_type: None,
                intensity_rule: IntensityRule::MaxInWindow,
            },
        }
    }

    #[test]
    fn render_xic_html_creates_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test_xic.html");
        let data = sample_xic_data();
        render_xic_html(&data, &path, PlotlyMode::Cdn).unwrap();
        assert!(path.exists());
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("PEPTIDEK"));
        assert!(content.contains("plotly"));
        assert!(content.contains("y5"));
    }

    #[test]
    fn render_xic_html_contains_json_data() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test_xic2.html");
        let data = sample_xic_data();
        render_xic_html(&data, &path, PlotlyMode::Cdn).unwrap();
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("450.25"));
        assert!(content.contains("574.28"));
    }
}
