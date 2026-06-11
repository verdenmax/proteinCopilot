//! 3D MS2 annotation visualization — renders [`Xic3dData`] to a self-contained
//! HTML file with Plotly.js (3D overview + per-scan b/y annotated spectra).

use std::fs;
use std::path::Path;

use protein_copilot_xic::PlotlyMode;

use crate::error::ReportError;
use crate::xic3d_types::Xic3dData;

const PLOTLY_CDN: &str = "https://cdn.plot.ly/plotly-2.35.2.min.js";
const TEMPLATE: &str = include_str!("../templates/xic3d.html");

/// Default number of non-matched peaks per scan drawn in the 3D overview.
const DEFAULT_MAX_PEAKS_PER_SCAN_3D: usize = 200;

/// Renders [`Xic3dData`] into a standalone HTML file with Plotly.js charts:
/// a 3D MS2 overview plus per-scan b/y annotated spectra.
///
/// `max_peaks_per_scan_3d` caps the number of non-matched peaks drawn per scan
/// in the 3D overview (display-only declutter; matched b/y peaks are always
/// kept, and the per-scan annotated spectra always show all peaks). Defaults
/// to 200 when `None`.
pub fn render_xic_3d(
    data: &Xic3dData,
    output_path: &Path,
    plotly_mode: PlotlyMode,
    max_peaks_per_scan_3d: Option<usize>,
) -> Result<(), ReportError> {
    let json =
        serde_json::to_string(data).map_err(|e| ReportError::SerializationError(e.to_string()))?;
    let json = crate::escape_json_for_html(&json);

    // MVP: both modes load Plotly from the CDN. `Embedded` is accepted for
    // API parity with the other renderers but does not yet inline plotly.js.
    let plotly_src = match plotly_mode {
        PlotlyMode::Cdn | PlotlyMode::Embedded => PLOTLY_CDN.to_string(),
    };

    let max_peaks = max_peaks_per_scan_3d.unwrap_or(DEFAULT_MAX_PEAKS_PER_SCAN_3D);

    let escaped_peptide = data
        .peptide_sequence
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;");

    // Inject the JSON payload LAST: the static placeholders are replaced first
    // so that user-controlled data inside the JSON (e.g. a `source_file` that
    // happens to contain a placeholder token like `__PLOTLY_SRC__`) cannot be
    // clobbered by a later replacement.
    let html = TEMPLATE
        .replace("/*__MAX_PEAKS_3D__*/", &max_peaks.to_string())
        .replace("__PLOTLY_SRC__", &plotly_src)
        .replace("__PEPTIDE_PLACEHOLDER__", &escaped_peptide)
        .replace("/*__XIC3D_JSON__*/", &json);

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
    use protein_copilot_search_engine::matching::{
        generate_b_ions_with_charge, generate_y_ions_with_charge,
    };
    use protein_copilot_xic::{PlotlyMode, RawScan};

    fn sample_data() -> crate::xic3d_types::Xic3dData {
        let pep = "PEPTIDEK";
        let b = generate_b_ions_with_charge(pep, &[], 1);
        let y = generate_y_ions_with_charge(pep, &[], 1);
        let mut mz: Vec<f64> = Vec::new();
        mz.extend(b.iter().take(3).copied());
        mz.extend(y.iter().take(3).copied());
        mz.push(180.0);
        mz.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let intensity = vec![1000.0; mz.len()];
        let scans = vec![
            RawScan {
                scan_number: 100,
                retention_time_min: 19.9,
                mz_array: mz.clone(),
                intensity_array: intensity.clone(),
            },
            RawScan {
                scan_number: 101,
                retention_time_min: 20.0,
                mz_array: mz,
                intensity_array: intensity,
            },
        ];
        let tol = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        crate::xic3d_build::build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol, 101, "sample.mzML")
            .unwrap()
    }

    #[test]
    fn render_creates_html_with_injected_data() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("xic3d.html");
        render_xic_3d(&sample_data(), &out, PlotlyMode::Cdn, Some(150)).unwrap();
        let html = std::fs::read_to_string(&out).unwrap();
        assert!(html.contains("PEPTIDEK"), "peptide missing");
        assert!(html.contains("plotly-2.35.2"), "plotly cdn missing");
        assert!(html.contains("\"scan_number\":101"), "scan data missing");
        assert!(
            html.contains("var MAX_PEAKS_3D = 150;"),
            "max peaks const missing"
        );
        assert!(
            !html.contains("__XIC3D_JSON__"),
            "JSON placeholder not replaced"
        );
        assert!(
            !html.contains("__PLOTLY_SRC__"),
            "plotly placeholder not replaced"
        );
        assert!(
            !html.contains("__MAX_PEAKS_3D__"),
            "max peaks placeholder not replaced"
        );
        assert!(
            !html.contains("__PEPTIDE_PLACEHOLDER__"),
            "peptide placeholder not replaced"
        );
    }

    #[test]
    fn render_defaults_max_peaks_when_none() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("xic3d2.html");
        render_xic_3d(&sample_data(), &out, PlotlyMode::Cdn, None).unwrap();
        let html = std::fs::read_to_string(&out).unwrap();
        assert!(
            html.contains("var MAX_PEAKS_3D = 200;"),
            "default max peaks (200) missing"
        );
    }

    #[test]
    fn render_creates_parent_dirs() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("a").join("b").join("xic3d.html");
        render_xic_3d(&sample_data(), &out, PlotlyMode::Cdn, None).unwrap();
        assert!(out.exists());
    }

    #[test]
    fn render_escapes_angle_brackets_in_data() {
        let pep = "PEPTIDEK";
        let b = generate_b_ions_with_charge(pep, &[], 1);
        let mut mz: Vec<f64> = b.iter().take(2).copied().collect();
        mz.push(180.0);
        mz.sort_by(|a, c| a.partial_cmp(c).unwrap());
        let intensity = vec![1000.0; mz.len()];
        let scans = vec![RawScan {
            scan_number: 1,
            retention_time_min: 1.0,
            mz_array: mz,
            intensity_array: intensity,
        }];
        let tol = MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        };
        let data = crate::xic3d_build::build_xic3d_data(
            &scans,
            pep,
            2,
            460.0,
            &[],
            &tol,
            1,
            "a<script>b.mzML",
        )
        .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("x.html");
        render_xic_3d(&data, &out, PlotlyMode::Cdn, None).unwrap();
        let html = std::fs::read_to_string(&out).unwrap();
        assert!(
            html.contains("a\\u003cscript\\u003eb.mzML"),
            "angle brackets in data not escaped"
        );
    }
}
