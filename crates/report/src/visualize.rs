//! Spectrum annotation visualization — renders [`SpectrumAnnotation`] to a
//! self-contained HTML file with an interactive SVG spectrum plot, peptide
//! fragment-ion coverage diagram, and hover tooltips.

use std::fs;
use std::path::Path;

use protein_copilot_search_engine::annotate::SpectrumAnnotation;

use crate::error::ReportError;

/// Embedded HTML template (no external dependencies).
const TEMPLATE: &str = include_str!("../templates/annotation.html");

/// Renders a [`SpectrumAnnotation`] into a standalone HTML file.
///
/// The template contains a placeholder comment `/*__ANNOTATION_JSON__*/` that
/// is replaced with the serialized annotation data so the browser-side
/// JavaScript can render the spectrum plot.
pub fn render_annotation_html(
    annotation: &SpectrumAnnotation,
    output_path: &Path,
) -> Result<(), ReportError> {
    let json = serde_json::to_string(annotation)
        .map_err(|e| ReportError::SerializationError(e.to_string()))?;
    let json = crate::escape_json_for_html(&json);

    let html = TEMPLATE.replace(
        "/*__ANNOTATION_JSON__*/",
        &format!("window.__ANNOTATION_DATA__ = {};", json),
    );

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

    use protein_copilot_search_engine::annotate::{
        AnnotatedPeak, IonAnnotation, IonType, SpectrumAnnotation, TheoreticalIon,
    };

    fn sample_annotation() -> SpectrumAnnotation {
        SpectrumAnnotation {
            scan_number: 42,
            retention_time_min: 120.0,
            peptide_sequence: "PEPTIDE".to_string(),
            charge: 2,
            precursor_mz: 400.6932,
            theoretical_mz: 400.6940,
            delta_mass_ppm: -2.0,
            score: 0.75,
            matched_ions: 9,
            total_ions: 12,
            protein_accessions: vec!["P12345".to_string()],
            peaks: vec![
                AnnotatedPeak {
                    mz: 100.0,
                    intensity: 500.0,
                    annotation: None,
                },
                AnnotatedPeak {
                    mz: 200.0,
                    intensity: 1000.0,
                    annotation: Some(IonAnnotation {
                        ion_type: IonType::B,
                        ion_number: 3,
                        charge: 1,
                        theoretical_mz: 200.001,
                        delta_mz: -0.001,
                        delta_ppm: -5.0,
                    }),
                },
            ],
            b_ions: vec![TheoreticalIon {
                ion_type: IonType::B,
                number: 1,
                charge: 1,
                theoretical_mz: 98.06,
                matched: true,
                matched_mz: Some(98.061),
                delta_ppm: Some(10.2),
            }],
            y_ions: vec![TheoreticalIon {
                ion_type: IonType::Y,
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

    #[test]
    fn render_creates_html_with_expected_content() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("annotation.html");

        render_annotation_html(&sample_annotation(), &out).expect("render failed");

        let html = fs::read_to_string(&out).expect("read");
        assert!(
            html.contains("__ANNOTATION_DATA__"),
            "missing data injection"
        );
        assert!(html.contains("svgEl("), "missing SVG rendering logic");
        assert!(html.contains("PEPTIDE"), "missing peptide sequence");
        assert!(html.contains("P12345"), "missing protein accession");
    }

    #[test]
    fn render_creates_parent_directories() {
        let dir = tempfile::tempdir().expect("tempdir");
        let out = dir.path().join("sub").join("deep").join("annotation.html");

        render_annotation_html(&sample_annotation(), &out).expect("render failed");

        assert!(out.exists());
    }
}
