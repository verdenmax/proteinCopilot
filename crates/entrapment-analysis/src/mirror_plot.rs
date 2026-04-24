//! Mirror Plot Renderer
//!
//! Generates a standalone HTML file containing a Plotly.js bar chart that
//! visualises fragment-ion provenance.  Peaks are coloured by their
//! [`IonOrigin`] classification:
//!
//! | Origin       | Colour  | Hex       |
//! |--------------|---------|-----------|
//! | TrapOnly     | blue    | `#1f77b4` |
//! | TargetOnly   | red     | `#d62728` |
//! | Shared       | purple  | `#9467bd` |
//! | Unassigned   | gray    | `#7f7f7f` |
//!
//! Ion labels (e.g. `b3+1`) are rendered as text annotations above each bar.

use crate::provenance::{FragmentProvenance, IonOrigin};

/// Render a mirror plot as a standalone HTML file.
///
/// The plot shows all observed fragment-ion peaks as upward bars coloured by
/// their provenance classification.  Ion labels are annotated above each bar.
pub fn render_mirror_plot(
    provenance: &FragmentProvenance,
    output_path: &std::path::Path,
) -> Result<(), std::io::Error> {
    let html = generate_mirror_html(provenance);
    std::fs::write(output_path, html)
}

/// Generate the HTML string for a mirror plot.
///
/// Returns a self-contained HTML document that loads Plotly.js from a CDN
/// and renders a bar chart of all annotated peaks.
pub fn generate_mirror_html(provenance: &FragmentProvenance) -> String {
    let trap_seq = &provenance.trap_sequence;
    let target_seq = &provenance.target_sequence;

    // Group peaks by origin for separate Plotly traces.
    let mut trap_only_mz = Vec::new();
    let mut trap_only_int = Vec::new();
    let mut trap_only_labels = Vec::new();

    let mut target_only_mz = Vec::new();
    let mut target_only_int = Vec::new();
    let mut target_only_labels = Vec::new();

    let mut shared_mz = Vec::new();
    let mut shared_int = Vec::new();
    let mut shared_labels = Vec::new();

    let mut unassigned_mz = Vec::new();
    let mut unassigned_int = Vec::new();
    let mut unassigned_labels: Vec<String> = Vec::new();

    for peak in &provenance.annotated_peaks {
        let label = peak
            .trap_ion_label
            .as_deref()
            .or(peak.target_ion_label.as_deref())
            .unwrap_or("")
            .to_string();

        match peak.origin {
            IonOrigin::TrapOnly => {
                trap_only_mz.push(peak.mz_observed);
                trap_only_int.push(peak.intensity);
                trap_only_labels.push(label);
            }
            IonOrigin::TargetOnly => {
                target_only_mz.push(peak.mz_observed);
                target_only_int.push(peak.intensity);
                target_only_labels.push(label);
            }
            IonOrigin::Shared => {
                shared_mz.push(peak.mz_observed);
                shared_int.push(peak.intensity);
                shared_labels.push(label);
            }
            IonOrigin::Unassigned => {
                unassigned_mz.push(peak.mz_observed);
                unassigned_int.push(peak.intensity);
                unassigned_labels.push(label);
            }
        }
    }

    let traces = build_traces(&[
        TraceData {
            name: "TrapOnly",
            color: "#1f77b4",
            mz: &trap_only_mz,
            intensity: &trap_only_int,
            labels: &trap_only_labels,
        },
        TraceData {
            name: "TargetOnly",
            color: "#d62728",
            mz: &target_only_mz,
            intensity: &target_only_int,
            labels: &target_only_labels,
        },
        TraceData {
            name: "Shared",
            color: "#9467bd",
            mz: &shared_mz,
            intensity: &shared_int,
            labels: &shared_labels,
        },
        TraceData {
            name: "Unassigned",
            color: "#7f7f7f",
            mz: &unassigned_mz,
            intensity: &unassigned_int,
            labels: &unassigned_labels,
        },
    ]);

    // Build text annotations for ion labels.
    let annotations = build_annotations(provenance);

    format!(
        r#"<!DOCTYPE html>
<html>
<head>
    <title>Mirror Plot: {trap_seq} vs {target_seq}</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
</head>
<body>
    <div id="mirror-plot" style="width:100%;height:600px;"></div>
    <script>
        var traces = [{traces}];
        var layout = {{
            title: 'Fragment Ion Provenance: {trap_seq} vs {target_seq}',
            xaxis: {{ title: 'm/z' }},
            yaxis: {{ title: 'Intensity' }},
            barmode: 'overlay',
            bargap: 0,
            showlegend: true,
            annotations: [{annotations}]
        }};
        Plotly.newPlot('mirror-plot', traces, layout);
    </script>
</body>
</html>"#,
        trap_seq = escape_js(trap_seq),
        target_seq = escape_js(target_seq),
        traces = traces,
        annotations = annotations,
    )
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Data for a single Plotly trace.
struct TraceData<'a> {
    name: &'a str,
    color: &'a str,
    mz: &'a [f64],
    intensity: &'a [f64],
    labels: &'a [String],
}

/// Build JavaScript array elements for Plotly traces.
fn build_traces(groups: &[TraceData<'_>]) -> String {
    let mut parts = Vec::new();
    for group in groups {
        if group.mz.is_empty() {
            continue;
        }
        let x_vals = format_f64_array(group.mz);
        let y_vals = format_f64_array(group.intensity);
        let text_vals = format_string_array(group.labels);
        parts.push(format!(
            r#"{{
                type: 'bar',
                name: '{name}',
                x: [{x}],
                y: [{y}],
                text: [{text}],
                width: 2,
                marker: {{ color: '{color}' }}
            }}"#,
            name = group.name,
            x = x_vals,
            y = y_vals,
            text = text_vals,
            color = group.color,
        ));
    }
    parts.join(",\n")
}

/// Build Plotly annotation objects for ion labels above bars.
fn build_annotations(provenance: &FragmentProvenance) -> String {
    let mut parts = Vec::new();
    for peak in &provenance.annotated_peaks {
        let label = peak
            .trap_ion_label
            .as_deref()
            .or(peak.target_ion_label.as_deref())
            .unwrap_or("");
        if label.is_empty() {
            continue;
        }
        parts.push(format!(
            r#"{{
                x: {mz},
                y: {intensity},
                text: '{label}',
                showarrow: false,
                yshift: 10,
                font: {{ size: 10 }}
            }}"#,
            mz = peak.mz_observed,
            intensity = peak.intensity,
            label = escape_js(label),
        ));
    }
    parts.join(",\n")
}

/// Format a slice of `f64` as comma-separated values.
fn format_f64_array(vals: &[f64]) -> String {
    vals.iter()
        .map(|v| format!("{v}"))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Format a slice of strings as comma-separated quoted JavaScript strings.
fn format_string_array(vals: &[String]) -> String {
    vals.iter()
        .map(|s| format!("'{}'", escape_js(s)))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Minimal JavaScript string escaping.
fn escape_js(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('\'', "\\'")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace("</", "<\\/")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::{AnnotatedPeak, FragmentProvenance, IonOrigin};

    fn sample_provenance() -> FragmentProvenance {
        FragmentProvenance {
            trap_sequence: "PEPTIDE".into(),
            target_sequence: "PAPTIDE".into(),
            annotated_peaks: vec![
                AnnotatedPeak {
                    mz_observed: 200.1,
                    intensity: 1000.0,
                    origin: IonOrigin::TrapOnly,
                    trap_ion_label: Some("b2+1".into()),
                    target_ion_label: None,
                },
                AnnotatedPeak {
                    mz_observed: 300.2,
                    intensity: 800.0,
                    origin: IonOrigin::TargetOnly,
                    trap_ion_label: None,
                    target_ion_label: Some("b3+1".into()),
                },
                AnnotatedPeak {
                    mz_observed: 400.3,
                    intensity: 1200.0,
                    origin: IonOrigin::Shared,
                    trap_ion_label: Some("y4+1".into()),
                    target_ion_label: Some("y4+1".into()),
                },
                AnnotatedPeak {
                    mz_observed: 500.4,
                    intensity: 200.0,
                    origin: IonOrigin::Unassigned,
                    trap_ion_label: None,
                    target_ion_label: None,
                },
            ],
            trap_matched_count: 1,
            target_matched_count: 1,
            shared_count: 1,
            unassigned_count: 1,
            shared_ratio: 0.3333,
            is_chimeric: true,
        }
    }

    #[test]
    fn test_generate_mirror_html_contains_plotly() {
        let html = generate_mirror_html(&sample_provenance());
        assert!(html.contains("plotly"));
        assert!(html.contains("mirror-plot"));
    }

    #[test]
    fn test_generate_mirror_html_contains_sequences() {
        let html = generate_mirror_html(&sample_provenance());
        assert!(html.contains("PEPTIDE"));
        assert!(html.contains("PAPTIDE"));
    }

    #[test]
    fn test_generate_mirror_html_contains_colors() {
        let html = generate_mirror_html(&sample_provenance());
        assert!(html.contains("#1f77b4") || html.contains("blue")); // TrapOnly
        assert!(html.contains("#d62728") || html.contains("red")); // TargetOnly
    }

    #[test]
    fn test_generate_mirror_html_contains_ion_labels() {
        let html = generate_mirror_html(&sample_provenance());
        assert!(html.contains("b2+1"));
        assert!(html.contains("b3+1"));
        assert!(html.contains("y4+1"));
    }

    #[test]
    fn test_render_mirror_plot_writes_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mirror.html");
        render_mirror_plot(&sample_provenance(), &path).unwrap();
        assert!(path.exists());
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("plotly"));
    }

    #[test]
    fn test_generate_mirror_html_empty_peaks() {
        let prov = FragmentProvenance {
            trap_sequence: "ABC".into(),
            target_sequence: "DEF".into(),
            annotated_peaks: vec![],
            trap_matched_count: 0,
            target_matched_count: 0,
            shared_count: 0,
            unassigned_count: 0,
            shared_ratio: 0.0,
            is_chimeric: false,
        };
        let html = generate_mirror_html(&prov);
        assert!(html.contains("plotly"));
        assert!(html.contains("ABC"));
    }
}
