//! Per-PSM HTML Report Renderer for Multi-Target Fragment Provenance.
//!
//! Generates self-contained HTML reports that visualise which fragment ions
//! of a trap PSM originate from co-eluting target peptides.  Each report
//! contains:
//!
//! 1. **Candidate table** — co-eluting target peptides with m/z, charge, RT, etc.
//! 2. **Mirror spectrum** — Plotly.js bar chart: trap peaks up, target peaks down.
//! 3. **Attribution table** — per-peak m/z, intensity, origin, target matches.
//!
//! A summary report lists all analysed PSMs with links to per-PSM reports.

use std::fmt::Write;
use std::path::Path;

use crate::types::{
    LabelForm, MultiAnnotatedPeak, MultiTargetProvenance, TargetIonMatch,
};

// ---------------------------------------------------------------------------
// Colour constants
// ---------------------------------------------------------------------------

/// Colours assigned to target candidates (cycled if > 10 candidates).
const CANDIDATE_COLORS: &[&str] = &[
    "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b",
    "#e377c2", "#bcbd22", "#17becf", "#1f77b4", "#aec7e8",
];

/// Colour for trap-only peaks.
const TRAP_COLOR: &str = "#1f77b4";

/// Colour for shared (trap + target) peaks.
const SHARED_COLOR: &str = "#9467bd";

/// Colour for unassigned peaks.
const UNASSIGNED_COLOR: &str = "#7f7f7f";

// ---------------------------------------------------------------------------
// Helper: classify peak origin
// ---------------------------------------------------------------------------

/// Origin category for a single peak in the multi-target context.
enum PeakOrigin {
    TrapOnly,
    Shared,
    TargetOnly,
    Unassigned,
}

fn classify_peak(peak: &MultiAnnotatedPeak) -> PeakOrigin {
    let has_trap = peak.trap_ion.is_some();
    let has_target = !peak.target_matches.is_empty();
    match (has_trap, has_target) {
        (true, true) => PeakOrigin::Shared,
        (true, false) => PeakOrigin::TrapOnly,
        (false, true) => PeakOrigin::TargetOnly,
        (false, false) => PeakOrigin::Unassigned,
    }
}

fn origin_label(origin: &PeakOrigin) -> &'static str {
    match origin {
        PeakOrigin::TrapOnly => "TrapOnly",
        PeakOrigin::Shared => "Shared",
        PeakOrigin::TargetOnly => "TargetOnly",
        PeakOrigin::Unassigned => "Unassigned",
    }
}

fn origin_css_class(origin: &PeakOrigin) -> &'static str {
    match origin {
        PeakOrigin::TrapOnly => "origin-trap",
        PeakOrigin::Shared => "origin-shared",
        PeakOrigin::TargetOnly => "origin-target",
        PeakOrigin::Unassigned => "origin-unassigned",
    }
}

fn candidate_color(index: usize) -> &'static str {
    CANDIDATE_COLORS[index % CANDIDATE_COLORS.len()]
}

fn label_form_str(lf: &LabelForm) -> &'static str {
    match lf {
        LabelForm::Light => "Light",
        LabelForm::Heavy { .. } => "Heavy",
    }
}

// ---------------------------------------------------------------------------
// Escape helper
// ---------------------------------------------------------------------------

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

// ---------------------------------------------------------------------------
// CSS
// ---------------------------------------------------------------------------

fn css_block() -> &'static str {
    r#"<style>
body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background: #fafafa; color: #333; }
.header { background: #2c3e50; color: white; padding: 20px; border-radius: 8px; margin-bottom: 20px; }
.header h1 { margin: 0 0 8px 0; font-size: 1.3em; }
.header .meta { font-size: 0.9em; opacity: 0.9; }
.section { background: white; border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin-bottom: 20px; }
.section h2 { margin-top: 0; font-size: 1.1em; border-bottom: 2px solid #3498db; padding-bottom: 6px; }
table { border-collapse: collapse; width: 100%; font-size: 0.85em; }
th, td { padding: 6px 10px; border: 1px solid #ddd; text-align: left; }
th { background: #f0f0f0; font-weight: 600; }
tr:nth-child(even) { background: #fafafa; }
.color-dot { display: inline-block; width: 12px; height: 12px; border-radius: 50%; margin-right: 4px; vertical-align: middle; }
.origin-trap { color: #1f77b4; font-weight: 600; }
.origin-shared { color: #9467bd; font-weight: 600; }
.origin-target { color: #d62728; font-weight: 600; }
.origin-unassigned { color: #7f7f7f; }
.chimeric { background: #fff3cd !important; }
.footer { text-align: center; color: #888; font-size: 0.8em; padding: 12px; }
#mirror-plot { width: 100%; height: 450px; }
</style>"#
}

// ---------------------------------------------------------------------------
// Per-PSM report
// ---------------------------------------------------------------------------

/// Generate a complete, self-contained HTML string for a multi-target
/// provenance report of a single trap PSM.
pub fn generate_multi_provenance_html(prov: &MultiTargetProvenance) -> String {
    let mut html = String::with_capacity(16384);

    // -- doctype + head --
    let _ = write!(
        html,
        r#"<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Provenance: {trap}</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
{css}
</head><body>
"#,
        trap = html_escape(&prov.trap_peptide),
        css = css_block(),
    );

    // -- header --
    let total_peaks = prov.annotated_peaks.len();
    let _ = write!(
        html,
        r#"<div class="header">
<h1>Fragment Ion Provenance Report</h1>
<div class="meta">Trap peptide: <b>{trap}</b> &nbsp;|&nbsp; Scan: <b>{scan}</b> &nbsp;|&nbsp; Candidates: <b>{ncand}</b> &nbsp;|&nbsp; Peaks: <b>{npeaks}</b></div>
</div>
"#,
        trap = html_escape(&prov.trap_peptide),
        scan = prov.scan_number,
        ncand = prov.candidates.len(),
        npeaks = total_peaks,
    );

    // -- candidate table --
    write_candidate_table(&mut html, prov);

    // -- mirror spectrum --
    write_mirror_spectrum(&mut html, prov);

    // -- attribution table --
    write_attribution_table(&mut html, prov);

    // -- footer --
    let _ = write!(
        html,
        r#"<div class="footer">
TrapOnly: {to} &nbsp;|&nbsp; Shared: {sh} &nbsp;|&nbsp; TargetOnly: {tgt} &nbsp;|&nbsp; Unassigned: {ua} &nbsp;|&nbsp; Total: {total}
</div>
</body></html>"#,
        to = prov.trap_only_count,
        sh = prov.shared_count,
        tgt = prov.target_only_count,
        ua = prov.unassigned_count,
        total = total_peaks,
    );

    html
}

/// Write the multi-target provenance report to an HTML file.
pub fn render_multi_provenance_report(
    prov: &MultiTargetProvenance,
    output_path: &Path,
) -> Result<(), std::io::Error> {
    let html = generate_multi_provenance_html(prov);
    std::fs::write(output_path, html)
}

// ---------------------------------------------------------------------------
// Candidate table section
// ---------------------------------------------------------------------------

fn write_candidate_table(html: &mut String, prov: &MultiTargetProvenance) {
    let _ = write!(
        html,
        r#"<div class="section">
<h2>Co-Eluting Target Candidates</h2>
<table>
<tr><th>#</th><th>Peptide</th><th>Protein</th><th>Label</th><th>m/z</th><th>Charge</th><th>RT Range (min)</th><th>Matched Ions</th></tr>
"#,
    );

    for (i, cand) in prov.candidates.iter().enumerate() {
        let matched = count_target_matches_for_candidate(prov, i);
        let color = candidate_color(i);
        let proteins = cand.protein_ids.join("; ");
        let _ = writeln!(
            html,
            r#"<tr><td><span class="color-dot" style="background:{color}"></span>{idx}</td><td>{pep}</td><td>{prot}</td><td>{label}</td><td>{mz:.4}</td><td>{z}+</td><td>{rt0:.2} – {rt1:.2}</td><td>{matched}</td></tr>"#,
            color = color,
            idx = i,
            pep = html_escape(&cand.peptide),
            prot = html_escape(&proteins),
            label = label_form_str(&cand.label_form),
            mz = cand.precursor_mz,
            z = cand.charge,
            rt0 = cand.rt_start,
            rt1 = cand.rt_stop,
            matched = matched,
        );
    }

    let _ = write!(html, "</table>\n</div>\n");
}

fn count_target_matches_for_candidate(prov: &MultiTargetProvenance, candidate_index: usize) -> usize {
    prov.annotated_peaks
        .iter()
        .filter(|p| p.target_matches.iter().any(|m| m.candidate_index == candidate_index))
        .count()
}

// ---------------------------------------------------------------------------
// Mirror spectrum section (Plotly)
// ---------------------------------------------------------------------------

fn write_mirror_spectrum(html: &mut String, prov: &MultiTargetProvenance) {
    let _ = write!(
        html,
        r#"<div class="section">
<h2>Mirror Spectrum</h2>
<div id="mirror-plot"></div>
<script>
"#,
    );

    // Build traces:
    // 1. Trap-matched peaks (positive Y): TrapOnly = blue, Shared = purple
    // 2. Target-matched peaks (negative Y): coloured by candidate index
    // 3. Unassigned (positive, gray, dimmed)

    // -- Trap trace (TrapOnly peaks) --
    write_plotly_trace(
        html,
        "trap_only",
        "Trap Only",
        TRAP_COLOR,
        &collect_peaks_by_origin(prov, |o| matches!(o, PeakOrigin::TrapOnly)),
        true,
    );

    // -- Shared trace (positive Y, purple) --
    write_plotly_trace(
        html,
        "shared",
        "Shared",
        SHARED_COLOR,
        &collect_peaks_by_origin(prov, |o| matches!(o, PeakOrigin::Shared)),
        true,
    );

    // -- Unassigned trace (positive Y, gray) --
    write_plotly_trace(
        html,
        "unassigned",
        "Unassigned",
        UNASSIGNED_COLOR,
        &collect_peaks_by_origin(prov, |o| matches!(o, PeakOrigin::Unassigned)),
        true,
    );

    // -- Target-only traces (negative Y, one per candidate) --
    let target_only_peaks: Vec<&MultiAnnotatedPeak> = prov
        .annotated_peaks
        .iter()
        .filter(|p| matches!(classify_peak(p), PeakOrigin::TargetOnly))
        .collect();

    // Also include shared peaks in negative direction per candidate
    let shared_peaks: Vec<&MultiAnnotatedPeak> = prov
        .annotated_peaks
        .iter()
        .filter(|p| matches!(classify_peak(p), PeakOrigin::Shared))
        .collect();

    for (ci, cand) in prov.candidates.iter().enumerate() {
        let color = candidate_color(ci);
        let mut mzs = Vec::new();
        let mut ints = Vec::new();
        let mut labels = Vec::new();

        // Target-only peaks matching this candidate
        for peak in &target_only_peaks {
            for m in &peak.target_matches {
                if m.candidate_index == ci {
                    mzs.push(peak.mz_observed);
                    ints.push(-peak.intensity); // negative Y
                    labels.push(m.ion_label.clone());
                }
            }
        }

        // Shared peaks matching this candidate (show in negative too)
        for peak in &shared_peaks {
            for m in &peak.target_matches {
                if m.candidate_index == ci {
                    mzs.push(peak.mz_observed);
                    ints.push(-peak.intensity);
                    labels.push(m.ion_label.clone());
                }
            }
        }

        if mzs.is_empty() {
            continue;
        }

        let trace_name = format!("Target: {}", cand.peptide);
        let var_name = format!("target_{ci}");
        write_plotly_trace_raw(html, &var_name, &trace_name, color, &mzs, &ints, &labels);
    }

    // -- Layout + render --
    let _ = write!(
        html,
        r#"
var traces = [trace_trap_only, trace_shared, trace_unassigned"#,
    );

    for ci in 0..prov.candidates.len() {
        let var_name = format!("trace_target_{ci}");
        // Only add if the variable was created
        let _ = write!(html, r#",
  typeof {var} !== 'undefined' ? {var} : null"#, var = var_name);
    }

    let _ = write!(
        html,
        r#"].filter(function(t) {{ return t !== null; }});

var layout = {{
  title: 'Fragment Ion Mirror Spectrum',
  xaxis: {{ title: 'm/z' }},
  yaxis: {{ title: 'Intensity' }},
  barmode: 'overlay',
  hovermode: 'closest',
  showlegend: true,
  legend: {{ x: 1, xanchor: 'right', y: 1 }}
}};

Plotly.newPlot('mirror-plot', traces, layout, {{responsive: true}});
</script>
</div>
"#,
    );
}

/// Collect peaks matching a given origin predicate.
struct PlotPeak {
    mz: f64,
    intensity: f64,
    label: String,
}

fn collect_peaks_by_origin<F>(prov: &MultiTargetProvenance, pred: F) -> Vec<PlotPeak>
where
    F: Fn(PeakOrigin) -> bool,
{
    prov.annotated_peaks
        .iter()
        .filter(|p| pred(classify_peak(p)))
        .map(|p| PlotPeak {
            mz: p.mz_observed,
            intensity: p.intensity,
            label: p.trap_ion.clone().unwrap_or_default(),
        })
        .collect()
}

fn write_plotly_trace(
    html: &mut String,
    var_suffix: &str,
    name: &str,
    color: &str,
    peaks: &[PlotPeak],
    positive: bool,
) {
    let mzs: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
    let ints: Vec<f64> = peaks
        .iter()
        .map(|p| if positive { p.intensity } else { -p.intensity })
        .collect();
    let labels: Vec<String> = peaks.iter().map(|p| p.label.clone()).collect();
    write_plotly_trace_raw(html, var_suffix, name, color, &mzs, &ints, &labels);
}

fn write_plotly_trace_raw(
    html: &mut String,
    var_suffix: &str,
    name: &str,
    color: &str,
    mzs: &[f64],
    ints: &[f64],
    labels: &[String],
) {
    let mz_json: Vec<String> = mzs.iter().map(|v| format!("{v}")).collect();
    let int_json: Vec<String> = ints.iter().map(|v| format!("{v}")).collect();
    let label_json: Vec<String> = labels.iter().map(|l| format!("\"{}\"", html_escape(l))).collect();

    let _ = write!(
        html,
        r#"var trace_{var} = {{
  x: [{mzs}],
  y: [{ints}],
  text: [{labels}],
  name: '{name}',
  type: 'bar',
  marker: {{ color: '{color}' }},
  hovertemplate: '%{{text}}<br>m/z: %{{x:.4f}}<br>Intensity: %{{y:.0f}}<extra></extra>'
}};
"#,
        var = var_suffix,
        mzs = mz_json.join(", "),
        ints = int_json.join(", "),
        labels = label_json.join(", "),
        name = html_escape(name),
        color = color,
    );
}

// ---------------------------------------------------------------------------
// Attribution table section
// ---------------------------------------------------------------------------

fn write_attribution_table(html: &mut String, prov: &MultiTargetProvenance) {
    let _ = write!(
        html,
        r#"<div class="section">
<h2>Fragment Ion Attribution</h2>
<table>
<tr><th>m/z</th><th>Intensity</th><th>Trap Ion</th><th>Origin</th><th>Target Matches</th></tr>
"#,
    );

    for peak in &prov.annotated_peaks {
        let origin = classify_peak(peak);
        let origin_str = origin_label(&origin);
        let css = origin_css_class(&origin);

        let trap_label = peak.trap_ion.as_deref().unwrap_or("–");

        let target_desc = if peak.target_matches.is_empty() {
            "–".to_string()
        } else {
            peak.target_matches
                .iter()
                .map(|m| format_target_match(m, prov))
                .collect::<Vec<_>>()
                .join("; ")
        };

        let _ = writeln!(
            html,
            r#"<tr><td>{mz:.3}</td><td>{int:.0}</td><td>{trap}</td><td class="{css}">{origin}</td><td>{target}</td></tr>"#,
            mz = peak.mz_observed,
            int = peak.intensity,
            trap = html_escape(trap_label),
            css = css,
            origin = origin_str,
            target = html_escape(&target_desc),
        );
    }

    let _ = write!(html, "</table>\n</div>\n");
}

fn format_target_match(m: &TargetIonMatch, prov: &MultiTargetProvenance) -> String {
    let pep = prov
        .candidates
        .get(m.candidate_index)
        .map(|c| c.peptide.as_str())
        .unwrap_or("?");
    format!(
        "[C{idx}:{pep}] {ion} (Δ{ppm:+.1}ppm)",
        idx = m.candidate_index,
        pep = pep,
        ion = m.ion_label,
        ppm = m.delta_ppm,
    )
}

// ---------------------------------------------------------------------------
// Summary report
// ---------------------------------------------------------------------------

/// Generate an HTML summary report listing all multi-target provenance results.
///
/// Each row links to a per-PSM report (expected at `{scan_number}.html`).
/// Rows where shared fraction exceeds 30% are highlighted as chimeric.
pub fn generate_provenance_summary_html(results: &[MultiTargetProvenance]) -> String {
    let mut html = String::with_capacity(8192);

    let _ = write!(
        html,
        r#"<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Provenance Summary</title>
{css}
</head><body>
<div class="header">
<h1>Provenance Summary</h1>
<div class="meta">Total PSMs analysed: <b>{n}</b></div>
</div>
<div class="section">
<h2>All Analysed Trap PSMs</h2>
<table>
<tr><th>Peptide</th><th>Scan</th><th>#Candidates</th><th>TrapOnly</th><th>Shared</th><th>TargetOnly</th><th>Unassigned</th><th>Report</th></tr>
"#,
        css = css_block(),
        n = results.len(),
    );

    for prov in results {
        let total = prov.trap_only_count + prov.shared_count + prov.target_only_count + prov.unassigned_count;
        let shared_frac = if total > 0 {
            f64::from(prov.shared_count) / f64::from(total)
        } else {
            0.0
        };
        let row_class = if shared_frac > 0.30 { " class=\"chimeric\"" } else { "" };

        let scan = prov.scan_number;
        let _ = writeln!(
            html,
            r#"<tr{cls}><td>{pep}</td><td>{scan}</td><td>{ncand}</td><td>{to}</td><td>{sh}</td><td>{tgt}</td><td>{ua}</td><td><a href="{scan}.html">view</a></td></tr>"#,
            cls = row_class,
            pep = html_escape(&prov.trap_peptide),
            scan = scan,
            ncand = prov.candidates.len(),
            to = prov.trap_only_count,
            sh = prov.shared_count,
            tgt = prov.target_only_count,
            ua = prov.unassigned_count,
        );
    }

    let _ = write!(html, "</table>\n</div>\n</body></html>");

    html
}

/// Write the summary report to an HTML file.
pub fn render_provenance_summary(
    results: &[MultiTargetProvenance],
    output_path: &Path,
) -> Result<(), std::io::Error> {
    let html = generate_provenance_summary_html(results);
    std::fs::write(output_path, html)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::CoElutingCandidate;

    fn make_test_provenance() -> MultiTargetProvenance {
        MultiTargetProvenance {
            trap_peptide: "STTTGHLIYK".to_string(),
            scan_number: 12345,
            candidates: vec![CoElutingCandidate {
                peptide: "STTSGHLVYK".to_string(),
                protein_ids: vec!["sp|P12345|EF1A_HUMAN".to_string()],
                precursor_mz: 548.12,
                charge: 2,
                rt_start: 34.5,
                rt_stop: 35.8,
                label_form: LabelForm::Light,
                modifications: vec![],
            }],
            annotated_peaks: vec![
                MultiAnnotatedPeak {
                    mz_observed: 285.155,
                    intensity: 45230.0,
                    trap_ion: Some("b3+1".to_string()),
                    target_matches: vec![TargetIonMatch {
                        candidate_index: 0,
                        ion_label: "b3+1".to_string(),
                        delta_ppm: -2.1,
                    }],
                },
                MultiAnnotatedPeak {
                    mz_observed: 386.203,
                    intensity: 72100.0,
                    trap_ion: Some("b4+1".to_string()),
                    target_matches: vec![],
                },
                MultiAnnotatedPeak {
                    mz_observed: 512.334,
                    intensity: 8200.0,
                    trap_ion: None,
                    target_matches: vec![],
                },
            ],
            trap_only_count: 1,
            target_only_count: 0,
            shared_count: 1,
            unassigned_count: 1,
        }
    }

    #[test]
    fn test_generate_html_contains_sections() {
        let prov = make_test_provenance();
        let html = generate_multi_provenance_html(&prov);
        assert!(html.contains("Fragment Ion Provenance Report"));
        assert!(html.contains("STTTGHLIYK"));
        assert!(html.contains("STTSGHLVYK"));
        assert!(html.contains("285.155"));
        assert!(html.contains("plotly"));
        assert!(html.contains("b3+1"));
        assert!(html.contains("TrapOnly"));
        assert!(html.contains("Shared"));
    }

    #[test]
    fn test_write_html_file() {
        let prov = make_test_provenance();
        let dir = std::env::temp_dir().join("test_multi_report");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_report.html");
        render_multi_provenance_report(&prov, &path).unwrap();
        assert!(path.exists());
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<!DOCTYPE html>"));
        assert!(content.contains("STTTGHLIYK"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_summary_report() {
        let provs = vec![make_test_provenance()];
        let html = generate_provenance_summary_html(&provs);
        assert!(html.contains("Provenance Summary"));
        assert!(html.contains("STTTGHLIYK"));
        assert!(html.contains("<!DOCTYPE html>"));
    }

    #[test]
    fn test_classify_peak_origins() {
        let shared = MultiAnnotatedPeak {
            mz_observed: 100.0,
            intensity: 1000.0,
            trap_ion: Some("b2+1".to_string()),
            target_matches: vec![TargetIonMatch {
                candidate_index: 0,
                ion_label: "b2+1".to_string(),
                delta_ppm: 1.0,
            }],
        };
        assert!(matches!(classify_peak(&shared), PeakOrigin::Shared));

        let trap_only = MultiAnnotatedPeak {
            mz_observed: 200.0,
            intensity: 500.0,
            trap_ion: Some("y3+1".to_string()),
            target_matches: vec![],
        };
        assert!(matches!(classify_peak(&trap_only), PeakOrigin::TrapOnly));

        let target_only = MultiAnnotatedPeak {
            mz_observed: 300.0,
            intensity: 800.0,
            trap_ion: None,
            target_matches: vec![TargetIonMatch {
                candidate_index: 0,
                ion_label: "y4+1".to_string(),
                delta_ppm: -0.5,
            }],
        };
        assert!(matches!(classify_peak(&target_only), PeakOrigin::TargetOnly));

        let unassigned = MultiAnnotatedPeak {
            mz_observed: 400.0,
            intensity: 200.0,
            trap_ion: None,
            target_matches: vec![],
        };
        assert!(matches!(classify_peak(&unassigned), PeakOrigin::Unassigned));
    }

    #[test]
    fn test_html_escape() {
        assert_eq!(html_escape("<b>"), "&lt;b&gt;");
        assert_eq!(html_escape("a&b"), "a&amp;b");
        assert_eq!(html_escape("\"q\""), "&quot;q&quot;");
    }

    #[test]
    fn test_candidate_color_wraps() {
        // Should wrap around at 10
        assert_eq!(candidate_color(0), CANDIDATE_COLORS[0]);
        assert_eq!(candidate_color(10), CANDIDATE_COLORS[0]);
        assert_eq!(candidate_color(3), CANDIDATE_COLORS[3]);
    }

    #[test]
    fn test_summary_chimeric_highlight() {
        // Create provenance with >30% shared
        let prov = MultiTargetProvenance {
            trap_peptide: "CHIMERIC".to_string(),
            scan_number: 999,
            candidates: vec![],
            annotated_peaks: vec![],
            trap_only_count: 1,
            target_only_count: 0,
            shared_count: 5, // 5/10 = 50% → chimeric
            unassigned_count: 4,
        };
        let html = generate_provenance_summary_html(&[prov]);
        assert!(html.contains("chimeric"));
    }

    #[test]
    fn test_empty_candidates() {
        let prov = MultiTargetProvenance {
            trap_peptide: "EMPTYK".to_string(),
            scan_number: 1,
            candidates: vec![],
            annotated_peaks: vec![],
            trap_only_count: 0,
            target_only_count: 0,
            shared_count: 0,
            unassigned_count: 0,
        };
        let html = generate_multi_provenance_html(&prov);
        assert!(html.contains("EMPTYK"));
        assert!(html.contains("Co-Eluting Target Candidates"));
    }

    #[test]
    fn test_render_summary_file() {
        let provs = vec![make_test_provenance()];
        let dir = std::env::temp_dir().join("test_summary_report");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("summary.html");
        render_provenance_summary(&provs, &path).unwrap();
        assert!(path.exists());
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("Provenance Summary"));
        std::fs::remove_dir_all(&dir).ok();
    }
}
