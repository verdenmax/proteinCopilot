//! Per-PSM HTML Report Renderer for Multi-Target Fragment Provenance.
//!
//! Generates self-contained HTML reports that visualise which fragment ions
//! of a trap PSM originate from co-eluting target peptides.  Each report
//! contains:
//!
//! 1. **Trap PSM info** — peptide, precursor m/z (light+heavy), charge, scan, file.
//! 2. **Candidate table** — co-eluting target peptides with precursor m/z, file, etc.
//! 3. **Light mirror spectrum** — normalized, matched peaks bold (trap ↑ vs light targets ↓).
//! 4. **Heavy mirror spectrum** — normalized, matched peaks bold (trap ↑ vs heavy targets ↓).
//! 5. **Trap ion attribution** — only trap fragment ions, showing target matches.
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

/// Check if a candidate is a Light form.
fn is_light(lf: &LabelForm) -> bool {
    matches!(lf, LabelForm::Light)
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
.header .meta { font-size: 0.9em; opacity: 0.9; line-height: 1.8; }
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
.mirror-plot { width: 100%; height: 420px; }
.info-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 4px 24px; }
.info-grid .label { color: rgba(255,255,255,0.7); font-size: 0.85em; }
.info-grid .value { font-weight: 600; }
</style>"#
}

// ---------------------------------------------------------------------------
// Per-PSM report
// ---------------------------------------------------------------------------

/// Generate a complete, self-contained HTML string for a multi-target
/// provenance report of a single trap PSM.
pub fn generate_multi_provenance_html(prov: &MultiTargetProvenance) -> String {
    let mut html = String::with_capacity(32768);

    // -- doctype + head --
    let _ = write!(
        html,
        r#"<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Provenance: {trap} (scan {scan})</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
{css}
</head><body>
"#,
        trap = html_escape(&prov.trap_peptide),
        scan = prov.scan_number,
        css = css_block(),
    );

    // -- header with extended info --
    write_header(&mut html, prov);

    // -- candidate table --
    write_candidate_table(&mut html, prov);

    // -- mirror spectra: light and heavy --
    let has_light_candidates = prov.candidates.iter().any(|c| is_light(&c.label_form));
    let has_heavy_candidates = prov.candidates.iter().any(|c| !is_light(&c.label_form));

    if has_light_candidates {
        write_mirror_spectrum(&mut html, prov, MirrorKind::Light);
    }
    if has_heavy_candidates {
        write_mirror_spectrum(&mut html, prov, MirrorKind::Heavy);
    }
    // Fallback: if no candidates of either kind, show combined mirror
    if !has_light_candidates && !has_heavy_candidates {
        write_mirror_spectrum(&mut html, prov, MirrorKind::Light);
    }

    // -- attribution table (trap ions only) --
    write_attribution_table(&mut html, prov);

    // -- footer --
    let total_peaks = prov.annotated_peaks.len();
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
// Header section (enhanced)
// ---------------------------------------------------------------------------

fn write_header(html: &mut String, prov: &MultiTargetProvenance) {
    let heavy_mz_str = match prov.trap_precursor_mz_heavy {
        Some(mz) => format!("{mz:.4}"),
        None => "N/A".to_string(),
    };

    let _ = write!(
        html,
        r#"<div class="header">
<h1>Fragment Ion Provenance Report</h1>
<div class="meta">
<div class="info-grid">
<div><span class="label">Trap Peptide:</span> <span class="value">{trap}</span></div>
<div><span class="label">Spectrum File:</span> <span class="value">{file}</span></div>
<div><span class="label">Scan:</span> <span class="value">{scan}</span></div>
<div><span class="label">Charge:</span> <span class="value">{charge}+</span></div>
<div><span class="label">Precursor m/z (Light):</span> <span class="value">{mz_light:.4}</span></div>
<div><span class="label">Precursor m/z (Heavy):</span> <span class="value">{mz_heavy}</span></div>
<div><span class="label">Candidates:</span> <span class="value">{ncand}</span></div>
<div><span class="label">Peaks:</span> <span class="value">{npeaks}</span></div>
</div>
</div>
</div>
"#,
        trap = html_escape(&prov.trap_peptide),
        file = html_escape(&prov.spectrum_file),
        scan = prov.scan_number,
        charge = prov.trap_charge,
        mz_light = prov.trap_precursor_mz,
        mz_heavy = heavy_mz_str,
        ncand = prov.candidates.len(),
        npeaks = prov.annotated_peaks.len(),
    );
}

// ---------------------------------------------------------------------------
// Candidate table section (enhanced)
// ---------------------------------------------------------------------------

fn write_candidate_table(html: &mut String, prov: &MultiTargetProvenance) {
    let _ = write!(
        html,
        r#"<div class="section">
<h2>Co-Eluting Target Candidates</h2>
<table>
<tr><th>#</th><th>Peptide</th><th>Protein</th><th>Label</th><th>Precursor m/z</th><th>Charge</th><th>RT Range (min)</th><th>Spectrum File</th><th>Matched Ions</th></tr>
"#,
    );

    for (i, cand) in prov.candidates.iter().enumerate() {
        let matched = count_target_matches_for_candidate(prov, i);
        let color = candidate_color(i);
        let proteins = cand.protein_ids.join("; ");
        // For Heavy candidates, show both light and heavy m/z
        let mz_display = match &cand.label_form {
            LabelForm::Heavy {
                precursor_mz_heavy, ..
            } => format!("{:.4} (H: {:.4})", cand.precursor_mz, precursor_mz_heavy),
            LabelForm::Light => format!("{:.4}", cand.precursor_mz),
        };
        let _ = writeln!(
            html,
            r#"<tr><td><span class="color-dot" style="background:{color}"></span>{idx}</td><td>{pep}</td><td>{prot}</td><td>{label}</td><td>{mz}</td><td>{z}+</td><td>{rt0:.2} – {rt1:.2}</td><td>{file}</td><td>{matched}</td></tr>"#,
            color = color,
            idx = i,
            pep = html_escape(&cand.peptide),
            prot = html_escape(&proteins),
            label = label_form_str(&cand.label_form),
            mz = mz_display,
            z = cand.charge,
            rt0 = cand.rt_start,
            rt1 = cand.rt_stop,
            file = html_escape(&prov.spectrum_file),
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
// Mirror spectrum section (Plotly) — Light and Heavy variants
// ---------------------------------------------------------------------------

/// Which mirror to render.
#[derive(Clone, Copy, PartialEq, Eq)]
enum MirrorKind {
    Light,
    Heavy,
}

fn write_mirror_spectrum(html: &mut String, prov: &MultiTargetProvenance, kind: MirrorKind) {
    let kind_label = match kind {
        MirrorKind::Light => "Light",
        MirrorKind::Heavy => "Heavy",
    };
    let plot_id = match kind {
        MirrorKind::Light => "mirror-plot-light",
        MirrorKind::Heavy => "mirror-plot-heavy",
    };

    // Collect candidate indices of the requested kind.
    let candidate_indices: Vec<usize> = prov
        .candidates
        .iter()
        .enumerate()
        .filter(|(_, c)| match kind {
            MirrorKind::Light => is_light(&c.label_form),
            MirrorKind::Heavy => !is_light(&c.label_form),
        })
        .map(|(i, _)| i)
        .collect();

    // Compute max intensity for normalization.
    let max_intensity = prov
        .annotated_peaks
        .iter()
        .map(|p| p.intensity)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let _ = write!(
        html,
        r#"<div class="section">
<h2>Mirror Spectrum — {kind} Targets (Scan {scan})</h2>
<div id="{plot_id}" class="mirror-plot"></div>
<script>
(function() {{
var maxI = {max_intensity};
function norm(v) {{ return v / maxI * 100.0; }}
"#,
        kind = kind_label,
        scan = prov.scan_number,
        plot_id = plot_id,
        max_intensity = max_intensity,
    );

    // Determine which peaks have target matches in the current kind.
    let matching_target_indices: std::collections::HashSet<usize> =
        candidate_indices.iter().copied().collect();

    // Classify peaks for the top half (trap).
    // Peaks matched by this kind's targets = shared; only trap = trap-only.
    let mut trap_only_mz = Vec::new();
    let mut trap_only_int = Vec::new();
    let mut trap_only_labels = Vec::new();
    let mut shared_mz = Vec::new();
    let mut shared_int = Vec::new();
    let mut shared_labels = Vec::new();
    let mut unassigned_mz = Vec::new();
    let mut unassigned_int = Vec::new();

    for peak in &prov.annotated_peaks {
        let norm_int = peak.intensity / max_intensity * 100.0;
        let has_trap = peak.trap_ion.is_some();
        let has_kind_target = peak
            .target_matches
            .iter()
            .any(|m| matching_target_indices.contains(&m.candidate_index));

        if has_trap && has_kind_target {
            shared_mz.push(peak.mz_observed);
            shared_int.push(norm_int);
            shared_labels.push(peak.trap_ion.clone().unwrap_or_default());
        } else if has_trap {
            trap_only_mz.push(peak.mz_observed);
            trap_only_int.push(norm_int);
            trap_only_labels.push(peak.trap_ion.clone().unwrap_or_default());
        } else if !has_trap && !peak.target_matches.is_empty() && !has_kind_target {
            // Target-only for the OTHER kind — show as unassigned in this mirror
            unassigned_mz.push(peak.mz_observed);
            unassigned_int.push(norm_int);
        } else if !has_trap && peak.target_matches.is_empty() {
            unassigned_mz.push(peak.mz_observed);
            unassigned_int.push(norm_int);
        }
        // target-only for THIS kind shown only in bottom half
    }

    // Top traces — normalized, bold (wider bar) for shared
    write_normalized_trace(
        html,
        "trap_only",
        "Trap Only",
        TRAP_COLOR,
        &trap_only_mz,
        &trap_only_int,
        &trap_only_labels,
        false,
    );
    write_normalized_trace(
        html,
        "shared_up",
        "Shared (Trap↑)",
        SHARED_COLOR,
        &shared_mz,
        &shared_int,
        &shared_labels,
        true, // bold
    );
    write_normalized_trace(
        html,
        "unassigned",
        "Unassigned",
        UNASSIGNED_COLOR,
        &unassigned_mz,
        &unassigned_int,
        &Vec::new(),
        false,
    );

    // Bottom traces — one per candidate of this kind
    let mut target_trace_vars = Vec::new();
    for &ci in &candidate_indices {
        let cand = &prov.candidates[ci];
        let color = candidate_color(ci);
        let mut mzs = Vec::new();
        let mut ints = Vec::new();
        let mut labels = Vec::new();
        let mut bold_flags = Vec::new();

        for peak in &prov.annotated_peaks {
            for m in &peak.target_matches {
                if m.candidate_index == ci {
                    let norm_int = peak.intensity / max_intensity * 100.0;
                    mzs.push(peak.mz_observed);
                    ints.push(-norm_int); // negative for bottom
                    labels.push(m.ion_label.clone());
                    bold_flags.push(peak.trap_ion.is_some()); // bold if shared
                }
            }
        }

        if mzs.is_empty() {
            continue;
        }

        let var_name = format!("target_{ci}");
        let trace_name = format!("[C{ci}] {pep} ({label})",
            ci = ci,
            pep = cand.peptide,
            label = label_form_str(&cand.label_form),
        );

        // Write trace with per-bar bold outline for shared peaks.
        write_target_trace_with_bold(html, &var_name, &trace_name, color, &mzs, &ints, &labels, &bold_flags);
        target_trace_vars.push(format!("trace_{var_name}"));
    }

    // Assemble and render
    let _ = write!(
        html,
        r#"
var traces = [trace_trap_only, trace_shared_up, trace_unassigned"#,
    );
    for v in &target_trace_vars {
        let _ = write!(html, ", {v}");
    }

    let _ = write!(
        html,
        r#"];

var layout = {{
  title: '{kind} Mirror — Scan {scan}',
  xaxis: {{ title: 'm/z' }},
  yaxis: {{ title: 'Relative Intensity (%)' }},
  barmode: 'overlay',
  hovermode: 'closest',
  showlegend: true,
  legend: {{ x: 1, xanchor: 'right', y: 1 }},
  shapes: [{{ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 0, y1: 0, line: {{ color: '#333', width: 1 }} }}]
}};

Plotly.newPlot('{plot_id}', traces, layout, {{responsive: true}});
}})();
</script>
</div>
"#,
        kind = kind_label,
        scan = prov.scan_number,
        plot_id = plot_id,
    );
}

/// Write a normalized bar trace for the top (trap) half.
fn write_normalized_trace(
    html: &mut String,
    var_suffix: &str,
    name: &str,
    color: &str,
    mzs: &[f64],
    ints: &[f64],
    labels: &[String],
    bold: bool,
) {
    let mz_json: Vec<String> = mzs.iter().map(|v| format!("{v}")).collect();
    let int_json: Vec<String> = ints.iter().map(|v| format!("{v}")).collect();
    let label_json: Vec<String> = if labels.is_empty() {
        ints.iter().map(|_| "\"\"".to_string()).collect()
    } else {
        labels.iter().map(|l| format!("\"{}\"", html_escape(l))).collect()
    };

    let line_width = if bold { 2.5 } else { 0.0 };
    let line_color = if bold { "#2c3e50" } else { color };

    let _ = write!(
        html,
        r#"var trace_{var} = {{
  x: [{mzs}],
  y: [{ints}],
  text: [{labels}],
  name: '{name}',
  type: 'bar',
  marker: {{ color: '{color}', line: {{ width: {lw}, color: '{lc}' }} }},
  hovertemplate: '%{{text}}<br>m/z: %{{x:.4f}}<br>Rel.Int: %{{y:.1f}}%<extra></extra>'
}};
"#,
        var = var_suffix,
        mzs = mz_json.join(", "),
        ints = int_json.join(", "),
        labels = label_json.join(", "),
        name = html_escape(name),
        color = color,
        lw = line_width,
        lc = line_color,
    );
}

/// Write a target trace with per-bar bold outline (for shared peaks).
fn write_target_trace_with_bold(
    html: &mut String,
    var_suffix: &str,
    name: &str,
    color: &str,
    mzs: &[f64],
    ints: &[f64],
    labels: &[String],
    bold_flags: &[bool],
) {
    let mz_json: Vec<String> = mzs.iter().map(|v| format!("{v}")).collect();
    let int_json: Vec<String> = ints.iter().map(|v| format!("{v}")).collect();
    let label_json: Vec<String> = labels.iter().map(|l| format!("\"{}\"", html_escape(l))).collect();
    let line_widths: Vec<String> = bold_flags.iter().map(|&b| if b { "2.5".to_string() } else { "0".to_string() }).collect();

    let _ = write!(
        html,
        r#"var trace_{var} = {{
  x: [{mzs}],
  y: [{ints}],
  text: [{labels}],
  name: '{name}',
  type: 'bar',
  marker: {{ color: '{color}', line: {{ width: [{lws}], color: '#2c3e50' }} }},
  hovertemplate: '%{{text}}<br>m/z: %{{x:.4f}}<br>Rel.Int: %{{y:.1f}}%<extra></extra>'
}};
"#,
        var = var_suffix,
        mzs = mz_json.join(", "),
        ints = int_json.join(", "),
        labels = label_json.join(", "),
        name = html_escape(name),
        color = color,
        lws = line_widths.join(", "),
    );
}

// ---------------------------------------------------------------------------
// Attribution table — trap ions only
// ---------------------------------------------------------------------------

fn write_attribution_table(html: &mut String, prov: &MultiTargetProvenance) {
    let _ = write!(
        html,
        r#"<div class="section">
<h2>Trap Fragment Ion Attribution</h2>
<p style="color:#666;font-size:0.85em;margin-top:0">Only peaks matching trap fragment ions are shown. Target matches indicate potential chimeric contamination.</p>
<table>
<tr><th>m/z</th><th>Rel. Intensity (%)</th><th>Trap Ion</th><th>Target Matches</th></tr>
"#,
    );

    let max_intensity = prov
        .annotated_peaks
        .iter()
        .map(|p| p.intensity)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    // Only show peaks that match a trap fragment ion.
    for peak in &prov.annotated_peaks {
        let trap_label = match &peak.trap_ion {
            Some(lbl) => lbl,
            None => continue, // skip non-trap peaks
        };

        let rel_int = peak.intensity / max_intensity * 100.0;
        let origin = classify_peak(peak);
        let css = origin_css_class(&origin);

        let target_desc = if peak.target_matches.is_empty() {
            "— (trap unique)".to_string()
        } else {
            peak.target_matches
                .iter()
                .map(|m| format_target_match(m, prov))
                .collect::<Vec<_>>()
                .join("<br>")
        };

        let _ = writeln!(
            html,
            r#"<tr><td>{mz:.4}</td><td>{rel:.1}</td><td class="{css}"><b>{trap}</b></td><td>{target}</td></tr>"#,
            mz = peak.mz_observed,
            rel = rel_int,
            css = css,
            trap = html_escape(trap_label),
            target = target_desc,
        );
    }

    let _ = write!(html, "</table>\n</div>\n");
}

fn format_target_match(m: &TargetIonMatch, prov: &MultiTargetProvenance) -> String {
    let cand = prov.candidates.get(m.candidate_index);
    let pep = cand.map(|c| c.peptide.as_str()).unwrap_or("?");
    let label = cand
        .map(|c| label_form_str(&c.label_form))
        .unwrap_or("?");
    format!(
        "<span class=\"color-dot\" style=\"background:{color}\"></span>[C{idx}] {pep} ({label}) → <b>{ion}</b> (Δ{ppm:+.1} ppm)",
        color = candidate_color(m.candidate_index),
        idx = m.candidate_index,
        pep = html_escape(pep),
        label = label,
        ion = html_escape(&m.ion_label),
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
<tr><th>Peptide</th><th>File</th><th>Scan</th><th>m/z (L)</th><th>Charge</th><th>#Cand</th><th>TrapOnly</th><th>Shared</th><th>TargetOnly</th><th>Unassigned</th><th>Report</th></tr>
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
        let report_filename = format!(
            "provenance/{}_{}_scan{}.html",
            prov.spectrum_file, prov.trap_peptide, scan
        );
        let _ = writeln!(
            html,
            r#"<tr{cls}><td>{pep}</td><td>{file}</td><td>{scan}</td><td>{mz:.4}</td><td>{z}+</td><td>{ncand}</td><td>{to}</td><td>{sh}</td><td>{tgt}</td><td>{ua}</td><td><a href="{report}">view</a></td></tr>"#,
            cls = row_class,
            pep = html_escape(&prov.trap_peptide),
            file = html_escape(&prov.spectrum_file),
            scan = scan,
            mz = prov.trap_precursor_mz,
            z = prov.trap_charge,
            ncand = prov.candidates.len(),
            to = prov.trap_only_count,
            sh = prov.shared_count,
            tgt = prov.target_only_count,
            ua = prov.unassigned_count,
            report = report_filename,
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
            trap_precursor_mz: 547.789,
            trap_precursor_mz_heavy: Some(556.803),
            trap_charge: 2,
            spectrum_file: "550_600_2Da_Rep1".to_string(),
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
        assert!(html.contains("Trap Fragment Ion Attribution"));
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
            trap_precursor_mz: 500.0,
            trap_precursor_mz_heavy: None,
            trap_charge: 2,
            spectrum_file: "test_run".to_string(),
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
            trap_precursor_mz: 300.0,
            trap_precursor_mz_heavy: None,
            trap_charge: 2,
            spectrum_file: "test_run".to_string(),
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
