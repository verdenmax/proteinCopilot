# Unified Annotation + XIC HTML with Interactive SILAC Controls

**Date**: 2026-04-11
**Status**: Approved
**Author**: ProteinCopilot Team

---

## 1. Problem Statement

Currently, spectrum annotation and XIC chromatograms are generated as **separate HTML files**. Users must:
1. Open two files side by side to correlate annotation with XIC
2. Manually specify SILAC `label_type` in the MCP tool call to see heavy traces
3. Re-run the tool with different parameters to try different SILAC configurations

This design unifies both views into a single HTML and adds client-side interactive SILAC controls.

## 2. Goals

- **G1**: Single HTML file with annotation (top) + XIC (bottom) in vertical layout
- **G2**: Interactive SILAC preset switching (None / Standard K8R10 / Medium K4R6 / Custom)
- **G3**: Client-side XIC recomputation from embedded raw scan data — no backend round-trip
- **G4**: Light/Heavy display modes: Overlay, Light-only, Heavy-only, Split
- **G5**: Backward compatible — existing `annotation.html` and `xic.html` templates unchanged

## 3. Non-Goals

- Server-side interactive SILAC (WebSocket / API calls from HTML)
- Multi-file XIC comparison (e.g., two different runs side by side)
- Quantification ratio calculation (L/H ratios)
- Plotly-based spectrum annotation (keep existing SVG implementation)

## 4. Layout

Vertical stack, single-page, no tabs:

```
┌──────────────────────────────────────┐
│ Info Panel (shared metadata)         │
│ Scan: 1234  Charge: 2+  RT: 2.0 min │
│ Peptide: PEPTIDEK  Score: 0.750      │
├──────────────────────────────────────┤
│ Fragment Ion Coverage (SVG brackets) │
│     y₇ y₆ y₅ y₄ y₃ y₂ y₁          │
│  P · E · P · T · I · D · E · K      │
│     b₁ b₂ b₃ b₄ b₅ b₆ b₇          │
├──────────────────────────────────────┤
│ Spectrum Annotation (SVG peaks)      │
│  ║  ║    ║║  ║  ║║║  ║   ║  ║       │
│  └──┴────┴┴──┴──┴┴┴──┴───┴──┘       │
│            m/z →                     │
├──────────────────────────────────────┤
│ ⚙️ SILAC Controls                    │
│ Preset: [Standard SILAC ▼]          │
│ Display: [Overlay ▼]                │
│ Badge: K+8.014 R+10.008             │
├──────────────────────────────────────┤
│ MS1 Precursor XIC (Plotly)           │
│  ── Light  ╌╌ Heavy  ┆ Target RT    │
├──────────────────────────────────────┤
│ MS2 Fragment Ion XIC (Plotly)        │
│  6 top ions × light/heavy traces    │
├──────────────────────────────────────┤
│ Parameters footer                    │
└──────────────────────────────────────┘
```

The Fragment Ion Coverage panel uses the existing SVG bracket notation:
- **b-ions** (below): vertical stub down `|` + horizontal left `_` (N-terminal)
- **y-ions** (above): vertical stub up `|` + horizontal right `¯` (C-terminal)
- Matched ions are colored (b: red, y: blue); unmatched are gray

## 5. Client-Side Interactive SILAC — Core Design

### 5.1 Approach

Embed raw MS1/MS2 peak arrays into the HTML alongside pre-computed XIC traces. JavaScript performs SILAC m/z calculation and peak extraction client-side when the user changes SILAC presets.

### 5.2 Embedded Data Structure

The HTML contains a single JSON blob with all required data:

```rust
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
    pub ms1_scans: Vec<RawScan>,
    pub ms2_scans: Vec<RawScan>,
}

/// A single raw scan's peak list.
#[derive(Debug, Clone, Serialize)]
pub struct RawScan {
    pub scan_number: u32,
    pub retention_time_sec: f64,
    /// m/z values (sorted ascending).
    pub mz_array: Vec<f64>,
    /// Intensity values (parallel to mz_array).
    pub intensity_array: Vec<f64>,
}

/// Metadata for one fragment ion enabling client-side SILAC calculation.
#[derive(Debug, Clone, Serialize)]
pub struct IonMetadataEntry {
    pub label: String,
    pub ion_type: IonType,
    pub ion_number: u32,
    pub charge: u32,
    pub light_mz: f64,
    /// Count of K residues in this fragment.
    pub k_count: u32,
    /// Count of R residues in this fragment.
    pub r_count: u32,
}

/// Peptide-level SILAC info.
#[derive(Debug, Clone, Serialize)]
pub struct PeptideInfo {
    pub sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub total_k: u32,
    pub total_r: u32,
}
```

### 5.3 Data Volume

| Component | Estimate | JSON Size |
|-----------|----------|-----------|
| MS1 raw peaks (±20 Da window around precursor) | ~10 scans × ~200 peaks | ~30 KB |
| MS2 raw peaks (full spectra) | ~10 scans × ~2000 peaks | ~300 KB |
| Pre-computed XIC traces | existing data | ~20 KB |
| Ion metadata | ~12 ions | ~2 KB |
| **Total overhead** | | **~350 KB** |

MS1 scans are trimmed to a narrow m/z window around the precursor (±20 Da) to keep data compact. This covers any realistic SILAC delta. MS2 scans are kept in full since fragment ions span the entire m/z range.

### 5.4 JavaScript Recomputation Algorithm

```javascript
function recomputeSilac(kDelta, rDelta) {
    const data = window.__UNIFIED_DATA__;
    const info = data.peptide_info;
    const tolerance_ppm = 20; // match backend default

    // 1. Heavy precursor m/z
    const precDelta = info.total_k * kDelta + info.total_r * rDelta;
    const heavyPrecMz = info.precursor_mz + precDelta / Math.abs(info.charge);

    // 2. Heavy fragment m/z for each ion
    const heavyIons = data.ion_metadata.map(ion => ({
        ...ion,
        heavy_mz: ion.light_mz + (ion.k_count * kDelta + ion.r_count * rDelta) / ion.charge
    }));

    // 3. Extract intensities from raw scans (binary search)
    function extractIntensity(targetMz, mzArr, intArr) {
        const tolDa = targetMz * tolerance_ppm * 1e-6;
        let lo = 0, hi = mzArr.length;
        while (lo < hi) { let mid = (lo + hi) >> 1; mzArr[mid] < targetMz - tolDa * 1.5 ? lo = mid + 1 : hi = mid; }
        let best = 0;
        for (let i = lo; i < mzArr.length && mzArr[i] <= targetMz + tolDa * 1.5; i++) {
            if (Math.abs(mzArr[i] - targetMz) / targetMz * 1e6 <= tolerance_ppm) {
                best = Math.max(best, intArr[i]);
            }
        }
        return best;
    }

    // 4. Build heavy XIC traces
    const ms1HeavyTrace = data.raw_scans.ms1_scans.map(scan => ({
        rt: scan.retention_time_sec,
        scan: scan.scan_number,
        intensity: extractIntensity(heavyPrecMz, scan.mz_array, scan.intensity_array)
    }));

    const ms2HeavyTraces = heavyIons.map(ion =>
        data.raw_scans.ms2_scans.map(scan => ({
            rt: scan.retention_time_sec,
            scan: scan.scan_number,
            intensity: extractIntensity(ion.heavy_mz, scan.mz_array, scan.intensity_array)
        }))
    );

    // 5. Update Plotly charts
    updatePlots(ms1HeavyTrace, ms2HeavyTraces, heavyIons);
}
```

### 5.5 SILAC Presets

| Preset | K Delta (Da) | R Delta (Da) | Description |
|--------|-------------|-------------|-------------|
| None | — | — | Light traces only, hide heavy |
| Standard SILAC (K8R10) | 8.014199 | 10.008269 | ¹³C₆¹⁵N₂-Lys + ¹³C₆¹⁵N₄-Arg |
| Medium SILAC (K4R6) | 4.025107 | 6.020129 | ²H₄-Lys + ¹³C₆-Arg |
| Custom | user input | user input | Free-form K/R delta entry |

### 5.6 Display Modes

| Mode | Description |
|------|-------------|
| Overlay | Light (solid lines) + Heavy (dashed lines, same color) on one plot |
| Light Only | Only light traces visible |
| Heavy Only | Only heavy traces visible |
| Split | Two vertically stacked subplots: light on top, heavy below |

## 6. MCP Tool Changes

### 6.1 Enhanced `annotate_spectrum`

New optional parameters added to `AnnotateSpectrumInput`:

```rust
/// Whether to include XIC chromatogram below annotation. Default: true for mzML, false for MGF.
include_xic: Option<bool>,
/// Number of DIA/DDA cycles before/after target scan. Default: 5.
n_cycles: Option<u32>,
/// Number of top fragment ions to display in XIC. Default: 6.
top_n_ions: Option<usize>,
/// SILAC label type for pre-computed heavy traces. Default: standard SILAC.
label_type: Option<LabelType>,
/// Embed raw scan data for interactive SILAC in browser. Default: true.
embed_raw_scans: Option<bool>,
```

**Default behavior change**: When `include_xic` is not specified and the input file is mzML, XIC is included automatically with standard SILAC pre-computation. MGF files skip XIC (no multi-scan data).

### 6.2 `extract_xic` — Unchanged

The standalone `extract_xic` tool continues to work as before for XIC-only use cases.

### 6.3 Output

When `include_xic=true`:
- Output file uses the unified template
- Default name: `output/annotation_scan{N}.html` (same as before)

When `include_xic=false`:
- Output file uses the existing annotation-only template
- Identical to current behavior

## 7. Rust Implementation Changes

### 7.1 New Files

| File | Purpose |
|------|---------|
| `crates/report/src/unified_visualize.rs` | `render_unified_html()` function |
| `crates/report/templates/unified.html` | Combined HTML template |
| `crates/report/src/types.rs` | `UnifiedViewData`, `RawScanData`, `IonMetadataEntry`, `PeptideInfo` structs |

### 7.2 Modified Files

| File | Change |
|------|--------|
| `crates/report/src/lib.rs` | Export new modules |
| `crates/xic/src/extract.rs` | New `extract_xic_with_raw()` returning `(XicData, RawScanData)` |
| `crates/xic/src/lib.rs` | Export `RawScanData`, `RawScan` types (or keep in report crate) |
| `crates/mcp-server/src/tools.rs` | `annotate_spectrum` handler gains XIC integration |

### 7.3 XIC Extract Changes

New function in `extract.rs`:

```rust
pub fn extract_xic_with_raw(
    file_path: &Path,
    target_scan: u32,
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    params: &ExtractionParams,
    ms1_mz_window: f64,  // ±Da for MS1 raw peak trimming (default: 20.0)
) -> Result<(XicData, RawScanData, Vec<IonMetadataEntry>), XicError>
```

During the Pass 1 scan loop, raw peak arrays are captured alongside XIC extraction. MS1 peaks are trimmed to `precursor_mz ± ms1_mz_window` to control data volume.

### 7.4 Ion Metadata Generation

After `build_target_ions()`, compute K/R counts for each fragment:

```rust
fn compute_ion_metadata(ions: &[TargetIon], peptide: &str) -> Vec<IonMetadataEntry> {
    let chars: Vec<char> = peptide.chars().collect();
    let n = chars.len();
    ions.iter().map(|ion| {
        let fragment_chars = match ion.ion_type {
            IonType::B => &chars[..(ion.ion_number as usize).min(n)],
            IonType::Y => &chars[n.saturating_sub(ion.ion_number as usize)..],
            IonType::Precursor => &chars[..],
        };
        IonMetadataEntry {
            label: ion.label.clone(),
            ion_type: ion.ion_type,
            ion_number: ion.ion_number,
            charge: ion.charge,
            light_mz: ion.mz,
            k_count: fragment_chars.iter().filter(|&&c| c == 'K').count() as u32,
            r_count: fragment_chars.iter().filter(|&&c| c == 'R').count() as u32,
        }
    }).collect()
}
```

## 8. Edge Cases

| Case | Handling |
|------|----------|
| Peptide has no K or R | Hide SILAC controls; show info message "No SILAC label sites (no K/R)" |
| MGF format input | `include_xic` defaults to `false`; output annotation-only HTML |
| DDA data (no MS1 scans) | MS1 XIC panel shows "No MS1 data available"; MS2 XIC works normally |
| XIC extraction fails | Log warning; output annotation-only HTML with note |
| `include_xic=false` | Output existing annotation.html template (full backward compat) |
| Raw scan data > 5 MB | Truncate MS2 scans to 1000 peaks each (top by intensity) |
| No matching MS2 scans in window | XIC panels show "No matching scans found" |
| Zero-intensity heavy traces | Show flat line (indicates heavy signal absent — valid information) |

## 9. Testing Strategy

- **Unit tests**: `IonMetadataEntry` K/R counting, `RawScan` serialization, JS extraction algorithm parity
- **Integration tests**: `render_unified_html()` produces valid HTML with all three data sections
- **Parity test**: JS `extractIntensity()` results match Rust `extract_intensity()` for same inputs
- **Template tests**: HTML with mock data renders without JS errors
- **Edge case tests**: No K/R, no MS1, MGF fallback, empty scans
