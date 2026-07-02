//! Shared utilities for spectrum I/O parsers.
//!
//! Contains helpers used by both MGF and mzML readers to avoid
//! code duplication.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use protein_copilot_core::spectrum::{MsLevel, Spectrum, SpectrumFormat, SpectrumSummary};

use crate::error::SpectrumIoError;

// ---------------------------------------------------------------------------
// On-load centroiding
// ---------------------------------------------------------------------------

/// Centroids a single profile-mode spectrum (local-maxima detection +
/// 3-point parabolic m/z refinement).
///
/// Algorithm adapted from ms2-met's `centroid_spectrum()`:
/// 1. Detect local maxima: `intensity[i] > intensity[i-1] && intensity[i] >= intensity[i+1]`
///    (asymmetric rule — plateaus resolve to the left-most interior index).
/// 2. Discard peaks below `max_intensity × rel_threshold`.
/// 3. Refine each surviving peak's m/z by 3-point parabolic interpolation:
///    `dx = 0.5 × (y_left − y_right) / (y_left − 2×y_center + y_right)`
///    Clamp `|dx| > 0.5` to 0 (degenerate fit guard).
///    Intensity is reported as the apex sample height (y_center).
///
/// Returns empty arrays when `len < 3`, `max_intensity ≤ 0`, or no peak
/// survives the threshold.
pub fn centroid_spectrum(
    mz: &[f64],
    intensity: &[f64],
    rel_threshold: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = mz.len();
    if n < 3 {
        return (Vec::new(), Vec::new());
    }

    // Step 1 — local maxima detection
    let max_intensity = intensity.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max_intensity <= 0.0 {
        return (Vec::new(), Vec::new());
    }
    let cutoff = max_intensity * rel_threshold;

    let mut peak_indices: Vec<usize> = Vec::with_capacity(n / 4);
    for i in 1..(n - 1) {
        let left = intensity[i - 1];
        let center = intensity[i];
        let right = intensity[i + 1];
        if center > left && center >= right && center >= cutoff {
            peak_indices.push(i);
        }
    }

    if peak_indices.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Step 2 — 3-point parabolic m/z refinement
    let len = peak_indices.len();
    let mut out_mz: Vec<f64> = Vec::with_capacity(len);
    let mut out_intensity: Vec<f64> = Vec::with_capacity(len);

    for &idx in &peak_indices {
        let y0 = intensity[idx - 1] as f64;
        let y1 = intensity[idx] as f64;
        let y2 = intensity[idx + 1] as f64;

        let denom = y0 - 2.0 * y1 + y2;
        let dx = if denom.abs() > 1e-12 {
            let d = 0.5 * (y0 - y2) / denom;
            if d.abs() > 0.5 { 0.0 } else { d }
        } else {
            0.0
        };

        let mz_center = mz[idx];
        let half_step = (mz[idx + 1] - mz[idx - 1]) * 0.5;
        let refined_mz = mz_center + dx * half_step;

        out_mz.push(refined_mz);
        out_intensity.push(y1);
    }

    (out_mz, out_intensity)
}

// ---------------------------------------------------------------------------
// File I/O
// ---------------------------------------------------------------------------

/// Opens a file with proper error mapping (NotFound vs other I/O errors).
pub(crate) fn open_buffered(path: &Path) -> Result<BufReader<File>, SpectrumIoError> {
    let file = File::open(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            SpectrumIoError::FileNotFound {
                path: path.to_path_buf(),
            }
        } else {
            SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            }
        }
    })?;
    Ok(BufReader::new(file))
}

// ---------------------------------------------------------------------------
// Peak array sorting
// ---------------------------------------------------------------------------

/// Sorts m/z and intensity arrays together by m/z ascending order.
///
/// No-op if already sorted. Handles NaN gracefully by treating it
/// as equal (NaN values will be caught by `Spectrum::new()` validation
/// after sorting).
pub(crate) fn sort_peaks_by_mz(mz: &mut Vec<f64>, intensity: &mut Vec<f64>) {
    // Guard against length-mismatched arrays: reordering by m/z indices would
    // panic (mz longer) or silently truncate (mz shorter) the intensity array,
    // producing mis-paired peaks that still pass downstream length checks. Leave
    // the arrays untouched so the caller's `Spectrum::new()` validation surfaces
    // a clean `ArrayLengthMismatch` (length is checked before sortedness).
    if mz.len() != intensity.len() {
        return;
    }
    if mz.windows(2).all(|w| w[0] <= w[1]) {
        return; // already sorted
    }
    let mut indices: Vec<usize> = (0..mz.len()).collect();
    // total_cmp handles NaN deterministically (sorts after all finite values).
    // NaN values will be caught by Spectrum::new() validation after sorting.
    indices.sort_by(|&a, &b| mz[a].total_cmp(&mz[b]));
    *mz = indices.iter().map(|&i| mz[i]).collect();
    *intensity = indices.iter().map(|&i| intensity[i]).collect();
}

// ---------------------------------------------------------------------------
// Summary accumulator
// ---------------------------------------------------------------------------

/// Accumulates statistics across spectra for building a [`SpectrumSummary`].
///
/// Used by both MGF and mzML readers to avoid duplicating the summary
/// aggregation logic.
pub(crate) struct SummaryAccumulator {
    total: u64,
    ms1_count: u64,
    ms2_count: u64,
    mz_min: f64,
    mz_max: f64,
    rt_min: f64,
    rt_max: f64,
    charge_dist: HashMap<i32, u64>,
    peak_counts: Vec<u32>,
    isolation_widths: Vec<f64>,
}

impl SummaryAccumulator {
    /// Creates a new accumulator with sentinel initial values.
    pub(crate) fn new() -> Self {
        Self {
            total: 0,
            ms1_count: 0,
            ms2_count: 0,
            mz_min: f64::MAX,
            mz_max: f64::MIN,
            rt_min: f64::MAX,
            rt_max: f64::MIN,
            charge_dist: HashMap::new(),
            peak_counts: Vec::new(),
            isolation_widths: Vec::new(),
        }
    }

    /// Observes a single spectrum and updates running statistics.
    pub(crate) fn observe(&mut self, s: &Spectrum) {
        self.total += 1;

        match s.ms_level {
            MsLevel::MS1 => self.ms1_count += 1,
            MsLevel::MS2 => self.ms2_count += 1,
            MsLevel::Other(_) => {}
        }

        if let Some(&first) = s.mz_array.first() {
            if first < self.mz_min {
                self.mz_min = first;
            }
        }
        if let Some(&last) = s.mz_array.last() {
            if last > self.mz_max {
                self.mz_max = last;
            }
        }

        if s.retention_time_min < self.rt_min {
            self.rt_min = s.retention_time_min;
        }
        if s.retention_time_min > self.rt_max {
            self.rt_max = s.retention_time_min;
        }

        for p in &s.precursors {
            if let Some(c) = p.charge {
                *self.charge_dist.entry(c).or_insert(0) += 1;
            }
            if let Some(iw) = &p.isolation_window {
                let width = iw.lower_offset + iw.upper_offset;
                if width.is_finite() && width > 0.0 {
                    self.isolation_widths.push(width);
                }
            }
        }

        self.peak_counts.push(s.num_peaks() as u32);
    }

    /// Finalizes the accumulated data into a validated [`SpectrumSummary`].
    pub(crate) fn into_summary(
        mut self,
        path: &Path,
        format: SpectrumFormat,
    ) -> Result<SpectrumSummary, SpectrumIoError> {
        // Handle empty file sentinels
        if self.total == 0 {
            tracing::warn!(
                path = %path.display(),
                "spectrum file contains 0 spectra — downstream analysis will produce no results"
            );
            self.mz_min = 0.0;
            self.mz_max = 0.0;
            self.rt_min = 0.0;
            self.rt_max = 0.0;
        }

        // Compute median peak count
        self.peak_counts.sort_unstable();
        let median_peaks = protein_copilot_core::util::compute_median_u32(&self.peak_counts);

        // Compute median isolation window width
        self.isolation_widths
            .sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_iw = if self.isolation_widths.is_empty() {
            None
        } else {
            Some(protein_copilot_core::util::compute_median(
                &self.isolation_widths,
            ))
        };

        let summary = SpectrumSummary {
            file_path: path.to_string_lossy().to_string(),
            format,
            total_spectra: self.total,
            ms1_count: self.ms1_count,
            ms2_count: self.ms2_count,
            mz_range: [self.mz_min, self.mz_max],
            rt_range_min: [self.rt_min, self.rt_max],
            precursor_charge_distribution: self.charge_dist,
            median_peaks_per_spectrum: median_peaks,
            median_isolation_window_da: median_iw,
        };
        summary
            .validate()
            .map_err(|e| SpectrumIoError::ValidationError {
                scan: 0,
                detail: format!("summary: {e}"),
            })?;
        Ok(summary)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sort_peaks_length_mismatch_does_not_panic_or_truncate() {
        // mz longer than intensity: reordering would panic (index OOB) without
        // the length guard. The function must leave both arrays untouched so the
        // caller's Spectrum::new() surfaces a clean ArrayLengthMismatch.
        let mut mz = vec![300.0, 100.0, 200.0];
        let mut intensity = vec![10.0, 20.0];
        sort_peaks_by_mz(&mut mz, &mut intensity);
        // Longer array must be unchanged (no truncation, no reordering).
        assert_eq!(mz, vec![300.0, 100.0, 200.0]);
        assert_eq!(intensity, vec![10.0, 20.0]);

        // mz shorter than intensity: reordering would silently truncate the
        // intensity array. Both arrays must remain untouched.
        let mut mz2 = vec![300.0, 100.0];
        let mut intensity2 = vec![10.0, 20.0, 30.0];
        sort_peaks_by_mz(&mut mz2, &mut intensity2);
        assert_eq!(mz2, vec![300.0, 100.0]);
        assert_eq!(intensity2, vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn sort_peaks_equal_length_still_sorts() {
        // Regression: equal-length arrays must still be co-sorted by m/z.
        let mut mz = vec![300.0, 100.0, 200.0];
        let mut intensity = vec![30.0, 10.0, 20.0];
        sort_peaks_by_mz(&mut mz, &mut intensity);
        assert_eq!(mz, vec![100.0, 200.0, 300.0]);
        assert_eq!(intensity, vec![10.0, 20.0, 30.0]);
    }

    // ── centroid_spectrum tests ─────────────────────────────────────────

    #[test]
    fn centroid_empty_with_short_input() {
        let mz = vec![100.0, 200.0];
        let intensity = vec![10.0, 20.0];
        let (out_mz, out_int) = centroid_spectrum(&mz, &intensity, 1e-3);
        assert!(out_mz.is_empty());
        assert!(out_int.is_empty());
    }

    #[test]
    fn centroid_empty_with_all_zero_intensity() {
        let mz = vec![100.0, 200.0, 300.0];
        let intensity = vec![0.0, 0.0, 0.0];
        let (out_mz, _out_int) = centroid_spectrum(&mz, &intensity, 1e-3);
        assert!(out_mz.is_empty());
    }

    #[test]
    fn centroid_detects_local_maximum() {
        let mz = vec![100.0, 200.0, 300.0];
        let intensity = vec![1.0, 10.0, 2.0];
        let (out_mz, out_int) = centroid_spectrum(&mz, &intensity, 1e-3);
        assert_eq!(out_mz.len(), 1);
        assert!((out_int[0] - 10.0).abs() < 1e-6);
    }

    #[test]
    fn centroid_refines_mz_by_parabolic_interpolation() {
        // Symmetric peak: m/z stays at center
        let mz = vec![100.0, 200.0, 300.0];
        let intensity = vec![1.0, 10.0, 1.0];
        let (out_mz, _out_int) = centroid_spectrum(&mz, &intensity, 0.0);
        assert!((out_mz[0] - 200.0).abs() < 1e-9);

        // Right-leaning peak
        let intensity2 = vec![1.0, 10.0, 5.0];
        let (out_mz2, _) = centroid_spectrum(&mz, &intensity2, 0.0);
        assert!(out_mz2[0] > 200.0);

        // Left-leaning peak
        let intensity3 = vec![5.0, 10.0, 1.0];
        let (out_mz3, _) = centroid_spectrum(&mz, &intensity3, 0.0);
        assert!(out_mz3[0] < 200.0);
    }

    #[test]
    fn centroid_plateau_resolves_leftmost() {
        let mz: Vec<f64> = (0..5).map(|i| (i as f64) * 100.0).collect();
        let intensity = vec![0.5, 10.0, 10.0, 10.0, 0.5];
        let (out_mz, out_int) = centroid_spectrum(&mz, &intensity, 1e-3);
        assert_eq!(out_mz.len(), 1);
        // Picks the leftmost interior index (1) of the plateau [1,2,3].
        // idx=1: y0=0.5, y1=10, y2=10 → dx=-0.5*(0.5-10)/(0.5-20+10) = -0.5*(-9.5)/(-9.5)=-0.5
        // Actually due to the asymmetric rule, only idx=1 is selected (center > left: 10>0.5; center >= right: 10>=10).
        // Parabola: y0=0.5, y1=10, y2=10 → denom = 0.5-20+10 = -9.5 → dx = 0.5*(0.5-10)/-9.5 = 0.5 → clamped ≤ 0.5, so dx=0.5
        // refined_mz = mz_center + dx * half_step = 100 + 0.5 * (200-0)/2 = 100 + 0.5*100 = 150
        assert!((out_mz[0] - 150.0).abs() < 1e-6, "got mz={}", out_mz[0]);
        assert!((out_int[0] - 10.0).abs() < 1e-6);
    }

    #[test]
    fn centroid_threshold_filters_low_intensity_peaks() {
        let mz: Vec<f64> = (0..10).map(|i| (i as f64) * 100.0).collect();
        let mut intensity = vec![1.0; 10];
        intensity[5] = 100.0;
        let (out_mz, out_int) = centroid_spectrum(&mz, &intensity, 0.01);
        assert_eq!(out_mz.len(), 1);
        assert!(out_int[0] > 50.0);
    }

    #[test]
    fn centroid_output_monotonic() {
        let mz: Vec<f64> = (0..20).map(|i| (i as f64) * 0.5).collect();
        let intensity: Vec<f64> = mz
            .iter()
            .map(|&x| {
                100.0 * (-(x - 3.0).powi(2) / 0.02).exp()
                    + 80.0 * (-(x - 7.0).powi(2) / 0.02).exp()
            })
            .collect();
        let (out_mz, _out_int) = centroid_spectrum(&mz, &intensity, 1e-3);
        assert!(
            out_mz.windows(2).all(|w| w[0] <= w[1]),
            "centroid output m/z must be monotonic"
        );
    }

    #[test]
    fn centroid_valley_not_detected() {
        // Valley (not a peak): center point is not a local maximum
        let mz = vec![100.0, 200.0, 300.0];
        let intensity = vec![100_000.0, 1.0, 100_000.0];
        let (out_mz, _out_int) = centroid_spectrum(&mz, &intensity, 1e-6);
        assert!(out_mz.is_empty());
    }

    #[test]
    fn centroid_threshold_boundary() {
        let mz = vec![100.0, 200.0, 300.0, 400.0, 500.0];
        let intensity = vec![1.0, 10.0, 5.0, 49.0, 2.0];
        let (out_mz, _) = centroid_spectrum(&mz, &intensity, 0.01);
        assert_eq!(out_mz.len(), 2);

        let (out_mz2, _) = centroid_spectrum(&mz, &intensity, 0.5);
        assert_eq!(out_mz2.len(), 1);
    }
}
