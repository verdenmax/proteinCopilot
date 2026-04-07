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

        if s.retention_time_sec < self.rt_min {
            self.rt_min = s.retention_time_sec;
        }
        if s.retention_time_sec > self.rt_max {
            self.rt_max = s.retention_time_sec;
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
            Some(protein_copilot_core::util::compute_median(&self.isolation_widths))
        };

        let summary = SpectrumSummary {
            file_path: path.to_string_lossy().to_string(),
            format,
            total_spectra: self.total,
            ms1_count: self.ms1_count,
            ms2_count: self.ms2_count,
            mz_range: [self.mz_min, self.mz_max],
            rt_range_sec: [self.rt_min, self.rt_max],
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
