//! Spectrum reader trait — unified interface for all spectrum file formats.

use std::path::Path;

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;

/// MS2 scan metadata returned by [`SpectrumReader::list_ms2_meta`].
#[derive(Debug, Clone)]
pub struct Ms2ScanMeta {
    pub scan_number: u32,
    pub rt_min: f64,
    pub isolation_window: Option<(f64, f64, f64)>,
}

/// Unified interface for reading spectrum files.
///
/// Each supported format (mgf, mzML) implements this trait.
/// Use [`crate::create_reader`] to obtain the appropriate reader
/// for a given file.
pub trait SpectrumReader: Send + Sync {
    /// Reads all spectra from the file.
    ///
    /// For large files, consider using [`Self::read_summary`] first to check
    /// data characteristics without loading all peak data into memory.
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError>;

    /// Computes a statistical summary of the spectrum file.
    ///
    /// Uses streaming parsing to avoid loading all spectra into memory
    /// simultaneously. This is the primary input for AI-driven parameter
    /// recommendation.
    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError>;

    /// Reads a single spectrum by scan number.
    ///
    /// Returns [`SpectrumIoError::ScanNotFound`] if the scan number
    /// does not exist in the file.
    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError>;

    /// Streams spectra one at a time, calling `handler` for each.
    ///
    /// The handler returns `Ok(true)` to continue or `Ok(false)` to stop early.
    /// Returns the number of spectra processed (including the one that stopped).
    ///
    /// This avoids loading all spectra into memory at once, which is important
    /// for large DIA files when only extracting specific ion chromatograms.
    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError>;

    /// Returns metadata for all MS2 scans: `(scan_number, rt_min, isolation_window)`.
    ///
    /// Default implementation reads all spectra (slow for large files).
    /// [`crate::IndexedMzMLReader`] overrides to read directly from the
    /// PCIX v2 disk cache index — zero I/O, sub-millisecond.
    fn list_ms2_meta(
        &self,
        path: &Path,
    ) -> Result<Vec<Ms2ScanMeta>, SpectrumIoError> {
        use protein_copilot_core::spectrum::MsLevel;

        let spectra = self.read_all(path)?;
        Ok(spectra
            .iter()
            .filter(|s| s.ms_level == MsLevel::MS2)
            .map(|s| {
                let iso = s.precursors.first().and_then(|p| {
                    p.isolation_window
                        .as_ref()
                        .map(|w| (w.target_mz, w.lower_offset, w.upper_offset))
                });
                Ms2ScanMeta {
                    scan_number: s.scan_number,
                    rt_min: s.retention_time_min,
                    isolation_window: iso,
                }
            })
            .collect())
    }

    /// Find the best MS2 scan matching a given RT and precursor m/z.
    ///
    /// Default implementation reads all spectra (slow for large files).
    /// [`crate::IndexedMzMLReader`] overrides with O(log N) binary search.
    ///
    /// Returns `(scan_number, rt_delta_min)` or `None`.
    fn find_by_rt(
        &self,
        path: &Path,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Result<Option<(u32, f64)>, SpectrumIoError> {
        use protein_copilot_core::spectrum::MsLevel;

        let spectra = self.read_all(path)?;

        let mut best: Option<(u32, f64)> = None;
        for spec in &spectra {
            if spec.ms_level != MsLevel::MS2 {
                continue;
            }
            let delta_min = spec.retention_time_min - rt_min;
            if delta_min.abs() > rt_tolerance_min {
                continue;
            }
            if let Some(p) = spec.precursors.first() {
                if let Some(w) = &p.isolation_window {
                    let low = w.target_mz - w.lower_offset;
                    let high = w.target_mz + w.upper_offset;
                    if precursor_mz < low || precursor_mz > high {
                        continue;
                    }
                }
            }
            match &best {
                None => best = Some((spec.scan_number, delta_min)),
                Some((_, bd)) => {
                    if delta_min.abs() < bd.abs() {
                        best = Some((spec.scan_number, delta_min));
                    }
                }
            }
        }
        Ok(best)
    }
}
