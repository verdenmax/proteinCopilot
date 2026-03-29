//! Spectrum reader trait — unified interface for all spectrum file formats.

use std::path::Path;

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;

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
}
