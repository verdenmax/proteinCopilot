//! mzML (PSI standard) spectrum reader.

use std::path::Path;

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Reader for mzML spectrum files.
///
/// Uses `quick-xml` for event-based streaming XML parsing. Binary data
/// arrays are decoded from base64 and optionally decompressed (zlib).
/// Supports both 32-bit and 64-bit float precision.
pub struct MzMLReader;

impl SpectrumReader for MzMLReader {
    fn read_all(&self, _path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        todo!("Task 1.2.3: implement mzML read_all")
    }

    fn read_summary(&self, _path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        todo!("Task 1.2.3: implement mzML read_summary")
    }

    fn read_spectrum(&self, _path: &Path, _scan: u32) -> Result<Spectrum, SpectrumIoError> {
        todo!("Task 1.2.3: implement mzML read_spectrum")
    }
}
