//! MGF (Mascot Generic Format) spectrum reader.

use std::path::Path;

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Reader for MGF (Mascot Generic Format) spectrum files.
///
/// MGF is a text-based format where each spectrum is delimited by
/// `BEGIN IONS` / `END IONS` blocks containing header key-value pairs
/// (PEPMASS, CHARGE, RTINSECONDS, SCANS) and m/z + intensity peak lines.
pub struct MgfReader;

impl SpectrumReader for MgfReader {
    fn read_all(&self, _path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        todo!("Task 1.2.2: implement MGF read_all")
    }

    fn read_summary(&self, _path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        todo!("Task 1.2.2: implement MGF read_summary")
    }

    fn read_spectrum(&self, _path: &Path, _scan: u32) -> Result<Spectrum, SpectrumIoError> {
        todo!("Task 1.2.2: implement MGF read_spectrum")
    }
}
