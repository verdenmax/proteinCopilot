//! Precursor extraction trait for DIA analysis.

use protein_copilot_core::spectrum::{IsolationWindow, PrecursorInfo, Spectrum};

/// Trait for precursor extraction algorithms.
///
/// Implementations analyze MS1 spectra within a given isolation window
/// to identify candidate precursor ions for DIA MS2 spectra.
pub trait PrecursorExtractor: Send + Sync {
    /// Extract candidate precursors from an MS1 spectrum within the given isolation window.
    ///
    /// Returns a list of identified precursor ions with m/z, charge, and intensity.
    fn extract(&self, ms1: &Spectrum, isolation_window: &IsolationWindow) -> Vec<PrecursorInfo>;
}
