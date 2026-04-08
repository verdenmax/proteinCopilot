//! XIC extraction core logic.
//!
//! Implements the 1.5-pass extraction algorithm:
//! - Pass 0: read target scan to get RT and isolation window
//! - Pass 1: stream all spectra, extracting intensities for target ions

use crate::IonType;

/// A target ion for XIC extraction.
#[derive(Debug, Clone)]
pub struct TargetIon {
    pub label: String,
    pub ion_type: IonType,
    pub ion_number: u32,
    pub charge: u32,
    pub mz: f64,
}

// Full extraction logic added in Task 3.
