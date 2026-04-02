//! # ProteinCopilot DIA Extraction
//!
//! Library crate for extracting candidate precursor ions from DIA
//! (Data-Independent Acquisition) mass spectrometry data.
//!
//! DIA spectra have wide isolation windows containing multiple co-fragmented
//! precursors. This module detects candidate precursors from MS1 spectra
//! using isotope pattern analysis and associates them with DIA MS2 spectra.

pub mod config;
pub mod correlation;
pub mod detection;
pub mod error;
pub mod extractor;
pub mod isotope;

pub use config::{DiaExtractionConfig, DiaExtractionResult, ExtractionStats};
pub use error::DiaExtractionError;
pub use extractor::PrecursorExtractor;

use protein_copilot_core::spectrum::{AcquisitionMode, Spectrum};
use std::collections::HashMap;

/// Extract candidate precursor ions from DIA mass spectrometry data.
///
/// Analyzes MS1 spectra to identify isotope patterns within DIA isolation windows,
/// then populates MS2 spectra with the extracted precursor candidates.
///
/// For DDA data, returns the original MS2 spectra unchanged.
pub fn extract_dia_precursors(
    spectra: &[Spectrum],
    extractor: &dyn PrecursorExtractor,
    config: &DiaExtractionConfig,
) -> Result<DiaExtractionResult, DiaExtractionError> {
    let (ms1_refs, ms2_refs) = detection::separate_by_ms_level(spectra);

    if ms2_refs.is_empty() {
        return Err(DiaExtractionError::NoMs2Spectra);
    }

    let detected_mode = config
        .acquisition_mode
        .unwrap_or_else(|| detection::detect_acquisition_mode(spectra, config.dia_threshold_da));

    if detected_mode == AcquisitionMode::DDA || detected_mode == AcquisitionMode::Unknown {
        let enhanced_spectra: Vec<Spectrum> = ms2_refs.iter().map(|s| (*s).clone()).collect();
        let stats = ExtractionStats {
            ms1_count: ms1_refs.len() as u32,
            ms2_count: ms2_refs.len() as u32,
            total_precursors_extracted: 0,
            avg_precursors_per_ms2: 0.0,
            charge_distribution: HashMap::new(),
        };
        return Ok(DiaExtractionResult {
            detected_mode,
            enhanced_spectra,
            stats,
        });
    }

    // DIA mode
    if ms1_refs.is_empty() {
        return Err(DiaExtractionError::NoMs1Spectra);
    }

    let ms1_indices = correlation::correlate_ms1_ms2(&ms1_refs, &ms2_refs);

    let mut enhanced_spectra = Vec::with_capacity(ms2_refs.len());
    for (ms2, ms1_idx) in ms2_refs.iter().zip(ms1_indices.iter()) {
        let mut cloned = (*ms2).clone();
        if let Some(idx) = ms1_idx {
            if let Some(iw) = ms2
                .precursors
                .first()
                .and_then(|p| p.isolation_window.as_ref())
            {
                let candidates = extractor.extract(ms1_refs[*idx], iw);
                cloned.precursors = candidates;
            }
        }
        enhanced_spectra.push(cloned);
    }

    let mut total_precursors_extracted: u32 = 0;
    let mut charge_distribution: HashMap<i32, u32> = HashMap::new();
    for s in &enhanced_spectra {
        total_precursors_extracted += s.precursors.len() as u32;
        for p in &s.precursors {
            if let Some(charge) = p.charge {
                *charge_distribution.entry(charge).or_insert(0) += 1;
            }
        }
    }

    let ms2_count = ms2_refs.len() as u32;
    let avg_precursors_per_ms2 = if ms2_count > 0 {
        total_precursors_extracted as f64 / ms2_count as f64
    } else {
        0.0
    };

    let stats = ExtractionStats {
        ms1_count: ms1_refs.len() as u32,
        ms2_count,
        total_precursors_extracted,
        avg_precursors_per_ms2,
        charge_distribution,
    };

    Ok(DiaExtractionResult {
        detected_mode,
        enhanced_spectra,
        stats,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotope::IsotopePatternExtractor;
    use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo};

    fn make_ms1(scan: u32, rt: f64, mz_array: Vec<f64>, intensity_array: Vec<f64>) -> Spectrum {
        Spectrum::new(scan, MsLevel::MS1, rt, vec![], mz_array, intensity_array).unwrap()
    }

    fn make_ms2_dia(
        scan: u32,
        rt: f64,
        target_mz: f64,
        lower: f64,
        upper: f64,
        source_scan: Option<u32>,
    ) -> Spectrum {
        Spectrum::new(
            scan,
            MsLevel::MS2,
            rt,
            vec![PrecursorInfo {
                mz: target_mz,
                charge: None,
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz,
                    lower_offset: lower,
                    upper_offset: upper,
                }),
                source_scan,
            }],
            vec![200.0, 300.0, 400.0],
            vec![500.0, 800.0, 300.0],
        )
        .unwrap()
    }

    fn make_ms2_dda(scan: u32, rt: f64, mz: f64) -> Spectrum {
        Spectrum::new(
            scan,
            MsLevel::MS2,
            rt,
            vec![PrecursorInfo {
                mz,
                charge: Some(2),
                intensity: Some(1000.0),
                isolation_window: Some(IsolationWindow {
                    target_mz: mz,
                    lower_offset: 1.0,
                    upper_offset: 1.0,
                }),
                source_scan: None,
            }],
            vec![200.0, 300.0],
            vec![500.0, 300.0],
        )
        .unwrap()
    }

    #[test]
    fn test_extract_dia_full_pipeline() {
        // MS1 with two isotope clusters
        let ms1 = make_ms1(
            1,
            10.0,
            vec![500.0, 500.502, 501.003, 600.0, 600.335, 600.669],
            vec![1000.0, 800.0, 400.0, 2000.0, 1500.0, 800.0],
        );

        // Two DIA MS2 spectra with wide isolation windows
        let ms2_a = make_ms2_dia(2, 10.1, 550.0, 100.0, 100.0, Some(1));
        let ms2_b = make_ms2_dia(3, 10.2, 600.0, 50.0, 50.0, Some(1));

        let spectra = vec![ms1, ms2_a, ms2_b];
        let extractor = IsotopePatternExtractor::default();
        let config = DiaExtractionConfig::default();

        let result = extract_dia_precursors(&spectra, &extractor, &config).unwrap();

        assert_eq!(result.detected_mode, AcquisitionMode::DIA);
        assert_eq!(result.enhanced_spectra.len(), 2);

        // Each spectrum should have at least 1 extracted precursor
        for s in &result.enhanced_spectra {
            assert!(
                !s.precursors.is_empty(),
                "Scan {} should have extracted precursors",
                s.scan_number
            );
        }
        assert!(result.stats.total_precursors_extracted > 0);
    }

    #[test]
    fn test_dda_passthrough() {
        let ms1 = make_ms1(1, 10.0, vec![500.0], vec![1000.0]);
        let ms2_a = make_ms2_dda(2, 10.1, 500.0);
        let ms2_b = make_ms2_dda(3, 10.2, 600.0);

        let spectra = vec![ms1, ms2_a.clone(), ms2_b.clone()];
        let extractor = IsotopePatternExtractor::default();
        let config = DiaExtractionConfig::default();

        let result = extract_dia_precursors(&spectra, &extractor, &config).unwrap();

        assert_eq!(result.detected_mode, AcquisitionMode::DDA);
        assert_eq!(result.enhanced_spectra.len(), 2);
        assert_eq!(result.stats.total_precursors_extracted, 0);

        // Spectra should be unchanged
        assert_eq!(result.enhanced_spectra[0].scan_number, ms2_a.scan_number);
        assert_eq!(result.enhanced_spectra[0].precursors, ms2_a.precursors);
        assert_eq!(result.enhanced_spectra[1].scan_number, ms2_b.scan_number);
        assert_eq!(result.enhanced_spectra[1].precursors, ms2_b.precursors);
    }

    #[test]
    fn test_pseudo_spectra_expansion() {
        let spectrum = Spectrum::new(
            10,
            MsLevel::MS2,
            20.0,
            vec![
                PrecursorInfo {
                    mz: 500.0,
                    charge: Some(2),
                    intensity: Some(1000.0),
                    isolation_window: None,
                    source_scan: None,
                },
                PrecursorInfo {
                    mz: 600.0,
                    charge: Some(3),
                    intensity: Some(800.0),
                    isolation_window: None,
                    source_scan: None,
                },
                PrecursorInfo {
                    mz: 700.0,
                    charge: Some(2),
                    intensity: Some(500.0),
                    isolation_window: None,
                    source_scan: None,
                },
            ],
            vec![100.0, 200.0],
            vec![500.0, 300.0],
        )
        .unwrap();

        let dia_result = DiaExtractionResult {
            detected_mode: AcquisitionMode::DIA,
            enhanced_spectra: vec![spectrum],
            stats: ExtractionStats {
                ms1_count: 1,
                ms2_count: 1,
                total_precursors_extracted: 3,
                avg_precursors_per_ms2: 3.0,
                charge_distribution: HashMap::from([(2, 2), (3, 1)]),
            },
        };

        let pseudo = dia_result.expand_to_pseudo_spectra();
        assert_eq!(pseudo.len(), 3);
        for ps in &pseudo {
            assert_eq!(ps.scan_number, 10);
            assert_eq!(ps.precursors.len(), 1);
        }
        assert_eq!(pseudo[0].precursors[0].mz, 500.0);
        assert_eq!(pseudo[1].precursors[0].mz, 600.0);
        assert_eq!(pseudo[2].precursors[0].mz, 700.0);
    }

    #[test]
    fn test_no_ms2_error() {
        let ms1 = make_ms1(1, 10.0, vec![500.0], vec![1000.0]);
        let spectra = vec![ms1];
        let extractor = IsotopePatternExtractor::default();
        let config = DiaExtractionConfig::default();

        let result = extract_dia_precursors(&spectra, &extractor, &config);
        assert!(matches!(result, Err(DiaExtractionError::NoMs2Spectra)));
    }

    #[test]
    fn test_dia_no_ms1_error() {
        // MS2 spectra with wide isolation windows, force DIA mode
        let ms2 = make_ms2_dia(1, 10.0, 550.0, 100.0, 100.0, None);
        let spectra = vec![ms2];
        let extractor = IsotopePatternExtractor::default();
        let config = DiaExtractionConfig {
            acquisition_mode: Some(AcquisitionMode::DIA),
            dia_threshold_da: 5.0,
        };

        let result = extract_dia_precursors(&spectra, &extractor, &config);
        assert!(matches!(result, Err(DiaExtractionError::NoMs1Spectra)));
    }
}
