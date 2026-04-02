//! DDA/DIA acquisition mode auto-detection.

use protein_copilot_core::spectrum::{AcquisitionMode, MsLevel, Spectrum};

/// Detects whether spectra were acquired in DDA or DIA mode based on
/// the median isolation window width of MS2 spectra.
///
/// Returns `AcquisitionMode::Unknown` when no MS2 spectra exist or none
/// carry isolation window information.
pub fn detect_acquisition_mode(spectra: &[Spectrum], threshold_da: f64) -> AcquisitionMode {
    let widths: Vec<f64> = spectra
        .iter()
        .filter(|s| s.ms_level == MsLevel::MS2)
        .flat_map(|s| &s.precursors)
        .filter_map(|p| {
            p.isolation_window
                .as_ref()
                .map(|w| w.lower_offset + w.upper_offset)
        })
        .collect();

    if widths.is_empty() {
        return AcquisitionMode::Unknown;
    }

    let median = median_f64(&widths);

    if median > threshold_da {
        AcquisitionMode::DIA
    } else {
        AcquisitionMode::DDA
    }
}

/// Separates spectra into MS1 and MS2 groups by reference.
pub fn separate_by_ms_level(spectra: &[Spectrum]) -> (Vec<&Spectrum>, Vec<&Spectrum>) {
    let mut ms1 = Vec::new();
    let mut ms2 = Vec::new();
    for s in spectra {
        match s.ms_level {
            MsLevel::MS1 => ms1.push(s),
            MsLevel::MS2 => ms2.push(s),
            MsLevel::Other(_) => {}
        }
    }
    (ms1, ms2)
}

/// Computes the median of a non-empty slice of `f64` values.
///
/// # Panics
///
/// Panics if `values` is empty.
fn median_f64(values: &[f64]) -> f64 {
    assert!(!values.is_empty(), "median_f64 requires non-empty input");
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo, Spectrum};

    fn make_ms2(scan: u32, lower: f64, upper: f64) -> Spectrum {
        Spectrum::new(
            scan,
            MsLevel::MS2,
            1.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: None,
                isolation_window: Some(IsolationWindow {
                    target_mz: 500.0,
                    lower_offset: lower,
                    upper_offset: upper,
                }),
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        )
        .unwrap()
    }

    fn make_ms2_no_window(scan: u32) -> Spectrum {
        Spectrum::new(
            scan,
            MsLevel::MS2,
            1.0,
            vec![PrecursorInfo {
                mz: 500.0,
                charge: Some(2),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            vec![100.0],
            vec![1000.0],
        )
        .unwrap()
    }

    fn make_ms1(scan: u32) -> Spectrum {
        Spectrum::new(scan, MsLevel::MS1, 1.0, vec![], vec![100.0], vec![1000.0]).unwrap()
    }

    #[test]
    fn test_detect_dda() {
        let spectra: Vec<Spectrum> = (1..=5).map(|i| make_ms2(i, 1.0, 1.0)).collect();
        assert_eq!(detect_acquisition_mode(&spectra, 5.0), AcquisitionMode::DDA);
    }

    #[test]
    fn test_detect_dia() {
        let spectra: Vec<Spectrum> = (1..=5).map(|i| make_ms2(i, 12.5, 12.5)).collect();
        assert_eq!(detect_acquisition_mode(&spectra, 5.0), AcquisitionMode::DIA);
    }

    #[test]
    fn test_detect_no_window() {
        let spectra: Vec<Spectrum> = (1..=3).map(|i| make_ms2_no_window(i)).collect();
        assert_eq!(
            detect_acquisition_mode(&spectra, 5.0),
            AcquisitionMode::Unknown
        );
    }

    #[test]
    fn test_detect_empty() {
        assert_eq!(detect_acquisition_mode(&[], 5.0), AcquisitionMode::Unknown);
    }

    #[test]
    fn test_detect_mixed() {
        // 2 narrow (2.0 Da) + 3 wide (25.0 Da) → median is 25.0 → DIA
        let mut spectra: Vec<Spectrum> = (1..=2).map(|i| make_ms2(i, 1.0, 1.0)).collect();
        spectra.extend((3..=5).map(|i| make_ms2(i, 12.5, 12.5)));
        assert_eq!(detect_acquisition_mode(&spectra, 5.0), AcquisitionMode::DIA);
    }

    #[test]
    fn test_separate_by_ms_level() {
        let spectra = vec![
            make_ms1(1),
            make_ms2(2, 1.0, 1.0),
            make_ms1(3),
            make_ms2(4, 1.0, 1.0),
        ];
        let (ms1, ms2) = separate_by_ms_level(&spectra);
        assert_eq!(ms1.len(), 2);
        assert_eq!(ms2.len(), 2);
        assert_eq!(ms1[0].scan_number, 1);
        assert_eq!(ms1[1].scan_number, 3);
        assert_eq!(ms2[0].scan_number, 2);
        assert_eq!(ms2[1].scan_number, 4);
    }
}
