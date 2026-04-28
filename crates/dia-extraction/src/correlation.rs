//! MS1 ↔ MS2 spectrum correlation logic.
//!
//! Correlates each MS2 spectrum to the most appropriate MS1 spectrum using a
//! three-level fallback strategy:
//! 1. **source_scan** — explicit reference from precursor metadata (most reliable)
//! 2. **scan order** — previous MS1 in acquisition order
//! 3. **retention time** — closest MS1 by absolute RT difference (last resort)

use protein_copilot_core::spectrum::Spectrum;

/// Correlates each MS2 spectrum to an MS1 spectrum using three-level fallback.
///
/// Returns a `Vec` with the same length as `ms2_spectra`. Each element is the
/// index into `ms1_spectra` of the correlated MS1, or `None` if no correlation
/// could be established (e.g., `ms1_spectra` is empty).
pub fn correlate_ms1_ms2(
    ms1_spectra: &[&Spectrum],
    ms2_spectra: &[&Spectrum],
) -> Vec<Option<usize>> {
    if ms1_spectra.is_empty() {
        return vec![None; ms2_spectra.len()];
    }

    let result: Vec<Option<usize>> = ms2_spectra
        .iter()
        .map(|ms2| correlate_single(ms1_spectra, ms2))
        .collect();

    tracing::info!(ms1 = ms1_spectra.len(), ms2 = ms2_spectra.len(), "MS1-MS2 correlation complete");

    result
}

/// Correlates a single MS2 spectrum to an MS1 using the three-level fallback.
///
/// Returns `(Option<usize>, &'static str)` — the MS1 index and the method name
/// ("source_scan", "scan_order", "rt_nearest", or "none").
pub fn correlate_single_with_method(
    ms1_spectra: &[&Spectrum],
    ms2: &Spectrum,
) -> (Option<usize>, &'static str) {
    // Level 1: source_scan from precursor metadata
    if let Some(source_scan) = ms2.precursors.first().and_then(|p| p.source_scan) {
        if let Some(idx) = ms1_spectra
            .iter()
            .position(|ms1| ms1.scan_number == source_scan)
        {
            return (Some(idx), "source_scan");
        }
    }

    // Level 2: largest MS1 scan_number < MS2 scan_number
    let by_scan_order = ms1_spectra
        .iter()
        .enumerate()
        .filter(|(_, ms1)| ms1.scan_number < ms2.scan_number)
        .max_by_key(|(_, ms1)| ms1.scan_number);

    if let Some((idx, _)) = by_scan_order {
        return (Some(idx), "scan_order");
    }

    // Level 3: closest MS1 by retention time (skip NaN RT values)
    let by_rt = ms1_spectra
        .iter()
        .enumerate()
        .filter(|(_, s)| s.retention_time_min.is_finite() && ms2.retention_time_min.is_finite())
        .min_by(|(_, a), (_, b)| {
            let da = (a.retention_time_min - ms2.retention_time_min).abs();
            let db = (b.retention_time_min - ms2.retention_time_min).abs();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(idx, _)| idx);

    match by_rt {
        Some(idx) => (Some(idx), "rt_nearest"),
        None => (None, "none"),
    }
}

/// Correlates a single MS2 spectrum to an MS1 using the three-level fallback.
fn correlate_single(ms1_spectra: &[&Spectrum], ms2: &Spectrum) -> Option<usize> {
    correlate_single_with_method(ms1_spectra, ms2).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo};

    fn make_ms1(scan_number: u32, rt: f64) -> Spectrum {
        Spectrum::new(scan_number, MsLevel::MS1, rt, vec![], vec![], vec![]).unwrap()
    }

    fn make_ms2(scan_number: u32, rt: f64, source_scan: Option<u32>) -> Spectrum {
        let precursor = PrecursorInfo {
            mz: 500.0,
            charge: Some(2),
            intensity: None,
            isolation_window: Some(IsolationWindow {
                target_mz: 500.0,
                lower_offset: 12.5,
                upper_offset: 12.5,
            }),
            source_scan,
        };
        Spectrum::new(
            scan_number,
            MsLevel::MS2,
            rt,
            vec![precursor],
            vec![],
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn test_correlate_by_source_scan() {
        let ms1 = make_ms1(5, 10.0);
        let ms2 = make_ms2(6, 10.5, Some(5));

        let result = correlate_ms1_ms2(&[&ms1], &[&ms2]);
        assert_eq!(result, vec![Some(0)]);
    }

    #[test]
    fn test_correlate_by_scan_order() {
        let ms1_a = make_ms1(1, 1.0);
        let ms1_b = make_ms1(5, 5.0);
        let ms1_c = make_ms1(15, 15.0);
        let ms2 = make_ms2(10, 10.0, None);

        let result = correlate_ms1_ms2(&[&ms1_a, &ms1_b, &ms1_c], &[&ms2]);
        // Largest MS1 scan < 10 is scan 5 at index 1
        assert_eq!(result, vec![Some(1)]);
    }

    #[test]
    fn test_correlate_by_rt() {
        // All MS1 scans have scan_number > MS2 scan, so scan-order fallback
        // won't match, triggering the RT fallback.
        let ms1_a = make_ms1(200, 10.0);
        let ms1_b = make_ms1(201, 20.0);
        let ms1_c = make_ms1(202, 30.0);
        let ms2 = make_ms2(100, 21.0, None);

        let result = correlate_ms1_ms2(&[&ms1_a, &ms1_b, &ms1_c], &[&ms2]);
        // Closest RT to 21.0 is 20.0 at index 1
        assert_eq!(result, vec![Some(1)]);
    }

    #[test]
    fn test_correlate_empty_ms1() {
        let ms2_a = make_ms2(10, 10.0, Some(5));
        let ms2_b = make_ms2(20, 20.0, None);

        let result = correlate_ms1_ms2(&[], &[&ms2_a, &ms2_b]);
        assert_eq!(result, vec![None, None]);
    }

    #[test]
    fn test_correlate_multiple_ms2() {
        let ms1_a = make_ms1(1, 1.0);
        let ms1_b = make_ms1(5, 5.0);
        let ms1_c = make_ms1(200, 30.0);

        // MS2 #1: matched by source_scan → ms1_c (scan 200, index 2)
        let ms2_source = make_ms2(201, 30.5, Some(200));
        // MS2 #2: matched by scan order → ms1_b (scan 5, largest < 10, index 1)
        let ms2_scan = make_ms2(10, 10.0, None);
        // MS2 #3: no MS1 with scan < 3, no source_scan → RT fallback
        //         closest RT to 2.0 is ms1_a at RT 1.0 (index 0)
        let ms2_rt = make_ms2(3, 2.0, None);

        let result = correlate_ms1_ms2(
            &[&ms1_a, &ms1_b, &ms1_c],
            &[&ms2_source, &ms2_scan, &ms2_rt],
        );
        assert_eq!(result, vec![Some(2), Some(1), Some(0)]);
    }
}
