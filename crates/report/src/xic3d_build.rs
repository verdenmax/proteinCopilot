//! Builder for [`Xic3dData`]: per-scan b/y annotation of every MS2 spectrum
//! in the target isolation window's ±n_cycles RT range.

use protein_copilot_core::search_params::{MassTolerance, Modification};
use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum};
use protein_copilot_search_engine::annotate::annotate_spectrum;
use protein_copilot_xic::RawScan;

use crate::error::ReportError;
use crate::xic3d_types::{Ms2ScanAnnotation, Xic3dData};

/// Builds [`Xic3dData`] from raw MS2 scans by annotating each spectrum against
/// the identified peptide's theoretical b/y ions.
///
/// `ms2_scans` are the full peak lists captured by `extract_xic_unified`
/// (target isolation window, ±n_cycles), already sorted by scan number.
/// Returns [`ReportError::EmptyMs2Window`] if `ms2_scans` is empty, or
/// [`ReportError::AnnotationError`] if a spectrum fails to annotate.
///
/// Fail-fast: any single scan that fails to annotate (e.g. an empty MS2
/// scan) aborts the whole build rather than being skipped. Numeric input
/// validation is delegated to `Spectrum::new` and `annotate_spectrum`.
#[allow(clippy::too_many_arguments)]
pub fn build_xic3d_data(
    ms2_scans: &[RawScan],
    peptide_sequence: &str,
    charge: i32,
    precursor_mz: f64,
    modifications: &[Modification],
    fragment_tolerance: &MassTolerance,
    target_scan: u32,
    source_file: &str,
) -> Result<Xic3dData, ReportError> {
    if ms2_scans.is_empty() {
        return Err(ReportError::EmptyMs2Window { scan: target_scan });
    }

    let mut scans = Vec::with_capacity(ms2_scans.len());
    for raw in ms2_scans {
        let spectrum = Spectrum::new(
            raw.scan_number,
            MsLevel::MS2,
            raw.retention_time_min,
            vec![PrecursorInfo {
                mz: precursor_mz,
                charge: Some(charge),
                intensity: None,
                isolation_window: None,
                source_scan: None,
            }],
            raw.mz_array.clone(),
            raw.intensity_array.clone(),
        )
        .map_err(|e| ReportError::AnnotationError {
            scan: raw.scan_number,
            detail: e.to_string(),
        })?;

        let annotation = annotate_spectrum(
            &spectrum,
            peptide_sequence,
            charge,
            fragment_tolerance,
            modifications,
            Vec::new(),
            false,
            false,
        )
        .map_err(|e| ReportError::AnnotationError {
            scan: raw.scan_number,
            detail: e.to_string(),
        })?;

        scans.push(Ms2ScanAnnotation {
            scan_number: raw.scan_number,
            retention_time_min: raw.retention_time_min,
            total_peaks: raw.mz_array.len(),
            is_target: raw.scan_number == target_scan,
            annotation,
        });
    }

    Ok(Xic3dData {
        peptide_sequence: peptide_sequence.to_string(),
        charge,
        precursor_mz,
        target_scan,
        source_file: source_file.to_string(),
        mz_tolerance: fragment_tolerance.clone(),
        scans,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use protein_copilot_core::search_params::ToleranceUnit;
    use protein_copilot_search_engine::matching::{
        generate_b_ions_with_charge, generate_y_ions_with_charge,
    };
    use protein_copilot_xic::RawScan;

    fn tol() -> MassTolerance {
        MassTolerance {
            value: 20.0,
            unit: ToleranceUnit::Ppm,
        }
    }

    /// A RawScan whose peaks include the first 3 b- and 3 y-ions (1+) of `pep`,
    /// plus two non-matching peaks. Peaks sorted ascending by m/z.
    fn raw_scan_with_matches(scan: u32, rt: f64, pep: &str) -> RawScan {
        let b = generate_b_ions_with_charge(pep, &[], 1);
        let y = generate_y_ions_with_charge(pep, &[], 1);
        let mut mz: Vec<f64> = Vec::new();
        mz.extend(b.iter().take(3).copied());
        mz.extend(y.iter().take(3).copied());
        mz.push(150.0);
        mz.push(1234.5);
        mz.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let intensity = vec![1000.0; mz.len()];
        RawScan {
            scan_number: scan,
            retention_time_min: rt,
            mz_array: mz,
            intensity_array: intensity,
        }
    }

    #[test]
    fn build_annotates_each_scan_and_marks_target() {
        let pep = "PEPTIDEK";
        let scans = vec![
            raw_scan_with_matches(100, 19.9, pep),
            raw_scan_with_matches(101, 20.0, pep),
            raw_scan_with_matches(102, 20.1, pep),
        ];
        let data = build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol(), 101, "s.mzML").unwrap();
        assert_eq!(data.scans.len(), 3);
        assert_eq!(data.scans[0].total_peaks, scans[0].mz_array.len());
        assert_eq!(data.scans.iter().filter(|s| s.is_target).count(), 1);
        assert!(data.scans[1].is_target);
        assert!(
            data.scans[1].annotation.matched_ions >= 6,
            "expected >=6 matched, got {}",
            data.scans[1].annotation.matched_ions
        );
    }

    #[test]
    fn build_marks_matched_peaks() {
        let pep = "PEPTIDEK";
        let scans = vec![raw_scan_with_matches(100, 20.0, pep)];
        let data = build_xic3d_data(&scans, pep, 2, 460.0, &[], &tol(), 100, "s.mzML").unwrap();
        let annotated = data.scans[0]
            .annotation
            .peaks
            .iter()
            .filter(|p| p.annotation.is_some())
            .count();
        assert!(
            annotated >= 6,
            "expected >=6 annotated peaks, got {annotated}"
        );
    }

    #[test]
    fn build_empty_window_errors() {
        let err =
            build_xic3d_data(&[], "PEPTIDEK", 2, 460.0, &[], &tol(), 100, "s.mzML").unwrap_err();
        assert!(matches!(err, ReportError::EmptyMs2Window { scan: 100 }));
    }

    #[test]
    fn build_scan_with_no_peaks_errors() {
        let scans = vec![RawScan {
            scan_number: 100,
            retention_time_min: 20.0,
            mz_array: vec![],
            intensity_array: vec![],
        }];
        let err =
            build_xic3d_data(&scans, "PEPTIDEK", 2, 460.0, &[], &tol(), 100, "s.mzML").unwrap_err();
        assert!(matches!(
            err,
            ReportError::AnnotationError { scan: 100, .. }
        ));
    }
}
