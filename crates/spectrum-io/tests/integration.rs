//! Integration tests for spectrum-io crate.
//!
//! Tests the full pipeline: detect_format → create_reader → read_all/read_summary
//! Also tests error handling for missing files, corrupt data, and format mismatches.

use std::path::{Path, PathBuf};

use protein_copilot_spectrum_io::{create_reader, detect_format, SpectrumReader};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
}

// ---------------------------------------------------------------------------
// Full pipeline: detect → create → read
// ---------------------------------------------------------------------------

#[test]
fn full_pipeline_mgf() {
    let path = fixtures_dir().join("small.mgf");
    let info = detect_format(&path).unwrap();
    assert_eq!(
        info.format,
        protein_copilot_core::spectrum::SpectrumFormat::Mgf
    );

    let reader = create_reader(&info);
    let spectra = reader.read_all(&path).unwrap();
    assert_eq!(spectra.len(), 10);

    // Validate every spectrum
    for s in &spectra {
        assert!(
            s.validate().is_ok(),
            "scan {} failed validation",
            s.scan_number
        );
    }

    let summary = reader.read_summary(&path).unwrap();
    assert_eq!(summary.total_spectra, 10);
    assert!(summary.validate().is_ok());
}

#[test]
fn full_pipeline_mzml() {
    let path = fixtures_dir().join("small.mzml");
    let info = detect_format(&path).unwrap();
    assert_eq!(
        info.format,
        protein_copilot_core::spectrum::SpectrumFormat::MzML
    );

    let reader = create_reader(&info);
    let spectra = reader.read_all(&path).unwrap();
    assert_eq!(spectra.len(), 10);

    for s in &spectra {
        assert!(
            s.validate().is_ok(),
            "scan {} failed validation",
            s.scan_number
        );
    }

    let summary = reader.read_summary(&path).unwrap();
    assert_eq!(summary.total_spectra, 10);
    assert!(summary.validate().is_ok());
}

// ---------------------------------------------------------------------------
// Cross-format consistency: mgf vs mzml should produce same data
// ---------------------------------------------------------------------------

#[test]
fn mgf_and_mzml_produce_same_spectra() {
    let mgf_path = fixtures_dir().join("small.mgf");
    let mzml_path = fixtures_dir().join("small.mzml");

    let mgf_info = detect_format(&mgf_path).unwrap();
    let mzml_info = detect_format(&mzml_path).unwrap();

    let mgf_spectra = create_reader(&mgf_info).read_all(&mgf_path).unwrap();
    let mzml_spectra = create_reader(&mzml_info).read_all(&mzml_path).unwrap();

    assert_eq!(mgf_spectra.len(), mzml_spectra.len());

    for (mgf, mzml) in mgf_spectra.iter().zip(mzml_spectra.iter()) {
        assert_eq!(mgf.scan_number, mzml.scan_number, "scan number mismatch");
        assert_eq!(
            mgf.num_peaks(),
            mzml.num_peaks(),
            "peak count mismatch for scan {}",
            mgf.scan_number
        );

        // Compare precursor m/z
        assert_eq!(mgf.precursors.len(), mzml.precursors.len());
        if let (Some(mp), Some(xp)) = (mgf.precursors.first(), mzml.precursors.first()) {
            assert!(
                (mp.mz - xp.mz).abs() < 1e-4,
                "precursor mz mismatch for scan {}: {} vs {}",
                mgf.scan_number,
                mp.mz,
                xp.mz
            );
            assert_eq!(
                mp.charge, xp.charge,
                "charge mismatch for scan {}",
                mgf.scan_number
            );
        }

        // Compare m/z arrays (should be identical since same data)
        for (i, (m, x)) in mgf.mz_array.iter().zip(mzml.mz_array.iter()).enumerate() {
            assert!(
                (m - x).abs() < 1e-3,
                "m/z mismatch at index {} for scan {}: {} vs {}",
                i,
                mgf.scan_number,
                m,
                x
            );
        }

        // Compare intensity arrays
        for (i, (m, x)) in mgf
            .intensity_array
            .iter()
            .zip(mzml.intensity_array.iter())
            .enumerate()
        {
            assert!(
                (m - x).abs() < 0.5,
                "intensity mismatch at index {} for scan {}: {} vs {}",
                i,
                mgf.scan_number,
                m,
                x
            );
        }
    }
}

#[test]
fn mgf_and_mzml_summaries_consistent() {
    let mgf_path = fixtures_dir().join("small.mgf");
    let mzml_path = fixtures_dir().join("small.mzml");

    let mgf_summary = create_reader(&detect_format(&mgf_path).unwrap())
        .read_summary(&mgf_path)
        .unwrap();
    let mzml_summary = create_reader(&detect_format(&mzml_path).unwrap())
        .read_summary(&mzml_path)
        .unwrap();

    assert_eq!(mgf_summary.total_spectra, mzml_summary.total_spectra);
    assert_eq!(mgf_summary.ms2_count, mzml_summary.ms2_count);
    assert_eq!(
        mgf_summary.median_peaks_per_spectrum,
        mzml_summary.median_peaks_per_spectrum
    );

    // Charge distributions should match
    for (charge, count) in &mgf_summary.precursor_charge_distribution {
        let mzml_count = mzml_summary
            .precursor_charge_distribution
            .get(charge)
            .unwrap_or(&0);
        assert_eq!(
            count, mzml_count,
            "charge {} distribution mismatch: mgf={} mzml={}",
            charge, count, mzml_count
        );
    }
}

// ---------------------------------------------------------------------------
// Error cases: missing files
// ---------------------------------------------------------------------------

#[test]
fn detect_format_nonexistent_file() {
    let err = detect_format(Path::new("/no/such/file.mgf")).unwrap_err();
    assert!(err.to_string().contains("not found"));
}

#[test]
fn read_all_nonexistent_mgf() {
    let reader = protein_copilot_spectrum_io::mgf::MgfReader;
    let err = reader.read_all(Path::new("/no/such/file.mgf")).unwrap_err();
    assert!(err.to_string().contains("not found"));
}

#[test]
fn read_all_nonexistent_mzml() {
    let reader = protein_copilot_spectrum_io::mzml::MzMLReader;
    let err = reader
        .read_all(Path::new("/no/such/file.mzml"))
        .unwrap_err();
    assert!(err.to_string().contains("not found"));
}

// ---------------------------------------------------------------------------
// Error cases: corrupt / unusual files
// ---------------------------------------------------------------------------

#[test]
fn mgf_corrupt_file_returns_zero_spectra() {
    // A file with no BEGIN IONS markers produces zero spectra (not an error)
    let reader = protein_copilot_spectrum_io::mgf::MgfReader;
    let spectra = reader
        .read_all(&fixtures_dir().join("corrupt.mgf"))
        .unwrap();
    assert_eq!(spectra.len(), 0);
}

#[test]
fn mgf_truncated_file_missing_end_ions() {
    // Truncated file missing END IONS — the incomplete block is dropped
    let reader = protein_copilot_spectrum_io::mgf::MgfReader;
    let spectra = reader
        .read_all(&fixtures_dir().join("truncated.mgf"))
        .unwrap();
    assert_eq!(spectra.len(), 0); // no complete spectrum
}

#[test]
fn mzml_no_binary_arrays() {
    // mzML spectrum with no binary data arrays — produces empty peak lists
    let reader = protein_copilot_spectrum_io::mzml::MzMLReader;
    let spectra = reader
        .read_all(&fixtures_dir().join("no_binary.mzml"))
        .unwrap();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].num_peaks(), 0);
}

#[test]
fn mzml_corrupt_xml() {
    let reader = protein_copilot_spectrum_io::mzml::MzMLReader;
    // Corrupt XML may parse partially or error — either is acceptable
    let result = reader.read_all(&fixtures_dir().join("corrupt.mzml"));
    // Should either return Ok with 0 spectra or Err
    match result {
        Ok(spectra) => assert_eq!(spectra.len(), 0),
        Err(e) => assert!(e.to_string().contains("XML") || e.to_string().contains("error")),
    }
}

// ---------------------------------------------------------------------------
// Error cases: scan not found
// ---------------------------------------------------------------------------

#[test]
fn read_spectrum_scan_not_found_mgf() {
    let reader = protein_copilot_spectrum_io::mgf::MgfReader;
    let err = reader
        .read_spectrum(&fixtures_dir().join("small.mgf"), 999)
        .unwrap_err();
    assert!(err.to_string().contains("999"));
}

#[test]
fn read_spectrum_scan_not_found_mzml() {
    let reader = protein_copilot_spectrum_io::mzml::MzMLReader;
    let err = reader
        .read_spectrum(&fixtures_dir().join("small.mzml"), 999)
        .unwrap_err();
    assert!(err.to_string().contains("999"));
}

// ---------------------------------------------------------------------------
// Unknown format detection
// ---------------------------------------------------------------------------

#[test]
fn detect_format_unknown_extension() {
    // Create a temp file with unsupported extension
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("data.raw");
    std::fs::write(&path, "fake content").unwrap();
    let err = detect_format(&path).unwrap_err();
    assert!(err.to_string().contains("detect format"));
}
