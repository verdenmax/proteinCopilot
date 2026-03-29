//! mzML (PSI standard) spectrum reader.
//!
//! Uses `quick-xml` for event-based streaming XML parsing. Binary data
//! arrays are decoded from base64 and optionally decompressed (zlib).
//! Supports both 32-bit and 64-bit float precision.
//!
//! # Relevant CV accessions
//!
//! | Accession | Meaning |
//! |-----------|---------|
//! | MS:1000511 | ms level |
//! | MS:1000016 | scan start time (seconds) |
//! | MS:1000514 | m/z array |
//! | MS:1000515 | intensity array |
//! | MS:1000521 | 32-bit float |
//! | MS:1000523 | 64-bit float |
//! | MS:1000574 | zlib compression |
//! | MS:1000576 | no compression |
//! | MS:1000744 | selected ion m/z |
//! | MS:1000041 | charge state |
//! | MS:1000042 | peak intensity |
//! | MS:1000827 | isolation window target m/z |
//! | MS:1000828 | isolation window lower offset |
//! | MS:1000829 | isolation window upper offset |

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use protein_copilot_core::spectrum::{
    IsolationWindow, MsLevel, PrecursorInfo, Spectrum, SpectrumFormat, SpectrumSummary,
};
use quick_xml::events::Event;
use quick_xml::Reader;

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Reader for mzML spectrum files.
pub struct MzMLReader;

// ---------------------------------------------------------------------------
// Internal parsing types
// ---------------------------------------------------------------------------

#[derive(Default)]
struct BinaryArrayMeta {
    is_mz: bool,
    is_intensity: bool,
    is_64bit: bool,
    is_zlib: bool,
}

#[derive(Default)]
struct SpectrumBuilder {
    scan_number: Option<u32>,
    ms_level: Option<u8>,
    rt_sec: Option<f64>,
    precursor_mz: Option<f64>,
    precursor_charge: Option<i32>,
    precursor_intensity: Option<f64>,
    isolation_target_mz: Option<f64>,
    isolation_lower: Option<f64>,
    isolation_upper: Option<f64>,
    mz_array: Vec<f64>,
    intensity_array: Vec<f64>,
}

impl SpectrumBuilder {
    fn build(self, _path: &Path) -> Result<Spectrum, SpectrumIoError> {
        let scan = self.scan_number.unwrap_or(1);
        let ms_level = match self.ms_level.unwrap_or(2) {
            1 => MsLevel::MS1,
            2 => MsLevel::MS2,
            n => MsLevel::Other(n),
        };

        let mut precursors = Vec::new();
        if let Some(mz) = self.precursor_mz {
            let isolation_window = match (
                self.isolation_target_mz,
                self.isolation_lower,
                self.isolation_upper,
            ) {
                (Some(t), Some(l), Some(u)) => Some(IsolationWindow {
                    target_mz: t,
                    lower_offset: l,
                    upper_offset: u,
                }),
                _ => None,
            };
            precursors.push(PrecursorInfo {
                mz,
                charge: self.precursor_charge,
                intensity: self.precursor_intensity,
                isolation_window,
            });
        }

        // Sort m/z + intensity arrays together by m/z ascending.
        // Some mzML files may have unsorted peaks.
        let mut mz_array = self.mz_array;
        let mut intensity_array = self.intensity_array;
        if !mz_array.windows(2).all(|w| w[0] <= w[1]) {
            let mut indices: Vec<usize> = (0..mz_array.len()).collect();
            indices.sort_by(|&a, &b| mz_array[a].partial_cmp(&mz_array[b]).unwrap());
            mz_array = indices.iter().map(|&i| mz_array[i]).collect();
            intensity_array = indices.iter().map(|&i| intensity_array[i]).collect();
        }

        Spectrum::new(
            scan,
            ms_level,
            self.rt_sec.unwrap_or(0.0),
            precursors,
            mz_array,
            intensity_array,
        )
        .map_err(|e| SpectrumIoError::ValidationError {
            scan,
            detail: e.to_string(),
        })
    }
}

// ---------------------------------------------------------------------------
// Binary data decoding
// ---------------------------------------------------------------------------

fn decode_binary_array(
    b64_text: &str,
    meta: &BinaryArrayMeta,
    path: &Path,
) -> Result<Vec<f64>, SpectrumIoError> {
    use base64::Engine;

    let raw = base64::engine::general_purpose::STANDARD
        .decode(b64_text.trim())
        .map_err(|e| SpectrumIoError::BinaryDecodeError {
            path: path.to_path_buf(),
            detail: format!("base64 decode failed: {e}"),
        })?;

    let bytes = if meta.is_zlib {
        use flate2::read::ZlibDecoder;
        use std::io::Read;
        let mut decoder = ZlibDecoder::new(&raw[..]);
        let mut decompressed = Vec::new();
        decoder
            .read_to_end(&mut decompressed)
            .map_err(|e| SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: format!("zlib decompress failed: {e}"),
            })?;
        decompressed
    } else {
        raw
    };

    if meta.is_64bit {
        if bytes.len() % 8 != 0 {
            return Err(SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: format!(
                    "64-bit array byte length {} not divisible by 8",
                    bytes.len()
                ),
            });
        }
        Ok(bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect())
    } else {
        // 32-bit float
        if bytes.len() % 4 != 0 {
            return Err(SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: format!(
                    "32-bit array byte length {} not divisible by 4",
                    bytes.len()
                ),
            });
        }
        Ok(bytes
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()) as f64)
            .collect())
    }
}

// ---------------------------------------------------------------------------
// XML attribute helpers
// ---------------------------------------------------------------------------

fn get_attr<'a>(e: &'a quick_xml::events::BytesStart<'a>, name: &[u8]) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.as_ref() == name)
        .and_then(|a| String::from_utf8(a.value.to_vec()).ok())
}

/// Extract scan number from spectrum id attribute (e.g., "scan=123").
fn parse_scan_from_id(id: &str) -> Option<u32> {
    id.split("scan=")
        .nth(1)
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse().ok())
}

// ---------------------------------------------------------------------------
// Core streaming parser
// ---------------------------------------------------------------------------

/// Streaming mzML parser. Handler returns `true` to continue, `false` to stop.
fn parse_mzml_streaming<R: BufRead, F>(
    xml_reader: &mut Reader<R>,
    path: &Path,
    mut handler: F,
) -> Result<u32, SpectrumIoError>
where
    F: FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
{
    let mut buf = Vec::new();
    let mut count: u32 = 0;

    // Parser state
    let mut in_spectrum = false;
    let mut in_precursor = false;
    let mut in_isolation_window = false;
    let mut in_selected_ion = false;
    let mut in_binary_data_array = false;
    let mut in_binary = false;
    let mut in_scan = false;

    let mut builder = SpectrumBuilder::default();
    let mut array_meta = BinaryArrayMeta::default();
    let mut binary_text = String::new();
    let mut fallback_scan: u32 = 0;

    loop {
        match xml_reader.read_event_into(&mut buf) {
            Ok(Event::Eof) => break,
            Ok(Event::Start(ref e)) => {
                let local = e.local_name();
                let tag = local.as_ref();

                match tag {
                    b"spectrum" => {
                        in_spectrum = true;
                        fallback_scan += 1;
                        builder = SpectrumBuilder::default();
                        if let Some(id) = get_attr(e, b"id") {
                            builder.scan_number = parse_scan_from_id(&id);
                        }
                        if builder.scan_number.is_none() {
                            builder.scan_number = Some(fallback_scan);
                        }
                    }
                    b"scan" if in_spectrum => {
                        in_scan = true;
                    }
                    b"precursor" if in_spectrum => {
                        in_precursor = true;
                    }
                    b"isolationWindow" if in_precursor => {
                        in_isolation_window = true;
                    }
                    b"selectedIon" if in_precursor => {
                        in_selected_ion = true;
                    }
                    b"binaryDataArray" if in_spectrum => {
                        in_binary_data_array = true;
                        array_meta = BinaryArrayMeta::default();
                        binary_text.clear();
                    }
                    b"binary" if in_binary_data_array => {
                        in_binary = true;
                        binary_text.clear();
                    }
                    _ => {}
                }
            }
            Ok(Event::Empty(ref e)) => {
                let local = e.local_name();
                let tag = local.as_ref();

                if tag == b"cvParam" {
                    let acc = get_attr(e, b"accession").unwrap_or_default();
                    let value = get_attr(e, b"value").unwrap_or_default();

                    match acc.as_str() {
                        "MS:1000511" if in_spectrum && !in_precursor => {
                            builder.ms_level = value.parse().ok();
                        }
                        "MS:1000016" if in_scan => {
                            builder.rt_sec = value.parse().ok();
                        }
                        "MS:1000827" if in_isolation_window => {
                            builder.isolation_target_mz = value.parse().ok();
                        }
                        "MS:1000828" if in_isolation_window => {
                            builder.isolation_lower = value.parse().ok();
                        }
                        "MS:1000829" if in_isolation_window => {
                            builder.isolation_upper = value.parse().ok();
                        }
                        "MS:1000744" if in_selected_ion => {
                            builder.precursor_mz = value.parse().ok();
                        }
                        "MS:1000041" if in_selected_ion => {
                            builder.precursor_charge = value.parse().ok();
                        }
                        "MS:1000042" if in_selected_ion => {
                            builder.precursor_intensity = value.parse().ok();
                        }
                        "MS:1000514" if in_binary_data_array => {
                            array_meta.is_mz = true;
                        }
                        "MS:1000515" if in_binary_data_array => {
                            array_meta.is_intensity = true;
                        }
                        "MS:1000523" if in_binary_data_array => {
                            array_meta.is_64bit = true;
                        }
                        "MS:1000521" if in_binary_data_array => {
                            array_meta.is_64bit = false;
                        }
                        "MS:1000574" if in_binary_data_array => {
                            array_meta.is_zlib = true;
                        }
                        "MS:1000576" if in_binary_data_array => {
                            array_meta.is_zlib = false;
                        }
                        _ => {}
                    }
                }
            }
            Ok(Event::End(ref e)) => {
                let local = e.local_name();
                let tag = local.as_ref();

                match tag {
                    b"spectrum" => {
                        in_spectrum = false;
                        let spectrum = builder.build(path)?;
                        let keep_going = handler(spectrum)?;
                        count += 1;
                        if !keep_going {
                            return Ok(count);
                        }
                        builder = SpectrumBuilder::default();
                    }
                    b"scan" => {
                        in_scan = false;
                    }
                    b"precursor" => {
                        in_precursor = false;
                    }
                    b"isolationWindow" => {
                        in_isolation_window = false;
                    }
                    b"selectedIon" => {
                        in_selected_ion = false;
                    }
                    b"binary" => {
                        in_binary = false;
                    }
                    b"binaryDataArray" => {
                        if !binary_text.is_empty() {
                            let decoded = decode_binary_array(&binary_text, &array_meta, path)?;
                            if array_meta.is_mz {
                                builder.mz_array = decoded;
                            } else if array_meta.is_intensity {
                                builder.intensity_array = decoded;
                            }
                        }
                        in_binary_data_array = false;
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(ref t)) => {
                if in_binary {
                    let text = t.unescape().map_err(|e| SpectrumIoError::XmlError {
                        path: path.to_path_buf(),
                        detail: format!("text unescape error: {e}"),
                    })?;
                    binary_text.push_str(&text);
                }
            }
            Err(e) => {
                return Err(SpectrumIoError::XmlError {
                    path: path.to_path_buf(),
                    detail: format!(
                        "XML parse error at position {}: {e}",
                        xml_reader.error_position()
                    ),
                });
            }
            _ => {}
        }
        buf.clear();
    }

    Ok(count)
}

// ---------------------------------------------------------------------------
// SpectrumReader implementation
// ---------------------------------------------------------------------------

impl SpectrumReader for MzMLReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        let file = std::fs::File::open(path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                SpectrumIoError::FileNotFound {
                    path: path.to_path_buf(),
                }
            } else {
                SpectrumIoError::IoError {
                    path: path.to_path_buf(),
                    source: e,
                }
            }
        })?;
        let buf_reader = std::io::BufReader::new(file);
        let mut xml_reader = Reader::from_reader(buf_reader);
        xml_reader.config_mut().trim_text(true);

        let mut spectra = Vec::new();
        parse_mzml_streaming(&mut xml_reader, path, |s| {
            spectra.push(s);
            Ok(true)
        })?;
        Ok(spectra)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        let file = std::fs::File::open(path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                SpectrumIoError::FileNotFound {
                    path: path.to_path_buf(),
                }
            } else {
                SpectrumIoError::IoError {
                    path: path.to_path_buf(),
                    source: e,
                }
            }
        })?;
        let buf_reader = std::io::BufReader::new(file);
        let mut xml_reader = Reader::from_reader(buf_reader);
        xml_reader.config_mut().trim_text(true);

        let mut total: u64 = 0;
        let mut ms1_count: u64 = 0;
        let mut ms2_count: u64 = 0;
        let mut mz_min = f64::MAX;
        let mut mz_max = f64::MIN;
        let mut rt_min = f64::MAX;
        let mut rt_max = f64::MIN;
        let mut charge_dist: HashMap<i32, u64> = HashMap::new();
        let mut peak_counts: Vec<u32> = Vec::new();

        parse_mzml_streaming(&mut xml_reader, path, |s| {
            total += 1;
            match s.ms_level {
                MsLevel::MS1 => ms1_count += 1,
                MsLevel::MS2 => ms2_count += 1,
                MsLevel::Other(_) => {}
            }

            if let Some(first) = s.mz_array.first() {
                if *first < mz_min {
                    mz_min = *first;
                }
            }
            if let Some(last) = s.mz_array.last() {
                if *last > mz_max {
                    mz_max = *last;
                }
            }
            if s.retention_time_sec < rt_min {
                rt_min = s.retention_time_sec;
            }
            if s.retention_time_sec > rt_max {
                rt_max = s.retention_time_sec;
            }
            for p in &s.precursors {
                if let Some(c) = p.charge {
                    *charge_dist.entry(c).or_insert(0) += 1;
                }
            }
            peak_counts.push(s.num_peaks() as u32);
            Ok(true)
        })?;

        if total == 0 {
            mz_min = 0.0;
            mz_max = 0.0;
            rt_min = 0.0;
            rt_max = 0.0;
        }

        peak_counts.sort_unstable();
        let median_peaks = if peak_counts.is_empty() {
            0
        } else {
            peak_counts[peak_counts.len() / 2]
        };

        let summary = SpectrumSummary {
            file_path: path.to_string_lossy().to_string(),
            format: SpectrumFormat::MzML,
            total_spectra: total,
            ms1_count,
            ms2_count,
            mz_range: (mz_min, mz_max),
            rt_range_sec: (rt_min, rt_max),
            precursor_charge_distribution: charge_dist,
            median_peaks_per_spectrum: median_peaks,
        };
        summary
            .validate()
            .map_err(|e| SpectrumIoError::ValidationError {
                scan: 0,
                detail: format!("summary: {e}"),
            })?;
        Ok(summary)
    }

    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        let file = std::fs::File::open(path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                SpectrumIoError::FileNotFound {
                    path: path.to_path_buf(),
                }
            } else {
                SpectrumIoError::IoError {
                    path: path.to_path_buf(),
                    source: e,
                }
            }
        })?;
        let buf_reader = std::io::BufReader::new(file);
        let mut xml_reader = Reader::from_reader(buf_reader);
        xml_reader.config_mut().trim_text(true);

        let mut found: Option<Spectrum> = None;
        parse_mzml_streaming(&mut xml_reader, path, |s| {
            if s.scan_number == scan {
                found = Some(s);
                Ok(false) // stop early
            } else {
                Ok(true)
            }
        })?;
        found.ok_or_else(|| SpectrumIoError::ScanNotFound {
            path: path.to_path_buf(),
            scan,
        })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("small.mzml")
    }

    // -- read_all -------------------------------------------------------

    #[test]
    fn read_all_parses_10_spectra() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        assert_eq!(spectra.len(), 10);
    }

    #[test]
    fn read_all_first_spectrum_correct() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let s = &spectra[0];

        assert_eq!(s.scan_number, 1);
        assert_eq!(s.ms_level, MsLevel::MS2);
        assert!((s.retention_time_sec - 120.5).abs() < 0.01);
        assert_eq!(s.num_peaks(), 5);
        assert_eq!(s.precursors.len(), 1);

        let p = &s.precursors[0];
        assert!((p.mz - 471.2561).abs() < 1e-4);
        assert_eq!(p.charge, Some(2));
        assert!((p.intensity.unwrap() - 1500000.0).abs() < 1.0);

        // Isolation window
        let iw = p.isolation_window.as_ref().unwrap();
        assert!((iw.target_mz - 471.2561).abs() < 1e-4);
        assert!((iw.lower_offset - 1.0).abs() < 0.01);
        assert!((iw.upper_offset - 1.0).abs() < 0.01);
    }

    #[test]
    fn read_all_binary_decode_correct() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let s = &spectra[0];

        // Verify decoded m/z values match fixture data
        assert!((s.mz_array[0] - 100.0510).abs() < 1e-3);
        assert!((s.mz_array[4] - 300.2100).abs() < 1e-3);
        assert!((s.intensity_array[0] - 1200.5).abs() < 0.1);
        assert!((s.intensity_array[3] - 15000.0).abs() < 0.1);
    }

    #[test]
    fn read_all_zlib_compressed_spectrum() {
        // Scan 9 uses zlib compression for its intensity array
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let s = &spectra[8]; // index 8 = scan 9

        assert_eq!(s.scan_number, 9);
        assert_eq!(s.num_peaks(), 8);
        assert!((s.intensity_array[4] - 9500.0).abs() < 0.1);
    }

    #[test]
    fn read_all_peak_counts_match_mgf() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let peaks: Vec<usize> = spectra.iter().map(|s| s.num_peaks()).collect();
        // Must match small.mgf fixture data
        assert_eq!(peaks, vec![5, 4, 6, 7, 3, 5, 6, 3, 8, 4]);
    }

    #[test]
    fn read_all_mz_arrays_sorted() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        for s in &spectra {
            for w in s.mz_array.windows(2) {
                assert!(w[0] <= w[1], "m/z not sorted in scan {}", s.scan_number);
            }
        }
    }

    #[test]
    fn read_all_charge_states() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        let charges: Vec<Option<i32>> = spectra
            .iter()
            .map(|s| s.precursors.first().and_then(|p| p.charge))
            .collect();
        assert_eq!(
            charges,
            vec![
                Some(2),
                Some(3),
                Some(2),
                Some(2),
                Some(3),
                Some(2),
                Some(3),
                Some(1),
                Some(4),
                Some(2),
            ]
        );
    }

    // -- read_summary ---------------------------------------------------

    #[test]
    fn read_summary_correct() {
        let reader = MzMLReader;
        let summary = reader.read_summary(&fixture_path()).unwrap();

        assert_eq!(summary.total_spectra, 10);
        assert_eq!(summary.ms2_count, 10);
        assert_eq!(summary.format, SpectrumFormat::MzML);
        assert!((summary.rt_range_sec.0 - 120.5).abs() < 0.1);
        assert!((summary.rt_range_sec.1 - 240.0).abs() < 0.1);

        let charge_2 = summary.precursor_charge_distribution.get(&2).unwrap_or(&0);
        assert_eq!(*charge_2, 5);
    }

    #[test]
    fn read_summary_validates() {
        let reader = MzMLReader;
        let summary = reader.read_summary(&fixture_path()).unwrap();
        assert!(summary.validate().is_ok());
    }

    // -- read_spectrum --------------------------------------------------

    #[test]
    fn read_spectrum_by_scan() {
        let reader = MzMLReader;
        let s = reader.read_spectrum(&fixture_path(), 7).unwrap();
        assert_eq!(s.scan_number, 7);
        assert_eq!(s.num_peaks(), 6);
        assert_eq!(s.precursors[0].charge, Some(3));
    }

    #[test]
    fn read_spectrum_not_found() {
        let reader = MzMLReader;
        let err = reader.read_spectrum(&fixture_path(), 999).unwrap_err();
        assert!(err.to_string().contains("999"));
    }

    // -- error cases ----------------------------------------------------

    #[test]
    fn read_all_file_not_found() {
        let reader = MzMLReader;
        let err = reader
            .read_all(Path::new("/nonexistent/file.mzml"))
            .unwrap_err();
        assert!(err.to_string().contains("not found"));
    }
}
