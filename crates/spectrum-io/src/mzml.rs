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

use std::io::BufRead;
use std::path::Path;

use protein_copilot_core::spectrum::{
    IsolationWindow, MsLevel, PrecursorInfo, Spectrum, SpectrumFormat,
    SpectrumRepresentation, SpectrumSummary,
};
use quick_xml::events::Event;
use quick_xml::Reader;

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Maximum peaks per spectrum. Prevents OOM from malformed mzML files.
/// 500k covers even wide DIA windows with high resolution.
const MAX_PEAKS_PER_SPECTRUM: usize = 500_000;

/// Reader for mzML spectrum files.
pub struct MzMLReader;

// ---------------------------------------------------------------------------
// Internal parsing types
// ---------------------------------------------------------------------------

#[derive(Default)]
pub(crate) struct BinaryArrayMeta {
    pub(crate) is_mz: bool,
    pub(crate) is_intensity: bool,
    pub(crate) is_64bit: bool,
    pub(crate) is_zlib: bool,
}

#[derive(Default)]
struct SpectrumBuilder {
    scan_number: Option<u32>,
    ms_level: Option<u8>,
    rt_min: Option<f64>,
    // Accumulated precursors (built at </precursor> end tag)
    precursors: Vec<PrecursorInfo>,
    // Temporary fields for the precursor currently being parsed
    cur_precursor_mz: Option<f64>,
    cur_precursor_charge: Option<i32>,
    cur_precursor_intensity: Option<f64>,
    cur_isolation_target_mz: Option<f64>,
    cur_isolation_lower: Option<f64>,
    cur_isolation_upper: Option<f64>,
    cur_precursor_source_scan: Option<u32>,
    mz_array: Vec<f64>,
    intensity_array: Vec<f64>,
    /// mzML CV term: MS:1000127 = centroid, MS:1000128 = profile
    representation: SpectrumRepresentation,
}

impl SpectrumBuilder {
    /// Build a `PrecursorInfo` from the current `cur_*` temporary fields,
    /// push it to `self.precursors`, then reset all temps.
    fn flush_precursor(&mut self) {
        if let Some(mz) = self.cur_precursor_mz.take() {
            let isolation_window = match (
                self.cur_isolation_target_mz.take(),
                self.cur_isolation_lower.take(),
                self.cur_isolation_upper.take(),
            ) {
                (Some(t), Some(l), Some(u)) => Some(IsolationWindow {
                    target_mz: t,
                    lower_offset: l,
                    upper_offset: u,
                }),
                _ => None,
            };
            self.precursors.push(PrecursorInfo {
                mz,
                charge: self.cur_precursor_charge.take(),
                intensity: self.cur_precursor_intensity.take(),
                isolation_window,
                source_scan: self.cur_precursor_source_scan.take(),
            });
        } else {
            // No m/z — discard partial precursor data
            self.cur_precursor_charge = None;
            self.cur_precursor_intensity = None;
            self.cur_isolation_target_mz = None;
            self.cur_isolation_lower = None;
            self.cur_isolation_upper = None;
            self.cur_precursor_source_scan = None;
        }
    }

    fn build(self, _path: &Path) -> Result<Spectrum, SpectrumIoError> {
        let scan = self.scan_number.unwrap_or(1);
        let ms_level = match self.ms_level.unwrap_or(2) {
            1 => MsLevel::MS1,
            2 => MsLevel::MS2,
            n => MsLevel::Other(n),
        };

        // Log when MS2+ spectrum has no retention time — RT-based lookups
        // (XIC, scan auto-matching) will use 0.0 as fallback.
        if self.rt_min.is_none() && !matches!(ms_level, MsLevel::MS1) {
            tracing::debug!(
                scan,
                "MS2 spectrum missing retention time, defaulting to 0.0"
            );
        }

        // Sort m/z + intensity arrays together by m/z ascending.
        let mut mz_array = self.mz_array;
        let mut intensity_array = self.intensity_array;
        crate::util::sort_peaks_by_mz(&mut mz_array, &mut intensity_array);

        // On-load centroiding: if the spectrum is explicitly profile-mode,
        // apply centroiding so downstream consumers always receive centroided
        // data.  Spectra already marked centroid or with unknown
        // representation are left as-is.
        if self.representation == SpectrumRepresentation::Profile {
            let (cent_mz, cent_int) = crate::util::centroid_spectrum(
                &mz_array,
                &intensity_array,
                1e-3,
            );
            if !cent_mz.is_empty() {
                mz_array = cent_mz;
                intensity_array = cent_int;
            }
        }

        let rep = if self.representation == SpectrumRepresentation::Profile {
            SpectrumRepresentation::Centroid
        } else {
            self.representation
        };
        Spectrum::new_with_rep(
            scan,
            ms_level,
            self.rt_min.unwrap_or(0.0),
            self.precursors,
            mz_array,
            intensity_array,
            rep,
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

pub(crate) fn decode_binary_array(
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
        // Guard against zlib bombs (zlib can amplify ≈1000×): bound the
        // decompression just past the largest legitimate array so a crafted
        // payload cannot exhaust memory before the per-spectrum peak cap below
        // is evaluated. The +1024 headroom keeps a max-size (MAX_PEAKS_PER_
        // SPECTRUM) array decoding fully without tripping the limit.
        let limit = (MAX_PEAKS_PER_SPECTRUM * 8 + 1024) as u64;
        let decoder = ZlibDecoder::new(&raw[..]);
        let mut decompressed = Vec::new();
        decoder
            .take(limit)
            .read_to_end(&mut decompressed)
            .map_err(|e| SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: format!("zlib decompress failed: {e}"),
            })?;
        if decompressed.len() as u64 == limit {
            return Err(SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: "decompressed binary array exceeds maximum size".to_string(),
            });
        }
        decompressed
    } else {
        raw
    };

    let values = if meta.is_64bit {
        if bytes.len() % 8 != 0 {
            return Err(SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: format!(
                    "64-bit array byte length {} not divisible by 8",
                    bytes.len()
                ),
            });
        }
        bytes
            .chunks_exact(8)
            .map(|c| {
                let arr: [u8; 8] = c.try_into().expect("chunks_exact(8) guarantees 8 bytes");
                f64::from_le_bytes(arr)
            })
            .collect::<Vec<f64>>()
    } else {
        if bytes.len() % 4 != 0 {
            return Err(SpectrumIoError::BinaryDecodeError {
                path: path.to_path_buf(),
                detail: format!(
                    "32-bit array byte length {} not divisible by 4",
                    bytes.len()
                ),
            });
        }
        bytes
            .chunks_exact(4)
            .map(|c| {
                let arr: [u8; 4] = c.try_into().expect("chunks_exact(4) guarantees 4 bytes");
                f32::from_le_bytes(arr) as f64
            })
            .collect::<Vec<f64>>()
    };

    if values.len() > MAX_PEAKS_PER_SPECTRUM {
        return Err(SpectrumIoError::BinaryDecodeError {
            path: path.to_path_buf(),
            detail: format!(
                "array has {} elements (max {MAX_PEAKS_PER_SPECTRUM}); file may be corrupt",
                values.len()
            ),
        });
    }

    Ok(values)
}

// ---------------------------------------------------------------------------
// XML attribute helpers
// ---------------------------------------------------------------------------

pub(crate) fn get_attr<'a>(
    e: &'a quick_xml::events::BytesStart<'a>,
    name: &[u8],
) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.as_ref() == name)
        .and_then(|a| String::from_utf8(a.value.to_vec()).ok())
}

/// Extract scan number from spectrum id attribute (e.g., "scan=123").
pub(crate) fn parse_scan_from_id(id: &str) -> Option<u32> {
    id.split("scan=")
        .nth(1)
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse().ok())
}

/// Extract scan number from a precursor `spectrumRef` attribute.
///
/// Typical format: `"controllerType=0 controllerNumber=1 scan=1234"`.
/// Falls back to parsing the whole string as a plain number.
pub(crate) fn parse_scan_from_spectrum_ref(spectrum_ref: &str) -> Option<u32> {
    if let Some(after) = spectrum_ref.split("scan=").nth(1) {
        // Take only leading digits after "scan="
        let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
        if !digits.is_empty() {
            return digits.parse().ok();
        }
    }
    // Fallback: try parsing the entire string as a plain number
    spectrum_ref.trim().parse().ok()
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
                        if let Some(spectrum_ref) = get_attr(e, b"spectrumRef") {
                            builder.cur_precursor_source_scan =
                                parse_scan_from_spectrum_ref(&spectrum_ref);
                        }
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
                            // Scan start time — check unit
                            // UO:0000010 = second, UO:0000031 = minute
                            // Internal convention: store in minutes
                            if let Ok(rt_val) = value.parse::<f64>() {
                                let unit_acc = get_attr(e, b"unitAccession").unwrap_or_default();
                                builder.rt_min = Some(if unit_acc == "UO:0000031" {
                                    rt_val // already minutes
                                } else if unit_acc == "UO:0000010" {
                                    rt_val / 60.0 // seconds → minutes
                                } else {
                                    if unit_acc.is_empty() {
                                        tracing::warn!(
                                            "MS:1000016 scan start time missing unitAccession; \
                                             assuming minutes (proteomics convention)"
                                        );
                                    }
                                    rt_val // assume minutes
                                });
                            }
                        }
                        "MS:1000827" if in_isolation_window => {
                            builder.cur_isolation_target_mz = value.parse().ok();
                        }
                        "MS:1000828" if in_isolation_window => {
                            builder.cur_isolation_lower = value.parse().ok();
                        }
                        "MS:1000829" if in_isolation_window => {
                            builder.cur_isolation_upper = value.parse().ok();
                        }
                        "MS:1000744" if in_selected_ion => {
                            builder.cur_precursor_mz = value.parse().ok();
                        }
                        "MS:1000041" if in_selected_ion => {
                            builder.cur_precursor_charge = value.parse().ok();
                        }
                        "MS:1000042" if in_selected_ion => {
                            builder.cur_precursor_intensity = value.parse().ok();
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
                        "MS:1000127" if in_spectrum => {
                            // centroid spectrum
                            builder.representation = SpectrumRepresentation::Centroid;
                        }
                        "MS:1000128" if in_spectrum => {
                            // profile spectrum
                            builder.representation = SpectrumRepresentation::Profile;
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
                        if count % 1000 == 0 {
                            tracing::info!(count, "streaming spectra");
                        }
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
                        builder.flush_precursor();
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

/// Opens a mzML file and creates a configured XML reader.
fn open_xml_reader(
    path: &Path,
) -> Result<Reader<std::io::BufReader<std::fs::File>>, SpectrumIoError> {
    let buf_reader = crate::util::open_buffered(path)?;
    let mut xml_reader = Reader::from_reader(buf_reader);
    xml_reader.config_mut().trim_text(true);
    Ok(xml_reader)
}

impl SpectrumReader for MzMLReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        let mut xml_reader = open_xml_reader(path)?;

        let mut spectra = Vec::new();
        parse_mzml_streaming(&mut xml_reader, path, |s| {
            spectra.push(s);
            Ok(true)
        })?;
        Ok(spectra)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        let mut xml_reader = open_xml_reader(path)?;

        let mut acc = crate::util::SummaryAccumulator::new();
        parse_mzml_streaming(&mut xml_reader, path, |s| {
            acc.observe(&s);
            Ok(true)
        })?;
        let summary = acc.into_summary(path, SpectrumFormat::MzML)?;
        tracing::info!(
            ms1 = summary.ms1_count,
            ms2 = summary.ms2_count,
            total = summary.total_spectra,
            "summary complete"
        );
        Ok(summary)
    }

    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        let mut xml_reader = open_xml_reader(path)?;

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

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        let mut xml_reader = open_xml_reader(path)?;
        parse_mzml_streaming(&mut xml_reader, path, handler)
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
        assert!((s.retention_time_min - 120.5 / 60.0).abs() < 0.01);
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
        assert!((summary.rt_range_min[0] - 120.5 / 60.0).abs() < 0.1);
        assert!((summary.rt_range_min[1] - 240.0 / 60.0).abs() < 0.1);

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

    // -- parse_scan_from_spectrum_ref -----------------------------------

    #[test]
    fn spectrum_ref_full_format() {
        assert_eq!(
            parse_scan_from_spectrum_ref("controllerType=0 controllerNumber=1 scan=1234"),
            Some(1234)
        );
    }

    #[test]
    fn spectrum_ref_scan_only() {
        assert_eq!(parse_scan_from_spectrum_ref("scan=42"), Some(42));
    }

    #[test]
    fn spectrum_ref_scan_after_other_key() {
        assert_eq!(
            parse_scan_from_spectrum_ref("spectrum=5 scan=100"),
            Some(100)
        );
    }

    #[test]
    fn spectrum_ref_empty() {
        assert_eq!(parse_scan_from_spectrum_ref(""), None);
    }

    #[test]
    fn spectrum_ref_no_scan() {
        assert_eq!(parse_scan_from_spectrum_ref("noscanhere"), None);
    }

    #[test]
    fn spectrum_ref_plain_number() {
        assert_eq!(parse_scan_from_spectrum_ref("5678"), Some(5678));
    }

    // -- multiple precursors per spectrum (regression) ---------------------

    /// Regression test: existing single-precursor fixture still produces
    /// exactly one precursor per spectrum with correct values.
    #[test]
    fn single_precursor_per_spectrum_preserved() {
        let reader = MzMLReader;
        let spectra = reader.read_all(&fixture_path()).unwrap();
        for s in &spectra {
            assert_eq!(
                s.precursors.len(),
                1,
                "scan {} should have exactly 1 precursor",
                s.scan_number
            );
        }
        // Spot-check first and last spectra
        let first = &spectra[0].precursors[0];
        assert!((first.mz - 471.2561).abs() < 1e-4);
        assert_eq!(first.charge, Some(2));

        let last = &spectra[9].precursors[0];
        assert!((last.mz - 445.23).abs() < 1e-4);
        assert_eq!(last.charge, Some(2));
    }

    /// Parsing an mzML with two <precursor> elements in one spectrum must
    /// produce a Spectrum with two entries in the `precursors` Vec.
    #[test]
    fn multiple_precursors_per_spectrum() {
        let mzml = r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="1" defaultDataProcessingRef="dp">
      <spectrum index="0" id="scan=1" defaultArrayLength="3">
        <cvParam cvRef="MS" accession="MS:1000511" value="2"/>
        <scanList count="1">
          <scan>
            <cvParam cvRef="MS" accession="MS:1000016" value="60.0" unitAccession="UO:0000010"/>
          </scan>
        </scanList>
        <precursorList count="2">
          <precursor spectrumRef="scan=100">
            <isolationWindow>
              <cvParam accession="MS:1000827" value="400.0"/>
              <cvParam accession="MS:1000828" value="0.5"/>
              <cvParam accession="MS:1000829" value="0.5"/>
            </isolationWindow>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam accession="MS:1000744" value="400.1234"/>
                <cvParam accession="MS:1000041" value="2"/>
                <cvParam accession="MS:1000042" value="50000.0"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
          <precursor spectrumRef="scan=101">
            <isolationWindow>
              <cvParam accession="MS:1000827" value="600.0"/>
              <cvParam accession="MS:1000828" value="1.0"/>
              <cvParam accession="MS:1000829" value="1.0"/>
            </isolationWindow>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam accession="MS:1000744" value="600.5678"/>
                <cvParam accession="MS:1000041" value="3"/>
                <cvParam accession="MS:1000042" value="75000.0"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
        </precursorList>
        <binaryDataArrayList count="2">
          <binaryDataArray>
            <cvParam accession="MS:1000514"/>
            <cvParam accession="MS:1000523"/>
            <cvParam accession="MS:1000576"/>
            <binary>AAAAAAAA+EAAAAAAAADwQAAAAAAAACRA</binary>
          </binaryDataArray>
          <binaryDataArray>
            <cvParam accession="MS:1000515"/>
            <cvParam accession="MS:1000523"/>
            <cvParam accession="MS:1000576"/>
            <binary>AAAAAAAA+EAAAAAAAADwQAAAAAAAACRA</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"#;

        let path = Path::new("multi_precursor.mzml");
        let mut xml_reader = Reader::from_str(mzml);
        xml_reader.config_mut().trim_text(true);

        let mut spectra = Vec::new();
        parse_mzml_streaming(&mut xml_reader, path, |s| {
            spectra.push(s);
            Ok(true)
        })
        .unwrap();

        assert_eq!(spectra.len(), 1);
        let s = &spectra[0];
        assert_eq!(s.precursors.len(), 2, "expected 2 precursors");

        let p0 = &s.precursors[0];
        assert!((p0.mz - 400.1234).abs() < 1e-4);
        assert_eq!(p0.charge, Some(2));
        assert!((p0.intensity.unwrap() - 50000.0).abs() < 1.0);
        assert_eq!(p0.source_scan, Some(100));
        let iw0 = p0.isolation_window.as_ref().unwrap();
        assert!((iw0.target_mz - 400.0).abs() < 1e-4);
        assert!((iw0.lower_offset - 0.5).abs() < 0.01);
        assert!((iw0.upper_offset - 0.5).abs() < 0.01);

        let p1 = &s.precursors[1];
        assert!((p1.mz - 600.5678).abs() < 1e-4);
        assert_eq!(p1.charge, Some(3));
        assert!((p1.intensity.unwrap() - 75000.0).abs() < 1.0);
        assert_eq!(p1.source_scan, Some(101));
        let iw1 = p1.isolation_window.as_ref().unwrap();
        assert!((iw1.target_mz - 600.0).abs() < 1e-4);
        assert!((iw1.lower_offset - 1.0).abs() < 0.01);
        assert!((iw1.upper_offset - 1.0).abs() < 0.01);
    }

    #[test]
    fn for_each_spectrum_streams_all() {
        let reader = MzMLReader;
        let path = fixture_path();
        let mut count = 0u32;
        let result = reader.for_each_spectrum(&path, &mut |_spec| {
            count += 1;
            Ok(true)
        });
        assert!(result.is_ok());
        let all = reader.read_all(&path).unwrap();
        assert_eq!(count, all.len() as u32);
    }

    #[test]
    fn for_each_spectrum_early_stop() {
        let reader = MzMLReader;
        let path = fixture_path();
        let mut count = 0u32;
        let _ = reader.for_each_spectrum(&path, &mut |_spec| {
            count += 1;
            Ok(count < 2)
        });
        assert_eq!(count, 2);
    }

    // -- zlib decompression bomb guard ----------------------------------

    #[test]
    fn decode_binary_array_rejects_zlib_bomb() {
        use base64::Engine;
        use flate2::write::ZlibEncoder;
        use flate2::Compression;
        use std::io::Write;

        // A payload that decompresses past the decoder limit
        // (MAX_PEAKS_PER_SPECTRUM*8 + 1024). All-zero bytes compress to a tiny
        // stream — a classic decompression bomb.
        let decompressed_size = (MAX_PEAKS_PER_SPECTRUM * 8 + 1024) + 4096;
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::best());
        encoder.write_all(&vec![0u8; decompressed_size]).unwrap();
        let compressed = encoder.finish().unwrap();
        // The compressed bomb is tiny — decoding it naïvely would balloon memory.
        assert!(
            compressed.len() < 100_000,
            "compressed bomb should be small, got {} bytes",
            compressed.len()
        );

        let b64 = base64::engine::general_purpose::STANDARD.encode(&compressed);
        let meta = BinaryArrayMeta {
            is_mz: true,
            is_intensity: false,
            is_64bit: true,
            is_zlib: true,
        };

        let err = decode_binary_array(&b64, &meta, Path::new("bomb.mzML")).unwrap_err();
        match err {
            SpectrumIoError::BinaryDecodeError { detail, .. } => assert!(
                detail.contains("exceeds maximum size"),
                "expected size-limit rejection, got: {detail}"
            ),
            other => panic!("expected BinaryDecodeError, got {other:?}"),
        }
    }

    #[test]
    fn decode_binary_array_accepts_legitimate_zlib() {
        use base64::Engine;
        use flate2::write::ZlibEncoder;
        use flate2::Compression;
        use std::io::Write;

        // A small, well-formed 64-bit array must still decode after the bomb
        // guard is in place (regression).
        let values: Vec<f64> = vec![100.0, 200.5, 300.25, 400.125];
        let mut raw = Vec::new();
        for v in &values {
            raw.extend_from_slice(&v.to_le_bytes());
        }
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&raw).unwrap();
        let compressed = encoder.finish().unwrap();
        let b64 = base64::engine::general_purpose::STANDARD.encode(&compressed);
        let meta = BinaryArrayMeta {
            is_mz: true,
            is_intensity: false,
            is_64bit: true,
            is_zlib: true,
        };
        let decoded = decode_binary_array(&b64, &meta, Path::new("ok.mzML")).unwrap();
        assert_eq!(decoded.len(), 4);
        assert!((decoded[1] - 200.5).abs() < 1e-9);
    }
}
