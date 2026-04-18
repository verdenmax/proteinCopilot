//! Indexed mzML reader with O(1) scan lookup via [`ScanIndex`].
//!
//! Uses byte offsets to seek directly to a `<spectrum>` node and parse
//! only the requested spectrum, avoiding a full file scan.

use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use protein_copilot_core::spectrum::{
    IsolationWindow, MsLevel, PrecursorInfo, Spectrum, SpectrumSummary,
};
use quick_xml::events::Event;
use quick_xml::Reader;

use crate::error::SpectrumIoError;
use crate::index::{build_index_by_byte_scan, ScanIndex};
use crate::mzml::{
    decode_binary_array, get_attr, parse_scan_from_id, parse_scan_from_spectrum_ref,
    BinaryArrayMeta, MzMLReader,
};
use crate::reader::SpectrumReader;

/// mzML reader backed by a [`ScanIndex`] for O(1) scan lookup.
///
/// Bulk operations (`read_all`, `read_summary`, `for_each_spectrum`) delegate
/// to the standard streaming [`MzMLReader`]; only `read_spectrum` uses the
/// seek-based path.
pub struct IndexedMzMLReader {
    index: ScanIndex,
    path: PathBuf,
}

impl IndexedMzMLReader {
    /// Opens an mzML file and builds a scan index.
    ///
    /// Uses a two-layer resolution strategy:
    /// 1. **Disk cache** (`.mzml.idx` sidecar, PCIX v2) — milliseconds, includes
    ///    full per-scan metadata (RT, ms_level, isolation_window)
    /// 2. **SIMD byte-scan** — seconds (scans full file with `memchr`), extracts
    ///    both offsets and metadata in one pass
    ///
    /// After layer 2 succeeds, the result is persisted to disk cache
    /// so future opens are instant.
    pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
        // Layer 1: Try disk cache
        if let Ok((file_size, file_mtime)) = crate::disk_cache::file_metadata(path) {
            match crate::disk_cache::load_index(path, file_size, file_mtime) {
                Ok(Some(cached_index)) => {
                    tracing::info!(
                        path = %path.display(),
                        scans = cached_index.len(),
                        "loaded index from disk cache"
                    );
                    return Ok(Self {
                        index: cached_index,
                        path: path.to_path_buf(),
                    });
                }
                Ok(None) => {} // Cache miss — continue to layer 2
                Err(e) => {
                    tracing::warn!(error = %e, "disk cache load error, continuing without cache");
                }
            }
        }

        // Layer 2: SIMD byte-scan (extracts full metadata: RT, ms_level,
        // isolation_window). We always prefer this over native <indexList>
        // because native index only provides byte offsets — no metadata —
        // which makes find_by_rt() unusable (all ms_level=0, all rt=0).
        let index = build_index_by_byte_scan(path)?;

        // Persist to disk cache for future opens (non-fatal if it fails)
        if let Ok((file_size, file_mtime)) = crate::disk_cache::file_metadata(path) {
            if let Err(e) = crate::disk_cache::save_index(path, &index, file_size, file_mtime) {
                tracing::warn!(error = %e, "failed to persist index cache (non-fatal)");
            }
        }

        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }

    /// Returns a reference to the underlying scan index.
    pub fn index(&self) -> &ScanIndex {
        &self.index
    }

    /// Seeks to `offset` in the file and parses the single `<spectrum>` node
    /// found there.
    fn read_spectrum_at_offset(&self, scan: u32, offset: u64) -> Result<Spectrum, SpectrumIoError> {
        let file = std::fs::File::open(&self.path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                SpectrumIoError::FileNotFound {
                    path: self.path.clone(),
                }
            } else {
                SpectrumIoError::IoError {
                    path: self.path.clone(),
                    source: e,
                }
            }
        })?;

        let mut buf_reader = BufReader::new(file);
        buf_reader
            .seek(SeekFrom::Start(offset))
            .map_err(|e| SpectrumIoError::IoError {
                path: self.path.clone(),
                source: e,
            })?;

        let mut xml_reader = Reader::from_reader(buf_reader);
        xml_reader.config_mut().trim_text(true);

        parse_single_spectrum(&mut xml_reader, &self.path, scan)
    }
}

impl SpectrumReader for IndexedMzMLReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        MzMLReader.read_all(path)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        MzMLReader.read_summary(path)
    }

    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        if path != self.path {
            let canonical_self =
                std::fs::canonicalize(&self.path).unwrap_or_else(|_| self.path.clone());
            let canonical_path = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
            if canonical_self != canonical_path {
                tracing::warn!(
                    "IndexedMzMLReader opened for {:?} but read_spectrum called with {:?}; using indexed file",
                    self.path, path
                );
            }
        }
        let offset = self
            .index
            .get_offset(scan)
            .ok_or_else(|| SpectrumIoError::ScanNotFound {
                path: self.path.clone(),
                scan,
            })?;
        self.read_spectrum_at_offset(scan, offset)
    }

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        MzMLReader.for_each_spectrum(path, handler)
    }

    fn find_by_rt(
        &self,
        _path: &Path,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Result<Option<(u32, f64)>, SpectrumIoError> {
        Ok(self.index.find_by_rt(rt_min, precursor_mz, rt_tolerance_min))
    }

    fn list_ms2_meta(
        &self,
        _path: &Path,
    ) -> Result<Vec<crate::reader::Ms2ScanMeta>, SpectrumIoError> {
        Ok(self
            .index
            .iter_meta()
            .filter(|(_, meta)| meta.ms_level == 2)
            .map(|(&scan, meta)| crate::reader::Ms2ScanMeta {
                scan_number: scan,
                rt_min: meta.rt_seconds / 60.0,
                isolation_window: meta.isolation_window,
            })
            .collect())
    }
}

// ---------------------------------------------------------------------------
// Single-spectrum parser
// ---------------------------------------------------------------------------

/// Parses exactly one `<spectrum>` node from the current XML reader position.
///
/// The reader must be positioned at (or just before) a `<spectrum>` start tag.
/// Returns the fully constructed [`Spectrum`] when `</spectrum>` is reached.
fn parse_single_spectrum<R: std::io::BufRead>(
    xml_reader: &mut Reader<R>,
    path: &Path,
    expected_scan: u32,
) -> Result<Spectrum, SpectrumIoError> {
    let mut buf = Vec::new();

    // State flags
    let mut in_spectrum = false;
    let mut in_precursor = false;
    let mut in_isolation_window = false;
    let mut in_selected_ion = false;
    let mut in_binary_data_array = false;
    let mut in_binary = false;
    let mut in_scan = false;

    // Spectrum fields
    let mut scan_number: Option<u32> = None;
    let mut ms_level: Option<u8> = None;
    let mut rt_min: Option<f64> = None;

    let mut precursors: Vec<PrecursorInfo> = Vec::new();
    let mut cur_precursor_mz: Option<f64> = None;
    let mut cur_precursor_charge: Option<i32> = None;
    let mut cur_precursor_intensity: Option<f64> = None;
    let mut cur_isolation_target_mz: Option<f64> = None;
    let mut cur_isolation_lower: Option<f64> = None;
    let mut cur_isolation_upper: Option<f64> = None;
    let mut cur_precursor_source_scan: Option<u32> = None;

    let mut mz_array: Vec<f64> = Vec::new();
    let mut intensity_array: Vec<f64> = Vec::new();
    let mut array_meta = BinaryArrayMeta::default();
    let mut binary_text = String::new();

    loop {
        match xml_reader.read_event_into(&mut buf) {
            Ok(Event::Eof) => {
                return Err(SpectrumIoError::XmlError {
                    path: path.to_path_buf(),
                    detail: format!(
                        "unexpected EOF while parsing spectrum at scan {expected_scan}"
                    ),
                });
            }
            Ok(Event::Start(ref e)) => {
                let local = e.local_name();
                let tag = local.as_ref();

                match tag {
                    b"spectrum" => {
                        in_spectrum = true;
                        if let Some(id) = get_attr(e, b"id") {
                            scan_number = parse_scan_from_id(&id);
                        }
                        if scan_number.is_none() {
                            scan_number = Some(expected_scan);
                        }
                    }
                    b"scan" if in_spectrum => {
                        in_scan = true;
                    }
                    b"precursor" if in_spectrum => {
                        in_precursor = true;
                        if let Some(spectrum_ref) = get_attr(e, b"spectrumRef") {
                            cur_precursor_source_scan = parse_scan_from_spectrum_ref(&spectrum_ref);
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
                            ms_level = value.parse().ok();
                        }
                        "MS:1000016" if in_scan => {
                            if let Ok(rt_val) = value.parse::<f64>() {
                                let unit_acc = get_attr(e, b"unitAccession").unwrap_or_default();
                                rt_min = Some(if unit_acc == "UO:0000031" {
                                    rt_val // already minutes
                                } else if unit_acc == "UO:0000010" {
                                    rt_val / 60.0 // seconds → minutes
                                } else {
                                    rt_val // assume minutes
                                });
                            }
                        }
                        "MS:1000827" if in_isolation_window => {
                            cur_isolation_target_mz = value.parse().ok();
                        }
                        "MS:1000828" if in_isolation_window => {
                            cur_isolation_lower = value.parse().ok();
                        }
                        "MS:1000829" if in_isolation_window => {
                            cur_isolation_upper = value.parse().ok();
                        }
                        "MS:1000744" if in_selected_ion => {
                            cur_precursor_mz = value.parse().ok();
                        }
                        "MS:1000041" if in_selected_ion => {
                            cur_precursor_charge = value.parse().ok();
                        }
                        "MS:1000042" if in_selected_ion => {
                            cur_precursor_intensity = value.parse().ok();
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
                        // Build and return the spectrum
                        let scan = scan_number.unwrap_or(expected_scan);
                        let level = match ms_level.unwrap_or(2) {
                            1 => MsLevel::MS1,
                            2 => MsLevel::MS2,
                            n => MsLevel::Other(n),
                        };

                        crate::util::sort_peaks_by_mz(&mut mz_array, &mut intensity_array);

                        return Spectrum::new(
                            scan,
                            level,
                            rt_min.unwrap_or(0.0),
                            precursors,
                            mz_array,
                            intensity_array,
                        )
                        .map_err(|e| SpectrumIoError::ValidationError {
                            scan,
                            detail: e.to_string(),
                        });
                    }
                    b"scan" => {
                        in_scan = false;
                    }
                    b"precursor" => {
                        in_precursor = false;
                        // Flush the current precursor
                        if let Some(mz) = cur_precursor_mz.take() {
                            let isolation_window = match (
                                cur_isolation_target_mz.take(),
                                cur_isolation_lower.take(),
                                cur_isolation_upper.take(),
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
                                charge: cur_precursor_charge.take(),
                                intensity: cur_precursor_intensity.take(),
                                isolation_window,
                                source_scan: cur_precursor_source_scan.take(),
                            });
                        } else {
                            cur_precursor_charge = None;
                            cur_precursor_intensity = None;
                            cur_isolation_target_mz = None;
                            cur_isolation_lower = None;
                            cur_isolation_upper = None;
                            cur_precursor_source_scan = None;
                        }
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
                                mz_array = decoded;
                            } else if array_meta.is_intensity {
                                intensity_array = decoded;
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::IndexSource;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml")
    }

    fn indexed_fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small_indexed.mzml")
    }

    #[test]
    fn open_plain_mzml_builds_scan_index() {
        // Remove any stale disk cache to test fresh index building
        let _ = std::fs::remove_file(crate::disk_cache::idx_path(&fixture_path()));
        let reader = IndexedMzMLReader::open(&fixture_path()).unwrap();
        assert_eq!(reader.index().len(), 10);
        assert_eq!(reader.index().source(), IndexSource::BuiltFromScan);
        // Clean up the .idx created as side effect
        let _ = std::fs::remove_file(crate::disk_cache::idx_path(&fixture_path()));
    }

    #[test]
    fn open_indexed_mzml_uses_byte_scan_with_metadata() {
        let path = indexed_fixture_path();
        if !path.exists() {
            return;
        }
        // Remove any stale disk cache to test fresh byte-scan path
        let _ = std::fs::remove_file(crate::disk_cache::idx_path(&path));
        let reader = IndexedMzMLReader::open(&path).unwrap();
        assert_eq!(reader.index().len(), 10);
        // Byte-scan is now always used (instead of native index) to ensure
        // full metadata (RT, ms_level, isolation_window) is available for find_by_rt.
        assert_eq!(reader.index().source(), IndexSource::BuiltFromScan);
        let _ = std::fs::remove_file(crate::disk_cache::idx_path(&path));
    }

    #[test]
    fn read_spectrum_by_index_returns_correct_scan() {
        let reader = IndexedMzMLReader::open(&fixture_path()).unwrap();
        let spec = reader.read_spectrum(&fixture_path(), 1).unwrap();
        assert_eq!(spec.scan_number, 1);
    }

    #[test]
    fn read_spectrum_by_index_scan_7_correct() {
        let reader = IndexedMzMLReader::open(&fixture_path()).unwrap();
        let spec = reader.read_spectrum(&fixture_path(), 7).unwrap();
        assert_eq!(spec.scan_number, 7);
        let standard = MzMLReader.read_spectrum(&fixture_path(), 7).unwrap();
        assert_eq!(spec.mz_array.len(), standard.mz_array.len());
        assert!((spec.retention_time_min - standard.retention_time_min).abs() < 0.01);
    }

    #[test]
    fn read_spectrum_not_found() {
        let reader = IndexedMzMLReader::open(&fixture_path()).unwrap();
        let err = reader.read_spectrum(&fixture_path(), 999).unwrap_err();
        assert!(err.to_string().contains("999"));
    }

    #[test]
    fn indexed_vs_standard_all_scans_match() {
        let indexed = IndexedMzMLReader::open(&fixture_path()).unwrap();
        for scan in 1..=10 {
            let idx_spec = indexed.read_spectrum(&fixture_path(), scan).unwrap();
            let std_spec = MzMLReader.read_spectrum(&fixture_path(), scan).unwrap();
            assert_eq!(idx_spec.scan_number, std_spec.scan_number);
            assert_eq!(idx_spec.mz_array.len(), std_spec.mz_array.len());
            assert_eq!(
                idx_spec.intensity_array.len(),
                std_spec.intensity_array.len()
            );
            assert!((idx_spec.retention_time_min - std_spec.retention_time_min).abs() < 0.001);
        }
    }

    #[test]
    fn open_creates_disk_cache() {
        let dir = tempfile::tempdir().unwrap();
        let src = fixture_path();
        let copy = dir.path().join("test.mzml");
        std::fs::copy(&src, &copy).unwrap();

        let idx_file = crate::disk_cache::idx_path(&copy);
        assert!(!idx_file.exists(), "idx should not exist before open");

        let _reader = IndexedMzMLReader::open(&copy).unwrap();
        assert!(idx_file.exists(), "idx should exist after open");
    }

    #[test]
    fn open_uses_disk_cache_on_second_call() {
        let dir = tempfile::tempdir().unwrap();
        let src = fixture_path();
        let copy = dir.path().join("test.mzml");
        std::fs::copy(&src, &copy).unwrap();

        // First open: builds index + saves cache
        let reader1 = IndexedMzMLReader::open(&copy).unwrap();

        // Second open: should load from disk cache
        let reader2 = IndexedMzMLReader::open(&copy).unwrap();
        assert_eq!(reader1.index().len(), reader2.index().len());

        // Verify scans match
        for scan in reader1.index().scan_numbers() {
            assert_eq!(
                reader1.index().get_offset(scan),
                reader2.index().get_offset(scan),
            );
        }
    }
}
