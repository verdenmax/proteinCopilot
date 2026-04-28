//! Scan index for O(1) spectrum lookup by scan number.
//!
//! Two construction paths:
//! 1. Parse native `<indexList>` from `<indexedmzML>` files (fast, reads only EOF)
//! 2. Build by scanning all `<spectrum>` tags and recording byte offsets (fallback)

use std::collections::HashMap;
use std::io::{BufRead, Read, Seek, SeekFrom};
use std::path::Path;

use crate::error::SpectrumIoError;

/// How the index was constructed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexSource {
    /// Parsed from the native `<indexList>` at the end of an `<indexedmzML>` file.
    NativeIndex,
    /// Built by scanning the file and recording `<spectrum>` byte offsets.
    BuiltFromScan,
}

/// Per-scan metadata stored in the index.
#[derive(Debug, Clone)]
pub struct ScanMeta {
    /// Byte offset of the `<spectrum>` opening tag in the mzML file.
    pub offset: u64,
    /// Retention time in seconds. 0.0 if unknown.
    pub rt_seconds: f64,
    /// MS level (1=MS1, 2=MS2). 0 if unknown.
    pub ms_level: u8,
    /// Isolation window: (target_mz, lower_offset, upper_offset). None for MS1 or unknown.
    pub isolation_window: Option<(f64, f64, f64)>,
}

/// Maps scan numbers to byte offsets and metadata within an mzML file.
///
/// Enables O(1) spectrum lookup by scan number and O(log N) lookup by
/// retention time via a pre-sorted RT index.
#[derive(Debug, Clone)]
pub struct ScanIndex {
    /// scan_number → metadata (offset, RT, ms_level, isolation_window).
    entries: HashMap<u32, ScanMeta>,
    /// How this index was built.
    source: IndexSource,
    /// Pre-sorted (rt_seconds, scan_number) pairs for binary search.
    rt_sorted: Vec<(f64, u32)>,
}

impl ScanIndex {
    /// Creates a ScanIndex from a legacy offset-only map.
    ///
    /// Metadata fields are set to defaults (rt=0, ms_level=0, no isolation).
    pub fn new(offsets: HashMap<u32, u64>, source: IndexSource) -> Self {
        let entries: HashMap<u32, ScanMeta> = offsets
            .into_iter()
            .map(|(scan, offset)| {
                (
                    scan,
                    ScanMeta {
                        offset,
                        rt_seconds: 0.0,
                        ms_level: 0,
                        isolation_window: None,
                    },
                )
            })
            .collect();
        let rt_sorted = build_rt_sorted(&entries);
        Self {
            entries,
            source,
            rt_sorted,
        }
    }

    /// Creates a ScanIndex from a full metadata map.
    pub fn from_meta(entries: HashMap<u32, ScanMeta>, source: IndexSource) -> Self {
        let rt_sorted = build_rt_sorted(&entries);
        Self {
            entries,
            source,
            rt_sorted,
        }
    }

    /// Returns the byte offset for a given scan number, or `None`.
    pub fn get_offset(&self, scan: u32) -> Option<u64> {
        self.entries.get(&scan).map(|m| m.offset)
    }

    /// Returns the full metadata for a given scan number, or `None`.
    pub fn get_meta(&self, scan: u32) -> Option<&ScanMeta> {
        self.entries.get(&scan)
    }

    /// Returns the number of indexed scans.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// How this index was constructed.
    pub fn source(&self) -> IndexSource {
        self.source
    }

    /// Returns offset references for disk cache serialization.
    pub fn iter_meta(&self) -> impl Iterator<Item = (&u32, &ScanMeta)> {
        self.entries.iter()
    }

    /// Returns a legacy offsets map (for backward compatibility).
    pub fn offsets(&self) -> HashMap<u32, u64> {
        self.entries
            .iter()
            .map(|(&scan, meta)| (scan, meta.offset))
            .collect()
    }

    /// Returns all indexed scan numbers, sorted ascending.
    pub fn scan_numbers(&self) -> Vec<u32> {
        let mut scans: Vec<u32> = self.entries.keys().copied().collect();
        scans.sort_unstable();
        scans
    }

    /// Returns the pre-sorted RT index for binary search.
    pub fn rt_sorted(&self) -> &[(f64, u32)] {
        &self.rt_sorted
    }

    /// Find the best MS2 scan matching a given RT and precursor m/z.
    ///
    /// Uses binary search on the pre-sorted RT index. O(log N + k) where
    /// k is the number of scans in the RT tolerance window.
    ///
    /// Returns `(scan_number, rt_delta_min)` or `None`.
    pub fn find_by_rt(
        &self,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Option<(u32, f64)> {
        let rt_sec = rt_min * 60.0;
        let tol_sec = rt_tolerance_min * 60.0;

        let start = self
            .rt_sorted
            .partition_point(|&(rt, _)| rt < rt_sec - tol_sec);

        let mut best: Option<(u32, f64)> = None;

        for &(rt, scan) in &self.rt_sorted[start..] {
            let delta_sec = rt - rt_sec;
            if delta_sec > tol_sec {
                break;
            }
            if delta_sec.abs() > tol_sec {
                continue;
            }

            let meta = match self.entries.get(&scan) {
                Some(m) => m,
                None => continue,
            };

            if meta.ms_level != 2 {
                continue;
            }

            if let Some((target, lower, upper)) = meta.isolation_window {
                let low = target - lower;
                let high = target + upper;
                if precursor_mz < low || precursor_mz > high {
                    continue;
                }
            }

            let delta_min = delta_sec / 60.0;
            match &best {
                None => best = Some((scan, delta_min)),
                Some((_, best_delta)) => {
                    if delta_min.abs() < best_delta.abs() {
                        best = Some((scan, delta_min));
                    }
                }
            }
        }

        best
    }
}

/// Build sorted (rt_seconds, scan_number) pairs from the entries map.
fn build_rt_sorted(entries: &HashMap<u32, ScanMeta>) -> Vec<(f64, u32)> {
    let mut sorted: Vec<(f64, u32)> = entries
        .iter()
        .map(|(&scan, meta)| (meta.rt_seconds, scan))
        .collect();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    sorted
}

/// Size of the tail chunk to read when searching for `<indexListOffset>`.
/// The indexListOffset element is always near the very end of the file.
const TAIL_READ_SIZE: usize = 4096;

/// Maximum allowed size (bytes) for the indexList XML region.
/// A legitimate `<indexList>` is typically a few hundred KB even for very
/// large files. 10 MB is a generous cap that prevents OOM from corrupted
/// `indexListOffset` values while still accommodating huge index lists.
const MAX_INDEX_READ_SIZE: u64 = 10 * 1024 * 1024;

/// Attempts to build a ScanIndex from the native `<indexList>` in an `<indexedmzML>` file.
///
/// Returns `Ok(Some(index))` if the file has a valid `<indexList>`,
/// `Ok(None)` if it's a plain `<mzML>` without an index,
/// `Err(...)` on I/O or parse errors.
pub fn build_index_from_native_mzml(path: &Path) -> Result<Option<ScanIndex>, SpectrumIoError> {
    use std::fs::File;

    let mut file = File::open(path).map_err(|e| {
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

    let file_len = file
        .metadata()
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?
        .len();

    if file_len == 0 {
        return Ok(None);
    }

    // Read the tail of the file to find <indexListOffset>
    let tail_start = file_len.saturating_sub(TAIL_READ_SIZE as u64);
    file.seek(SeekFrom::Start(tail_start))
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;

    let mut tail_bytes = Vec::new();
    file.read_to_end(&mut tail_bytes)
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;
    let tail = String::from_utf8_lossy(&tail_bytes);

    // Look for <indexListOffset>NNN</indexListOffset>
    let offset_str = match extract_between(&tail, "<indexListOffset>", "</indexListOffset>") {
        Some(s) => s,
        None => return Ok(None), // Not an indexedmzML file
    };

    let index_list_offset: u64 =
        offset_str
            .trim()
            .parse()
            .map_err(|_| SpectrumIoError::IndexParseError {
                path: path.to_path_buf(),
                detail: format!("invalid indexListOffset value: '{offset_str}'"),
            })?;

    // Sanity-check: offset must be within the file and the remaining
    // region must be reasonably small (prevents OOM on corrupted values).
    if index_list_offset >= file_len {
        return Err(SpectrumIoError::IndexParseError {
            path: path.to_path_buf(),
            detail: format!(
                "indexListOffset ({index_list_offset}) is beyond file size ({file_len})"
            ),
        });
    }
    let index_region_size = file_len - index_list_offset;
    if index_region_size > MAX_INDEX_READ_SIZE {
        return Err(SpectrumIoError::IndexParseError {
            path: path.to_path_buf(),
            detail: format!(
                "index region is too large ({index_region_size} bytes, max {MAX_INDEX_READ_SIZE}); \
                 indexListOffset may be corrupted"
            ),
        });
    }

    // Seek to the indexList and parse it
    file.seek(SeekFrom::Start(index_list_offset))
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;

    let mut index_xml = String::with_capacity(index_region_size as usize);
    file.read_to_string(&mut index_xml)
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;

    let offsets = parse_index_list(&index_xml, path)?;
    if offsets.is_empty() {
        return Ok(None);
    }

    Ok(Some(ScanIndex::new(offsets, IndexSource::NativeIndex)))
}

/// Extracts the text between two delimiters in a string.
fn extract_between<'a>(text: &'a str, start: &str, end: &str) -> Option<&'a str> {
    let start_pos = text.find(start)? + start.len();
    let end_pos = text[start_pos..].find(end)? + start_pos;
    Some(&text[start_pos..end_pos])
}

/// Parses `<indexList>` XML to extract spectrum offsets.
fn parse_index_list(xml: &str, path: &Path) -> Result<HashMap<u32, u64>, SpectrumIoError> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    let mut reader = Reader::from_str(xml);
    reader.config_mut().trim_text(true);

    let mut offsets = HashMap::new();
    let mut in_spectrum_index = false;
    let mut in_offset = false;
    let mut current_id_ref = String::new();
    let mut fallback_scan: u32 = 0;
    let mut buf = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Eof) => break,
            Ok(Event::Start(ref e)) => {
                let tag = e.local_name();
                match tag.as_ref() {
                    b"index" => {
                        if let Some(name) = e
                            .attributes()
                            .filter_map(|a| a.ok())
                            .find(|a| a.key.as_ref() == b"name")
                            .and_then(|a| String::from_utf8(a.value.to_vec()).ok())
                        {
                            in_spectrum_index = name == "spectrum";
                        }
                    }
                    b"offset" if in_spectrum_index => {
                        in_offset = true;
                        current_id_ref = e
                            .attributes()
                            .filter_map(|a| a.ok())
                            .find(|a| a.key.as_ref() == b"idRef")
                            .and_then(|a| String::from_utf8(a.value.to_vec()).ok())
                            .unwrap_or_default();
                        fallback_scan += 1;
                    }
                    _ => {}
                }
            }
            Ok(Event::Text(ref t)) if in_offset => {
                let text = t.unescape().map_err(|e| SpectrumIoError::IndexParseError {
                    path: path.to_path_buf(),
                    detail: format!("XML unescape error in indexList: {e}"),
                })?;
                if let Ok(byte_offset) = text.trim().parse::<u64>() {
                    let scan = parse_scan_from_id_ref(&current_id_ref).unwrap_or(fallback_scan);
                    if let Some(prev_offset) = offsets.insert(scan, byte_offset) {
                        tracing::warn!(
                            "duplicate scan {} in index: offset {} replaced by {}",
                            scan,
                            prev_offset,
                            byte_offset
                        );
                    }
                }
            }
            Ok(Event::End(ref e)) => match e.local_name().as_ref() {
                b"offset" => in_offset = false,
                b"index" => in_spectrum_index = false,
                b"indexList" => break,
                _ => {}
            },
            Err(e) => {
                return Err(SpectrumIoError::IndexParseError {
                    path: path.to_path_buf(),
                    detail: format!("XML error parsing indexList: {e}"),
                });
            }
            _ => {}
        }
        buf.clear();
    }

    Ok(offsets)
}

/// Extracts scan number from an idRef like "scan=123" or
/// "controllerType=0 controllerNumber=1 scan=123".
fn parse_scan_from_id_ref(id_ref: &str) -> Option<u32> {
    id_ref.split("scan=").nth(1).and_then(|s| {
        let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
        digits.parse().ok()
    })
}

/// Builds a ScanIndex by scanning the file and recording byte offsets of
/// `<spectrum` opening tags. Used as fallback when native index is absent.
///
/// This reads the file as raw bytes (not XML) for speed, looking for
/// `<spectrum ` or `<spectrum>` tag starts, then extracts scan numbers
/// from the `id` attribute.
pub fn build_index_by_scanning(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    use std::fs::File;
    use std::io::BufReader;

    let file = File::open(path).map_err(|e| {
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

    let mut reader = BufReader::new(file);
    let mut offsets = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut byte_pos: u64 = 0;

    let mut line = String::new();
    loop {
        let line_start = byte_pos;
        line.clear();
        let bytes_read = reader
            .read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;
        if bytes_read == 0 {
            break;
        }
        byte_pos += bytes_read as u64;

        let trimmed = line.trim();

        if trimmed.starts_with("<spectrum ") || trimmed.starts_with("<spectrum>") {
            fallback_scan += 1;
            let scan = extract_id_attr(trimmed)
                .and_then(|id| parse_scan_from_id_ref(&id))
                .unwrap_or(fallback_scan);
            if let Some(prev_offset) = offsets.insert(scan, line_start) {
                tracing::warn!(
                    "duplicate scan {} found while scanning: offset {} replaced by {}",
                    scan,
                    prev_offset,
                    line_start
                );
            }
        }
    }

    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
}

/// Extracts a scan number from raw `<spectrum …>` tag bytes without UTF-8 validation.
///
/// Searches for ` id="` (note leading space to avoid matching `nativeID=`),
/// extracts the attribute value, then parses `scan=N` from it.
/// Falls back to single-quote variant ` id='`.
/// At most `512` bytes are inspected to avoid scanning into peak data.
/// Returns `fallback_scan` if no id attribute or no `scan=` is found.
fn extract_scan_from_tag_bytes(tag_bytes: &[u8], fallback_scan: u32) -> u32 {
    let limit = tag_bytes.len().min(512);
    let region = &tag_bytes[..limit];

    // Try double-quote first, then single-quote
    let (after_id, closing_quote) = if let Some(pos) = memchr::memmem::find(region, b" id=\"") {
        (&region[pos + 5..], b'"')
    } else if let Some(pos) = memchr::memmem::find(region, b" id='") {
        (&region[pos + 5..], b'\'')
    } else {
        return fallback_scan;
    };

    // Find the closing quote to delimit the attribute value
    let end = match memchr::memchr(closing_quote, after_id) {
        Some(e) => e,
        None => return fallback_scan,
    };
    let id_value = &after_id[..end];

    // Look for "scan=" inside the id value
    if let Some(scan_pos) = memchr::memmem::find(id_value, b"scan=") {
        let digits_start = scan_pos + 5;
        let mut digits_end = digits_start;
        while digits_end < id_value.len() && id_value[digits_end].is_ascii_digit() {
            digits_end += 1;
        }
        if digits_end > digits_start {
            // SAFETY: we checked that all bytes are ASCII digits, so this is valid UTF-8
            if let Ok(s) = std::str::from_utf8(&id_value[digits_start..digits_end]) {
                if let Ok(n) = s.parse::<u32>() {
                    return n;
                }
            }
        }
    }

    fallback_scan
}

/// Extracts the `value="..."` attribute from a cvParam byte region.
///
/// Searches for `value="` within the first 200 bytes, parses the f64.
fn extract_cv_value(region: &[u8]) -> Option<f64> {
    let limit = region.len().min(200);
    let search = &region[..limit];
    let pos = memchr::memmem::find(search, b"value=\"")?;
    let after = &search[pos + 7..];
    let end = memchr::memchr(b'"', after)?;
    let val_bytes = &after[..end];
    std::str::from_utf8(val_bytes)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
}

/// Extracts RT, ms_level, and isolation window from raw XML bytes
/// following a `<spectrum ` tag.
///
/// Searches for well-known cvParam accession numbers in the raw bytes.
/// Stops at `<binaryDataArrayList` to avoid scanning into peak data.
fn extract_meta_from_region(region: &[u8]) -> (f64, u8, Option<(f64, f64, f64)>) {
    let mut rt_seconds: f64 = 0.0;
    let mut ms_level: u8 = 0;
    let mut iso_target: Option<f64> = None;
    let mut iso_lower: Option<f64> = None;
    let mut iso_upper: Option<f64> = None;

    // Stop scanning if we hit binaryDataArrayList (no metadata after this)
    let limit = memchr::memmem::find(region, b"<binaryDataArrayList").unwrap_or(region.len());
    let region = &region[..limit];

    // MS:1000016 — scan start time (RT)
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000016") {
        rt_seconds = extract_cv_value(&region[pos..]).unwrap_or(0.0);
        // Check if unit is minutes (UO:0000031) — convert to seconds
        let end = region.len().min(pos + 300);
        let after = &region[pos..end];
        if memchr::memmem::find(after, b"UO:0000031").is_some() {
            rt_seconds *= 60.0;
        }
    }

    // MS:1000511 — ms level (the value attribute contains the level number)
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000511") {
        ms_level = extract_cv_value(&region[pos..])
            .map(|v| v as u8)
            .unwrap_or(0);
    }

    // MS:1000827 — isolation window target m/z
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000827") {
        iso_target = extract_cv_value(&region[pos..]);
    }
    // MS:1000828 — isolation window lower offset
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000828") {
        iso_lower = extract_cv_value(&region[pos..]);
    }
    // MS:1000829 — isolation window upper offset
    if let Some(pos) = memchr::memmem::find(region, b"MS:1000829") {
        iso_upper = extract_cv_value(&region[pos..]);
    }

    let isolation_window = match (iso_target, iso_lower, iso_upper) {
        (Some(t), Some(l), Some(u)) => Some((t, l, u)),
        _ => None,
    };

    (rt_seconds, ms_level, isolation_window)
}

/// Builds a [`ScanIndex`] by byte-level scanning with SIMD-accelerated search.
///
/// Uses `memchr::memmem` to find `<spectrum ` needles in large buffered reads,
/// avoiding per-line `String` allocation and UTF-8 validation. This is expected
/// to be 5–10× faster than [`build_index_by_scanning`] on multi-GB mzML files.
///
/// Cross-buffer-boundary matches are handled by keeping `needle.len() - 1`
/// bytes of overlap between consecutive buffer fills.
pub fn build_index_by_byte_scan(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    use std::fs::File;
    use std::io::BufReader;

    let file = File::open(path).map_err(|e| {
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

    const BUF_CAPACITY: usize = 256 * 1024;
    let mut reader = BufReader::with_capacity(BUF_CAPACITY, file);
    let needle = b"<spectrum ";
    let mut entries: HashMap<u32, ScanMeta> = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut global_pos: u64 = 0;

    // We need enough bytes after a `<spectrum ` match to extract metadata up to
    // the `<binaryDataArrayList` tag.  In real DIA files the isolation-window
    // cvParams (MS:1000827/28/29) can appear 3000+ bytes from the tag start
    // (e.g. Thermo DIA mzML with long scanWindow / userParam blocks).
    // 8192 bytes safely covers all observed layouts.
    const TAG_MIN_CONTENT: usize = 8192;

    loop {
        let buf = reader.fill_buf().map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;
        if buf.is_empty() {
            break;
        }
        let buf_len = buf.len();
        let mut search_start = 0;

        while let Some(pos) = memchr::memmem::find(&buf[search_start..], needle) {
            let local_pos = search_start + pos;
            let remaining = buf_len - local_pos;

            // Only skip near-boundary tags when the buffer is at full capacity,
            // meaning more data may follow. When `buf_len < BUF_CAPACITY` we are
            // at the tail of the file — process with whatever bytes remain.
            if remaining < TAG_MIN_CONTENT && buf_len >= BUF_CAPACITY {
                break;
            }

            let abs_pos = global_pos + local_pos as u64;
            fallback_scan += 1;
            // Extract scan from tag bytes (limit to 512 bytes or end of buffer)
            let tag_end = (local_pos + 512).min(buf_len);
            let scan = extract_scan_from_tag_bytes(&buf[local_pos..tag_end], fallback_scan);

            let meta_end = (local_pos + TAG_MIN_CONTENT).min(buf_len);
            let (rt_seconds, ms_level, isolation_window) =
                extract_meta_from_region(&buf[local_pos..meta_end]);

            let meta = ScanMeta {
                offset: abs_pos,
                rt_seconds,
                ms_level,
                isolation_window,
            };

            if let Some(prev) = entries.insert(scan, meta) {
                tracing::warn!(
                    "duplicate scan {} found while byte-scanning: offset {} replaced by {}",
                    scan,
                    prev.offset,
                    abs_pos
                );
            }
            search_start = local_pos + needle.len();
        }

        // Keep overlap large enough to cover tags near the buffer boundary.
        // TAG_MIN_CONTENT bytes ensures any skipped tag is fully available
        // in the next fill.
        let overlap = TAG_MIN_CONTENT + needle.len();
        let consumed = if buf_len > overlap {
            buf_len - overlap
        } else {
            buf_len
        };
        global_pos += consumed as u64;
        reader.consume(consumed);
    }

    Ok(ScanIndex::from_meta(entries, IndexSource::BuiltFromScan))
}

/// Extracts the value of the `id` attribute from an XML tag string.
fn extract_id_attr(tag_text: &str) -> Option<String> {
    // Match " id=" with leading space to avoid suffix matches (e.g., "nativeID=")
    let search_dq = " id=\"";
    let search_sq = " id='";
    let after_id = if let Some(pos) = tag_text.find(search_dq) {
        &tag_text[pos + search_dq.len()..]
    } else if let Some(pos) = tag_text.find(search_sq) {
        &tag_text[pos + search_sq.len()..]
    } else {
        return None;
    };
    let end = after_id.find('"').or_else(|| after_id.find('\''))?;
    Some(after_id[..end].to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_index_basic_operations() {
        let mut map = HashMap::new();
        map.insert(1, 100);
        map.insert(5, 500);
        map.insert(10, 1000);
        let idx = ScanIndex::new(map, IndexSource::NativeIndex);

        assert_eq!(idx.len(), 3);
        assert!(!idx.is_empty());
        assert_eq!(idx.get_offset(1), Some(100));
        assert_eq!(idx.get_offset(5), Some(500));
        assert_eq!(idx.get_offset(10), Some(1000));
        assert_eq!(idx.get_offset(99), None);
        assert_eq!(idx.source(), IndexSource::NativeIndex);
        assert_eq!(idx.scan_numbers(), vec![1, 5, 10]);
    }

    #[test]
    fn scan_index_empty() {
        let idx = ScanIndex::new(HashMap::new(), IndexSource::BuiltFromScan);
        assert!(idx.is_empty());
        assert_eq!(idx.len(), 0);
        assert_eq!(idx.scan_numbers(), Vec::<u32>::new());
    }

    #[test]
    fn parse_scan_from_id_ref_standard() {
        assert_eq!(parse_scan_from_id_ref("scan=123"), Some(123));
    }

    #[test]
    fn parse_scan_from_id_ref_with_controller() {
        assert_eq!(
            parse_scan_from_id_ref("controllerType=0 controllerNumber=1 scan=456"),
            Some(456)
        );
    }

    #[test]
    fn parse_scan_from_id_ref_no_scan() {
        assert_eq!(parse_scan_from_id_ref("spectrum_123"), None);
    }

    #[test]
    fn extract_id_attr_double_quotes() {
        assert_eq!(
            extract_id_attr(r#"<spectrum index="0" id="scan=1" defaultArrayLength="4">"#),
            Some("scan=1".to_string())
        );
    }

    #[test]
    fn extract_id_attr_single_quotes() {
        assert_eq!(
            extract_id_attr("<spectrum index='0' id='scan=1'>"),
            Some("scan=1".to_string())
        );
    }

    #[test]
    fn extract_between_works() {
        let text = "prefix<indexListOffset>12345</indexListOffset>suffix";
        assert_eq!(
            extract_between(text, "<indexListOffset>", "</indexListOffset>"),
            Some("12345")
        );
    }

    #[test]
    fn extract_between_missing() {
        assert_eq!(extract_between("no tag here", "<a>", "</a>"), None);
    }

    #[test]
    fn build_index_from_native_plain_mzml_returns_none() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let result = build_index_from_native_mzml(&path).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn build_index_by_scanning_finds_spectra() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let idx = build_index_by_scanning(&path).unwrap();
        assert_eq!(idx.len(), 10); // small.mzml has 10 spectra
        assert_eq!(idx.source(), IndexSource::BuiltFromScan);
        assert!(idx.get_offset(1).is_some());
        assert!(idx.get_offset(10).is_some());
    }

    #[test]
    fn build_index_from_native_indexed_mzml() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small_indexed.mzml");
        if !path.exists() {
            generate_indexed_fixture(&path);
        }
        let result = build_index_from_native_mzml(&path).unwrap();
        let idx = result.expect("should find native index");
        assert_eq!(idx.len(), 10);
        assert_eq!(idx.source(), IndexSource::NativeIndex);
        for scan in 1..=10 {
            assert!(idx.get_offset(scan).is_some(), "missing scan {scan}");
        }
    }

    // ── byte-scan & extract_scan_from_tag_bytes tests ──────────────────

    #[test]
    fn byte_scan_matches_line_scan() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let line_idx = build_index_by_scanning(&path).unwrap();
        let byte_idx = build_index_by_byte_scan(&path).unwrap();

        assert_eq!(
            line_idx.len(),
            byte_idx.len(),
            "number of indexed scans must match"
        );
        assert_eq!(byte_idx.source(), IndexSource::BuiltFromScan);

        // Both must find the same scan numbers
        assert_eq!(line_idx.scan_numbers(), byte_idx.scan_numbers());

        // The byte scanner finds the exact `<spectrum ` position while
        // the line scanner records the start of the line (including
        // leading whitespace). The byte scanner offset should be ≥ the
        // line offset and within a small indentation delta.
        let file_bytes = std::fs::read(&path).unwrap();
        for scan in line_idx.scan_numbers() {
            let byte_off = byte_idx
                .get_offset(scan)
                .unwrap_or_else(|| panic!("byte_scan missing scan {scan}"));
            let line_off = line_idx.get_offset(scan).unwrap();
            assert!(
                byte_off >= line_off,
                "byte offset {byte_off} should be >= line offset {line_off} for scan {scan}"
            );
            assert!(
                byte_off - line_off < 64,
                "byte offset delta too large for scan {scan}: line={line_off}, byte={byte_off}"
            );
            // Verify the byte offset actually points at `<spectrum `
            assert_eq!(
                &file_bytes[byte_off as usize..byte_off as usize + 10],
                b"<spectrum ",
                "byte offset for scan {scan} should point at '<spectrum '"
            );
        }
    }

    #[test]
    fn byte_scan_indexed_mzml_matches() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small_indexed.mzml");
        if !path.exists() {
            generate_indexed_fixture(&path);
        }
        let idx = build_index_by_byte_scan(&path).unwrap();
        assert_eq!(idx.len(), 10, "small_indexed.mzml should have 10 spectra");
        assert_eq!(idx.source(), IndexSource::BuiltFromScan);
    }

    #[test]
    fn extract_scan_from_tag_bytes_standard() {
        let tag = b"<spectrum index=\"0\" id=\"scan=42\" defaultArrayLength=\"4\">";
        assert_eq!(extract_scan_from_tag_bytes(tag, 99), 42);
    }

    #[test]
    fn extract_scan_from_tag_bytes_with_controller() {
        let tag = b"<spectrum id=\"controllerType=0 controllerNumber=1 scan=123\">";
        assert_eq!(extract_scan_from_tag_bytes(tag, 99), 123);
    }

    #[test]
    fn extract_scan_from_tag_bytes_no_id() {
        let tag = b"<spectrum index=\"0\">";
        assert_eq!(extract_scan_from_tag_bytes(tag, 99), 99);
    }

    fn generate_indexed_fixture(output_path: &std::path::Path) {
        use std::io::Write;

        let source_path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let source = std::fs::read_to_string(&source_path).unwrap();

        let mzml_content = source
            .strip_prefix("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
            .unwrap_or(&source);

        let header = "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\">\n";

        let body = format!("{}{}", header, mzml_content);

        let mut offsets: Vec<(u32, usize)> = Vec::new();
        let mut search_start = 0;
        let mut fallback_scan = 0u32;
        while let Some(pos) = body[search_start..].find("<spectrum ") {
            let abs_pos = search_start + pos;
            fallback_scan += 1;
            let tag_end = body[abs_pos..].find('>').unwrap_or(200) + abs_pos;
            let tag_text = &body[abs_pos..tag_end];
            let scan = extract_id_attr(tag_text)
                .and_then(|id| parse_scan_from_id_ref(&id))
                .unwrap_or(fallback_scan);
            offsets.push((scan, abs_pos));
            search_start = abs_pos + 1;
        }

        let mut index_entries = String::new();
        for (scan, offset) in &offsets {
            index_entries.push_str(&format!(
                "      <offset idRef=\"scan={scan}\">{offset}</offset>\n"
            ));
        }

        let index_list_offset = body.len();
        let index_list = format!(
            "  <indexList count=\"{}\">\n    <index name=\"spectrum\">\n{}    </index>\n  </indexList>\n",
            offsets.len(),
            index_entries,
        );

        let footer = format!(
            "{}  <indexListOffset>{}</indexListOffset>\n</indexedmzML>\n",
            index_list, index_list_offset,
        );

        let mut out = std::fs::File::create(output_path).unwrap();
        write!(out, "{}{}", body, footer).unwrap();
    }

    fn make_rt_index() -> ScanIndex {
        let mut meta = HashMap::new();
        // MS1 scans (should be skipped by find_by_rt)
        meta.insert(
            1,
            ScanMeta {
                offset: 100,
                rt_seconds: 100.0 * 60.0,
                ms_level: 1,
                isolation_window: None,
            },
        );
        meta.insert(
            3,
            ScanMeta {
                offset: 300,
                rt_seconds: 200.0 * 60.0,
                ms_level: 1,
                isolation_window: None,
            },
        );
        // MS2 scans
        meta.insert(
            2,
            ScanMeta {
                offset: 200,
                rt_seconds: 100.0 * 60.0,
                ms_level: 2,
                isolation_window: Some((500.0, 12.5, 12.5)),
            },
        );
        meta.insert(
            4,
            ScanMeta {
                offset: 400,
                rt_seconds: 200.0 * 60.0,
                ms_level: 2,
                isolation_window: Some((600.0, 12.5, 12.5)),
            },
        );
        meta.insert(
            5,
            ScanMeta {
                offset: 500,
                rt_seconds: 300.0 * 60.0,
                ms_level: 2,
                isolation_window: Some((500.0, 25.0, 25.0)),
            },
        );
        meta.insert(
            6,
            ScanMeta {
                offset: 600,
                rt_seconds: 400.0 * 60.0,
                ms_level: 2,
                isolation_window: None,
            },
        );
        ScanIndex::from_meta(meta, IndexSource::BuiltFromScan)
    }

    #[test]
    fn find_by_rt_exact_match() {
        let idx = make_rt_index();
        let result = idx.find_by_rt(100.0, 500.0, 30.0);
        assert_eq!(result.unwrap().0, 2);
    }

    #[test]
    fn find_by_rt_skips_ms1() {
        let idx = make_rt_index();
        let result = idx.find_by_rt(100.0, 500.0, 30.0);
        assert_eq!(result.unwrap().0, 2);
    }

    #[test]
    fn find_by_rt_mz_outside_window() {
        let idx = make_rt_index();
        let result = idx.find_by_rt(100.0, 550.0, 30.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_by_rt_outside_tolerance() {
        let idx = make_rt_index();
        let result = idx.find_by_rt(150.0, 500.0, 30.0);
        assert!(result.is_none());
    }

    #[test]
    fn find_by_rt_dda_no_isolation_accepts_any_mz() {
        let idx = make_rt_index();
        let result = idx.find_by_rt(400.0, 999.0, 30.0);
        assert_eq!(result.unwrap().0, 6);
    }

    #[test]
    fn find_by_rt_picks_closest() {
        let mut meta = HashMap::new();
        meta.insert(
            1,
            ScanMeta {
                offset: 100,
                rt_seconds: 100.0 * 60.0,
                ms_level: 2,
                isolation_window: Some((500.0, 25.0, 25.0)),
            },
        );
        meta.insert(
            2,
            ScanMeta {
                offset: 200,
                rt_seconds: 105.0 * 60.0,
                ms_level: 2,
                isolation_window: Some((500.0, 25.0, 25.0)),
            },
        );
        let idx = ScanIndex::from_meta(meta, IndexSource::BuiltFromScan);
        let result = idx.find_by_rt(103.0, 500.0, 30.0);
        assert_eq!(result.unwrap().0, 2);
    }

    #[test]
    fn find_by_rt_empty_index() {
        let idx = ScanIndex::from_meta(HashMap::new(), IndexSource::BuiltFromScan);
        assert!(idx.find_by_rt(100.0, 500.0, 30.0).is_none());
    }

    #[test]
    fn scan_index_with_meta_basic() {
        let mut meta_map = HashMap::new();
        meta_map.insert(
            1,
            ScanMeta {
                offset: 100,
                rt_seconds: 120.5,
                ms_level: 2,
                isolation_window: Some((500.0, 1.0, 1.0)),
            },
        );
        meta_map.insert(
            5,
            ScanMeta {
                offset: 5000,
                rt_seconds: 300.0,
                ms_level: 1,
                isolation_window: None,
            },
        );
        let idx = ScanIndex::from_meta(meta_map, IndexSource::NativeIndex);

        assert_eq!(idx.len(), 2);
        assert_eq!(idx.get_offset(1), Some(100));
        assert_eq!(idx.get_offset(5), Some(5000));
        assert_eq!(idx.get_offset(99), None);

        let meta = idx.get_meta(1).unwrap();
        assert_eq!(meta.ms_level, 2);
        assert!((meta.rt_seconds - 120.5).abs() < 0.001);
        assert!(meta.isolation_window.is_some());
        assert!(idx.get_meta(99).is_none());
    }

    #[test]
    fn scan_index_from_meta_backward_compat() {
        let mut offsets = HashMap::new();
        offsets.insert(1u32, 100u64);
        offsets.insert(2, 200);
        let idx = ScanIndex::new(offsets, IndexSource::BuiltFromScan);
        assert_eq!(idx.len(), 2);
        assert_eq!(idx.get_offset(1), Some(100));
        let meta = idx.get_meta(1).unwrap();
        assert_eq!(meta.offset, 100);
        assert_eq!(meta.ms_level, 0);
        assert!((meta.rt_seconds).abs() < 0.001);
    }

    #[test]
    fn byte_scan_extracts_metadata() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let idx = build_index_by_byte_scan(&path).unwrap();

        assert_eq!(idx.len(), 10);

        // Scan 1: RT=120.5s, ms_level=2, isolation=(471.2561, 1.0, 1.0)
        let meta1 = idx.get_meta(1).expect("scan 1 should exist");
        assert_eq!(meta1.ms_level, 2, "scan 1 ms_level");
        assert!(
            (meta1.rt_seconds - 120.5).abs() < 0.1,
            "scan 1 RT: expected ~120.5, got {}",
            meta1.rt_seconds
        );
        let iw = meta1
            .isolation_window
            .expect("scan 1 should have isolation window");
        assert!(
            (iw.0 - 471.2561).abs() < 0.01,
            "scan 1 isolation target_mz: expected ~471.2561, got {}",
            iw.0
        );

        // Scan 2: RT=125.3s, ms_level=2, isolation=(523.7832, 1.0, 1.0)
        let meta2 = idx.get_meta(2).expect("scan 2 should exist");
        assert_eq!(meta2.ms_level, 2);
        assert!(
            (meta2.rt_seconds - 125.3).abs() < 0.1,
            "scan 2 RT: expected ~125.3, got {}",
            meta2.rt_seconds
        );
    }

    #[test]
    fn extract_cv_value_basic() {
        let region = b"accession=\"MS:1000016\" value=\"120.5\" unitCvRef=\"UO\"";
        assert!((extract_cv_value(region).unwrap() - 120.5).abs() < 0.001);
    }

    #[test]
    fn extract_cv_value_no_value() {
        let region = b"accession=\"MS:1000016\" name=\"scan start time\"";
        assert!(extract_cv_value(region).is_none());
    }

    #[test]
    fn extract_meta_from_region_basic() {
        let xml = br#"<spectrum index="0" id="scan=1">
        <cvParam accession="MS:1000511" value="2"/>
        <scan>
            <cvParam accession="MS:1000016" value="120.5" unitAccession="UO:0000010"/>
        </scan>
        <isolationWindow>
            <cvParam accession="MS:1000827" value="500.0"/>
            <cvParam accession="MS:1000828" value="1.0"/>
            <cvParam accession="MS:1000829" value="1.0"/>
        </isolationWindow>
        <binaryDataArrayList>"#;
        let (rt, ms, iso) = extract_meta_from_region(xml);
        assert!((rt - 120.5).abs() < 0.001);
        assert_eq!(ms, 2);
        let (t, l, u) = iso.unwrap();
        assert!((t - 500.0).abs() < 0.001);
        assert!((l - 1.0).abs() < 0.001);
        assert!((u - 1.0).abs() < 0.001);
    }

    #[test]
    fn extract_meta_minutes_conversion() {
        let xml = br#"<cvParam accession="MS:1000511" value="1"/>
        <cvParam accession="MS:1000016" value="10.5" unitAccession="UO:0000031"/>"#;
        let (rt, ms, _) = extract_meta_from_region(xml);
        assert!(
            (rt - 630.0).abs() < 0.1,
            "10.5 min should be 630 seconds, got {rt}"
        );
        assert_eq!(ms, 1);
    }
}
