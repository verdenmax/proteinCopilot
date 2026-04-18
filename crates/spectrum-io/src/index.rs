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

/// Maps scan numbers to byte offsets within an mzML file.
///
/// Enables O(1) spectrum lookup: seek to offset → parse single `<spectrum>` node.
#[derive(Debug, Clone)]
pub struct ScanIndex {
    /// scan_number → byte offset of the `<spectrum>` opening tag.
    offsets: HashMap<u32, u64>,
    /// How this index was built.
    source: IndexSource,
}

impl ScanIndex {
    /// Creates a new ScanIndex from a pre-built map.
    pub fn new(offsets: HashMap<u32, u64>, source: IndexSource) -> Self {
        Self { offsets, source }
    }

    /// Returns the byte offset for a given scan number, or `None`.
    pub fn get_offset(&self, scan: u32) -> Option<u64> {
        self.offsets.get(&scan).copied()
    }

    /// Returns the number of indexed scans.
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// How this index was constructed.
    pub fn source(&self) -> IndexSource {
        self.source
    }

    /// Returns a reference to the underlying offsets map.
    pub fn offsets(&self) -> &HashMap<u32, u64> {
        &self.offsets
    }

    /// Returns all indexed scan numbers, sorted ascending.
    pub fn scan_numbers(&self) -> Vec<u32> {
        let mut scans: Vec<u32> = self.offsets.keys().copied().collect();
        scans.sort_unstable();
        scans
    }
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
    let (after_id, closing_quote) =
        if let Some(pos) = memchr::memmem::find(region, b" id=\"") {
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

    let mut reader = BufReader::with_capacity(256 * 1024, file);
    let needle = b"<spectrum ";
    let mut offsets = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut global_pos: u64 = 0;

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
            let abs_pos = global_pos + local_pos as u64;
            fallback_scan += 1;
            // Extract scan from tag bytes (limit to 512 bytes or end of buffer)
            let tag_end = (local_pos + 512).min(buf_len);
            let scan = extract_scan_from_tag_bytes(&buf[local_pos..tag_end], fallback_scan);
            if let Some(prev_offset) = offsets.insert(scan, abs_pos) {
                tracing::warn!(
                    "duplicate scan {} found while byte-scanning: offset {} replaced by {}",
                    scan,
                    prev_offset,
                    abs_pos
                );
            }
            search_start = local_pos + needle.len();
        }

        // Keep overlap for cross-boundary matches
        let overlap = needle.len() - 1;
        let consumed = if buf_len > overlap {
            buf_len - overlap
        } else {
            buf_len
        };
        global_pos += consumed as u64;
        reader.consume(consumed);
    }

    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
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
}
