# mzML Scan Index & Spectrum Cache — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add O(1) scan lookup to spectrum-io by parsing mzML native `<indexList>` and building scan offset indices, then cache IndexedReaders in the MCP server to eliminate repeated file parsing.

**Architecture:** New `ScanIndex` type maps scan→byte_offset. `IndexedMzMLReader` holds a `ScanIndex` + `PathBuf` and implements `SpectrumReader`. For indexed reads, it seeks the file to the offset, creates a fresh quick-xml `Reader`, and parses exactly one `<spectrum>` node. For files without native `<indexList>`, a fallback builds the index by scanning the file once. MCP server caches `IndexedMzMLReader` instances in an LRU. Existing `MzMLReader`/`MgfReader` are untouched (backward compatible).

**Tech Stack:** Rust, quick-xml 0.37, std::io::Seek, lru crate (new dep)

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `crates/spectrum-io/src/index.rs` | Create | `ScanIndex`, `IndexSource` types + `build_index_from_native_mzml()` + `build_index_by_scanning()` |
| `crates/spectrum-io/src/indexed_mzml.rs` | Create | `IndexedMzMLReader` struct + `SpectrumReader` impl + `read_spectrum_at_offset()` |
| `crates/spectrum-io/src/indexed_mgf.rs` | Create | `IndexedMgfReader` struct + `SpectrumReader` impl |
| `crates/spectrum-io/src/lib.rs` | Modify | Add `pub mod index; pub mod indexed_mzml; pub mod indexed_mgf;` + `create_indexed_reader()` |
| `crates/spectrum-io/src/error.rs` | Modify | Add `IndexParseError` variant |
| `crates/spectrum-io/Cargo.toml` | Modify | — (no new deps needed; `std::io::Seek` is stdlib) |
| `crates/mcp-server/Cargo.toml` | Modify | Add `lru` dependency |
| `crates/mcp-server/src/tools.rs` | Modify | Add `reader_cache` field + `get_or_create_reader()` + migrate tools |
| `crates/spectrum-io/tests/fixtures/small_indexed.mzml` | Create | Test fixture with `<indexedmzML>` wrapper and `<indexList>` |

---

## Task 1: ScanIndex data types and native indexList parser

**Files:**
- Create: `crates/spectrum-io/src/index.rs`
- Modify: `crates/spectrum-io/src/lib.rs`
- Modify: `crates/spectrum-io/src/error.rs`

- [ ] **Step 1: Add IndexParseError variant to error.rs**

In `crates/spectrum-io/src/error.rs`, add a new variant to `SpectrumIoError` enum (after `ScanNotFound`):

```rust
    /// Index parsing error (mzML indexList).
    #[error("index parse error in {path}: {detail}")]
    IndexParseError {
        /// The file path.
        path: PathBuf,
        /// What went wrong.
        detail: String,
    },
```

Also add its conversion in the `From<SpectrumIoError> for CoreError` impl (after the `ScanNotFound` arm):

```rust
            SpectrumIoError::IndexParseError { path, detail } => {
                protein_copilot_core::error::CoreError::SpectrumParseError {
                    format: "mzML".to_string(),
                    detail: format!("{}: {detail}", path.display()),
                    suggestion: "The mzML index may be corrupted. Try re-converting the file."
                        .to_string(),
                }
            }
```

- [ ] **Step 2: Create index.rs with ScanIndex and IndexSource**

Create `crates/spectrum-io/src/index.rs`:

```rust
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

    /// Returns all indexed scan numbers, sorted ascending.
    pub fn scan_numbers(&self) -> Vec<u32> {
        let mut scans: Vec<u32> = self.offsets.keys().copied().collect();
        scans.sort_unstable();
        scans
    }
}
```

- [ ] **Step 3: Add `build_index_from_native_mzml()` to index.rs**

This function reads the end of the file to find `<indexListOffset>`, seeks there, and parses all `<offset>` entries. Append to `index.rs`:

```rust
/// Size of the tail chunk to read when searching for `<indexListOffset>`.
/// The indexListOffset element is always near the very end of the file.
const TAIL_READ_SIZE: usize = 4096;

/// Attempts to build a ScanIndex from the native `<indexList>` in an `<indexedmzML>` file.
///
/// Returns `Ok(Some(index))` if the file has a valid `<indexList>`,
/// `Ok(None)` if it's a plain `<mzML>` without an index,
/// `Err(...)` on I/O or parse errors.
pub fn build_index_from_native_mzml(path: &Path) -> Result<Option<ScanIndex>, SpectrumIoError> {
    use std::fs::File;

    let mut file = File::open(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            SpectrumIoError::FileNotFound { path: path.to_path_buf() }
        } else {
            SpectrumIoError::IoError { path: path.to_path_buf(), source: e }
        }
    })?;

    let file_len = file.metadata()
        .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?
        .len();

    if file_len == 0 {
        return Ok(None);
    }

    // Read the tail of the file to find <indexListOffset>
    let tail_start = file_len.saturating_sub(TAIL_READ_SIZE as u64);
    file.seek(SeekFrom::Start(tail_start))
        .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;

    let mut tail = String::new();
    file.read_to_string(&mut tail)
        .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;

    // Look for <indexListOffset>NNN</indexListOffset>
    let offset_str = match extract_between(&tail, "<indexListOffset>", "</indexListOffset>") {
        Some(s) => s,
        None => return Ok(None), // Not an indexedmzML file
    };

    let index_list_offset: u64 = offset_str.trim().parse().map_err(|_| {
        SpectrumIoError::IndexParseError {
            path: path.to_path_buf(),
            detail: format!("invalid indexListOffset value: '{offset_str}'"),
        }
    })?;

    // Seek to the indexList and parse it
    file.seek(SeekFrom::Start(index_list_offset))
        .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;

    let mut index_xml = String::new();
    file.read_to_string(&mut index_xml)
        .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;

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
///
/// Expects XML like:
/// ```xml
/// <indexList count="N">
///   <index name="spectrum">
///     <offset idRef="scan=1">4523</offset>
///     <offset idRef="controllerType=0 controllerNumber=1 scan=2">18904</offset>
///   </index>
/// </indexList>
/// ```
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
                        // Check if this is the spectrum index (not chromatogram)
                        if let Some(name) = e.attributes()
                            .filter_map(|a| a.ok())
                            .find(|a| a.key.as_ref() == b"name")
                            .and_then(|a| String::from_utf8(a.value.to_vec()).ok())
                        {
                            in_spectrum_index = name == "spectrum";
                        }
                    }
                    b"offset" if in_spectrum_index => {
                        in_offset = true;
                        current_id_ref = e.attributes()
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
                    let scan = parse_scan_from_id_ref(&current_id_ref)
                        .unwrap_or(fallback_scan);
                    offsets.insert(scan, byte_offset);
                }
            }
            Ok(Event::End(ref e)) => {
                match e.local_name().as_ref() {
                    b"offset" => in_offset = false,
                    b"index" => in_spectrum_index = false,
                    _ => {}
                }
            }
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
    id_ref
        .split("scan=")
        .nth(1)
        .and_then(|s| {
            let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
            digits.parse().ok()
        })
}
```

- [ ] **Step 4: Add `build_index_by_scanning()` fallback to index.rs**

This scans the entire file tracking `<spectrum>` byte positions. Append to `index.rs`:

```rust
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
            SpectrumIoError::FileNotFound { path: path.to_path_buf() }
        } else {
            SpectrumIoError::IoError { path: path.to_path_buf(), source: e }
        }
    })?;

    let mut reader = BufReader::new(file);
    let mut offsets = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut byte_pos: u64 = 0;

    // Read line-by-line (mzML is line-oriented in practice)
    let mut line = String::new();
    loop {
        let line_start = byte_pos;
        line.clear();
        let bytes_read = reader.read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;
        if bytes_read == 0 {
            break; // EOF
        }
        byte_pos += bytes_read as u64;

        let trimmed = line.trim();

        // Look for <spectrum lines — either <spectrum ... > or <spectrum .../>
        if trimmed.starts_with("<spectrum ") || trimmed.starts_with("<spectrum>") {
            fallback_scan += 1;

            // Try to extract scan from id="..." attribute
            let scan = extract_id_attr(trimmed)
                .and_then(|id| parse_scan_from_id_ref(&id))
                .unwrap_or(fallback_scan);

            offsets.insert(scan, line_start);
        }
    }

    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
}

/// Extracts the value of the `id` attribute from an XML tag string.
fn extract_id_attr(tag_text: &str) -> Option<String> {
    // Look for id="..." or id='...'
    let after_id = tag_text.split("id=\"").nth(1)
        .or_else(|| tag_text.split("id='").nth(1))?;
    let end = after_id.find('"').or_else(|| after_id.find('\''))?;
    Some(after_id[..end].to_string())
}
```

- [ ] **Step 5: Register module in lib.rs**

In `crates/spectrum-io/src/lib.rs`, add after `mod util;`:

```rust
pub mod index;
```

- [ ] **Step 6: Write tests for ScanIndex**

Append to `crates/spectrum-io/src/index.rs`:

```rust
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
        // small.mzml is a plain <mzML> file without <indexList>
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small.mzml");
        let result = build_index_from_native_mzml(&path).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn build_index_by_scanning_finds_spectra() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small.mzml");
        let idx = build_index_by_scanning(&path).unwrap();
        assert_eq!(idx.len(), 10); // small.mzml has 10 spectra
        assert_eq!(idx.source(), IndexSource::BuiltFromScan);
        // Scan 1 should exist
        assert!(idx.get_offset(1).is_some());
        // Scan 10 should exist
        assert!(idx.get_offset(10).is_some());
    }
}
```

- [ ] **Step 7: Run tests**

Run: `cargo test -p protein-copilot-spectrum-io -- index 2>&1 | tail -20`
Expected: all tests pass

- [ ] **Step 8: Commit**

```bash
git add crates/spectrum-io/src/index.rs crates/spectrum-io/src/error.rs crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add ScanIndex with native mzML indexList parser

Introduces ScanIndex type mapping scan number → byte offset.
Two construction paths:
- build_index_from_native_mzml(): reads <indexList> from <indexedmzML> EOF
- build_index_by_scanning(): line-by-line fallback for plain <mzML>

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 2: Create indexed mzML test fixture

**Files:**
- Create: `crates/spectrum-io/tests/fixtures/small_indexed.mzml`

- [ ] **Step 1: Create the indexed mzML fixture**

This is a copy of `small.mzml` but wrapped in `<indexedmzML>` with a real `<indexList>` at the end. To build it correctly, we need actual byte offsets. The simplest approach: write a small helper script that reads small.mzml, wraps it, and computes real offsets.

Create `crates/spectrum-io/tests/fixtures/small_indexed.mzml` by transforming the existing fixture. The approach:

1. Read `small.mzml` content
2. Compute the byte offset of each `<spectrum` tag after wrapping in `<indexedmzML>`
3. Write the wrapped file with the correct `<indexList>` and `<indexListOffset>`

Write a one-off Rust integration test that generates this fixture and verifies it. Place this logic in a `#[test]` that only runs when the fixture doesn't exist yet:

```rust
// In crates/spectrum-io/src/index.rs tests section, add:
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
    // Verify all 10 scans are indexed
    for scan in 1..=10 {
        assert!(idx.get_offset(scan).is_some(), "missing scan {scan}");
    }
}

/// Generates small_indexed.mzml from small.mzml by wrapping in <indexedmzML>.
fn generate_indexed_fixture(output_path: &std::path::Path) {
    use std::io::Write;

    let source_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/small.mzml");
    let source = std::fs::read_to_string(&source_path).unwrap();

    // Remove XML declaration if present (we'll add our own)
    let mzml_content = source
        .strip_prefix("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
        .unwrap_or(&source);

    let header = "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\">\n";

    // Build body = header + mzml content
    let body = format!("{}{}", header, mzml_content);

    // Find all <spectrum byte offsets in body
    let mut offsets: Vec<(u32, usize)> = Vec::new();
    let mut search_start = 0;
    let mut fallback_scan = 0u32;
    while let Some(pos) = body[search_start..].find("<spectrum ") {
        let abs_pos = search_start + pos;
        fallback_scan += 1;
        // Extract scan from id attribute
        let tag_end = body[abs_pos..].find('>').unwrap_or(200) + abs_pos;
        let tag_text = &body[abs_pos..tag_end];
        let scan = super::extract_id_attr(tag_text)
            .and_then(|id| super::parse_scan_from_id_ref(&id))
            .unwrap_or(fallback_scan);
        offsets.push((scan, abs_pos));
        search_start = abs_pos + 1;
    }

    // Build indexList
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
        index_list,
        index_list_offset,
    );

    let mut out = std::fs::File::create(output_path).unwrap();
    write!(out, "{}{}", body, footer).unwrap();
}
```

- [ ] **Step 2: Run the fixture generation test**

Run: `cargo test -p protein-copilot-spectrum-io -- build_index_from_native_indexed_mzml 2>&1 | tail -10`
Expected: PASS, fixture file created

- [ ] **Step 3: Verify fixture is valid**

Run: `head -5 crates/spectrum-io/tests/fixtures/small_indexed.mzml && echo "..." && tail -20 crates/spectrum-io/tests/fixtures/small_indexed.mzml`
Expected: starts with `<indexedmzML>`, ends with `<indexList>` + `<indexListOffset>`

- [ ] **Step 4: Commit**

```bash
git add crates/spectrum-io/tests/fixtures/small_indexed.mzml crates/spectrum-io/src/index.rs
git commit -m "test(spectrum-io): add indexed mzML test fixture

Auto-generated from small.mzml with correct byte offsets.
Wrapped in <indexedmzML> with <indexList> for native index testing.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 3: IndexedMzMLReader with seek-based spectrum reading

**Files:**
- Create: `crates/spectrum-io/src/indexed_mzml.rs`
- Modify: `crates/spectrum-io/src/lib.rs`

- [ ] **Step 1: Create indexed_mzml.rs**

Create `crates/spectrum-io/src/indexed_mzml.rs`:

```rust
//! Indexed mzML reader — O(1) spectrum lookup via scan index.
//!
//! Uses [`ScanIndex`] to seek directly to a spectrum's byte offset,
//! then parses just that one `<spectrum>` node. For operations that
//! need all spectra (read_all, read_summary, for_each_spectrum), delegates
//! to the standard [`MzMLReader`].

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use quick_xml::events::Event;
use quick_xml::Reader;

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;
use crate::index::{build_index_by_scanning, build_index_from_native_mzml, IndexSource, ScanIndex};
use crate::mzml::MzMLReader;
use crate::reader::SpectrumReader;

/// Indexed mzML reader with O(1) scan lookup.
///
/// On construction, builds or loads a [`ScanIndex`] mapping scan numbers
/// to byte offsets. Subsequent `read_spectrum()` calls seek to the offset
/// and parse a single `<spectrum>` node instead of scanning from the start.
pub struct IndexedMzMLReader {
    index: ScanIndex,
    path: PathBuf,
}

impl IndexedMzMLReader {
    /// Opens an mzML file and builds a scan index.
    ///
    /// First tries to parse the native `<indexList>` (fast, reads only EOF).
    /// Falls back to scanning the entire file to build the index.
    pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
        let index = if let Some(native) = build_index_from_native_mzml(path)? {
            tracing::debug!(
                path = %path.display(),
                scans = native.len(),
                "loaded native mzML index"
            );
            native
        } else {
            let built = build_index_by_scanning(path)?;
            tracing::debug!(
                path = %path.display(),
                scans = built.len(),
                "built scan index by file scanning"
            );
            built
        };

        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }

    /// Returns a reference to the underlying scan index.
    pub fn index(&self) -> &ScanIndex {
        &self.index
    }

    /// Reads a single spectrum by seeking to its byte offset.
    ///
    /// Opens the file, seeks to the offset, creates a new XML reader,
    /// and parses exactly one `<spectrum>` node.
    fn read_spectrum_at_offset(&self, scan: u32, offset: u64) -> Result<Spectrum, SpectrumIoError> {
        let file = File::open(&self.path).map_err(|e| SpectrumIoError::IoError {
            path: self.path.clone(),
            source: e,
        })?;
        let mut buf_reader = BufReader::new(file);
        buf_reader.seek(SeekFrom::Start(offset)).map_err(|e| SpectrumIoError::IoError {
            path: self.path.clone(),
            source: e,
        })?;

        let mut xml_reader = Reader::from_reader(buf_reader);
        xml_reader.config_mut().trim_text(true);

        // Parse exactly one <spectrum> node using the existing streaming parser.
        // We reuse parse_single_spectrum which reads until </spectrum>.
        parse_single_spectrum(&mut xml_reader, &self.path, scan)
    }
}

impl SpectrumReader for IndexedMzMLReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        // Delegate to standard reader — index doesn't help here
        MzMLReader.read_all(path)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        MzMLReader.read_summary(path)
    }

    fn read_spectrum(&self, _path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        match self.index.get_offset(scan) {
            Some(offset) => self.read_spectrum_at_offset(scan, offset),
            None => Err(SpectrumIoError::ScanNotFound {
                path: self.path.clone(),
                scan,
            }),
        }
    }

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        // Delegate to standard reader — streaming doesn't benefit from index
        MzMLReader.for_each_spectrum(path, handler)
    }
}
```

- [ ] **Step 2: Implement `parse_single_spectrum()`**

This function parses exactly one `<spectrum>` node from the current reader position. It reuses the same XML parsing logic as `parse_mzml_streaming` but stops after the first `</spectrum>`. Append to `indexed_mzml.rs`:

```rust
/// Parses a single `<spectrum>` node from the current reader position.
///
/// Expects the reader to be positioned at or just before a `<spectrum>` start tag.
/// Reads until `</spectrum>` and returns the parsed Spectrum.
fn parse_single_spectrum<R: std::io::BufRead>(
    xml_reader: &mut Reader<R>,
    path: &Path,
    expected_scan: u32,
) -> Result<Spectrum, SpectrumIoError> {
    use crate::mzml::{decode_binary_array_pub, get_attr_pub, parse_scan_from_id_pub};
    use protein_copilot_core::spectrum::{IsolationWindow, MsLevel, PrecursorInfo};

    let mut buf = Vec::new();
    let mut in_spectrum = false;
    let mut in_precursor = false;
    let mut in_isolation_window = false;
    let mut in_selected_ion = false;
    let mut in_binary_data_array = false;
    let mut in_binary = false;
    let mut in_scan = false;

    // SpectrumBuilder fields (inline to avoid coupling)
    let mut scan_number: Option<u32> = None;
    let mut ms_level: Option<u8> = None;
    let mut rt_sec: Option<f64> = None;
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

    struct BinaryMeta {
        is_mz: bool,
        is_intensity: bool,
        is_64bit: bool,
        is_zlib: bool,
    }
    let mut array_meta = BinaryMeta { is_mz: false, is_intensity: false, is_64bit: false, is_zlib: false };
    let mut binary_text = String::new();

    loop {
        match xml_reader.read_event_into(&mut buf) {
            Ok(Event::Eof) => {
                return Err(SpectrumIoError::ScanNotFound {
                    path: path.to_path_buf(),
                    scan: expected_scan,
                });
            }
            Ok(Event::Start(ref e)) => {
                let tag = e.local_name();
                match tag.as_ref() {
                    b"spectrum" if !in_spectrum => {
                        in_spectrum = true;
                        scan_number = get_attr_pub(e, b"id")
                            .and_then(|id| parse_scan_from_id_pub(&id));
                    }
                    b"scan" if in_spectrum => in_scan = true,
                    b"precursor" if in_spectrum => {
                        in_precursor = true;
                        if let Some(sref) = get_attr_pub(e, b"spectrumRef") {
                            cur_precursor_source_scan = crate::mzml::parse_scan_from_spectrum_ref_pub(&sref);
                        }
                    }
                    b"isolationWindow" if in_precursor => in_isolation_window = true,
                    b"selectedIon" if in_precursor => in_selected_ion = true,
                    b"binaryDataArray" if in_spectrum => {
                        in_binary_data_array = true;
                        array_meta = BinaryMeta { is_mz: false, is_intensity: false, is_64bit: false, is_zlib: false };
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
                if e.local_name().as_ref() == b"cvParam" {
                    let acc = get_attr_pub(e, b"accession").unwrap_or_default();
                    let value = get_attr_pub(e, b"value").unwrap_or_default();
                    match acc.as_str() {
                        "MS:1000511" if in_spectrum && !in_precursor => {
                            ms_level = value.parse().ok();
                        }
                        "MS:1000016" if in_scan => {
                            if let Ok(rv) = value.parse::<f64>() {
                                let unit = get_attr_pub(e, b"unitAccession").unwrap_or_default();
                                rt_sec = Some(if unit == "UO:0000031" { rv * 60.0 } else { rv });
                            }
                        }
                        "MS:1000827" if in_isolation_window => { cur_isolation_target_mz = value.parse().ok(); }
                        "MS:1000828" if in_isolation_window => { cur_isolation_lower = value.parse().ok(); }
                        "MS:1000829" if in_isolation_window => { cur_isolation_upper = value.parse().ok(); }
                        "MS:1000744" if in_selected_ion => { cur_precursor_mz = value.parse().ok(); }
                        "MS:1000041" if in_selected_ion => { cur_precursor_charge = value.parse().ok(); }
                        "MS:1000042" if in_selected_ion => { cur_precursor_intensity = value.parse().ok(); }
                        "MS:1000514" if in_binary_data_array => { array_meta.is_mz = true; }
                        "MS:1000515" if in_binary_data_array => { array_meta.is_intensity = true; }
                        "MS:1000523" if in_binary_data_array => { array_meta.is_64bit = true; }
                        "MS:1000521" if in_binary_data_array => { array_meta.is_64bit = false; }
                        "MS:1000574" if in_binary_data_array => { array_meta.is_zlib = true; }
                        "MS:1000576" if in_binary_data_array => { array_meta.is_zlib = false; }
                        _ => {}
                    }
                }
            }
            Ok(Event::End(ref e)) => {
                match e.local_name().as_ref() {
                    b"spectrum" => {
                        // Done! Build and return.
                        let scan = scan_number.unwrap_or(expected_scan);
                        let level = match ms_level.unwrap_or(2) {
                            1 => MsLevel::MS1,
                            2 => MsLevel::MS2,
                            n => MsLevel::Other(n),
                        };
                        crate::util::sort_peaks_by_mz(&mut mz_array, &mut intensity_array);
                        return Spectrum::new(
                            scan, level, rt_sec.unwrap_or(0.0),
                            precursors, mz_array, intensity_array,
                        ).map_err(|e| SpectrumIoError::ValidationError {
                            scan, detail: e.to_string(),
                        });
                    }
                    b"scan" => in_scan = false,
                    b"precursor" => {
                        in_precursor = false;
                        // Flush precursor
                        if let Some(mz) = cur_precursor_mz.take() {
                            let iw = match (cur_isolation_target_mz.take(), cur_isolation_lower.take(), cur_isolation_upper.take()) {
                                (Some(t), Some(l), Some(u)) => Some(IsolationWindow { target_mz: t, lower_offset: l, upper_offset: u }),
                                _ => None,
                            };
                            precursors.push(PrecursorInfo {
                                mz,
                                charge: cur_precursor_charge.take(),
                                intensity: cur_precursor_intensity.take(),
                                isolation_window: iw,
                                source_scan: cur_precursor_source_scan.take(),
                            });
                        }
                    }
                    b"isolationWindow" => in_isolation_window = false,
                    b"selectedIon" => in_selected_ion = false,
                    b"binary" => in_binary = false,
                    b"binaryDataArray" => {
                        if !binary_text.is_empty() {
                            let meta = crate::mzml::BinaryArrayMetaPub {
                                is_mz: array_meta.is_mz,
                                is_intensity: array_meta.is_intensity,
                                is_64bit: array_meta.is_64bit,
                                is_zlib: array_meta.is_zlib,
                            };
                            let decoded = decode_binary_array_pub(&binary_text, &meta, path)?;
                            if array_meta.is_mz { mz_array = decoded; }
                            else if array_meta.is_intensity { intensity_array = decoded; }
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
                    detail: format!("XML error at position {}: {e}", xml_reader.error_position()),
                });
            }
            _ => {}
        }
        buf.clear();
    }
}
```

**NOTE:** This step requires exposing some helper functions from `mzml.rs` as `pub(crate)`. The implementer must add `pub(crate)` visibility to `get_attr`, `parse_scan_from_id`, `parse_scan_from_spectrum_ref`, `decode_binary_array`, and `BinaryArrayMeta` in `mzml.rs`. Specifically:

- Rename `get_attr` → keep private, add `pub(crate) fn get_attr_pub(...)` that delegates
- Or simpler: change `fn get_attr` to `pub(crate) fn get_attr` and re-export with `_pub` suffix aliases
- Best approach: make the existing internal functions `pub(crate)` directly:
  - `get_attr` → `pub(crate) fn get_attr`
  - `parse_scan_from_id` → `pub(crate) fn parse_scan_from_id`
  - `parse_scan_from_spectrum_ref` → `pub(crate) fn parse_scan_from_spectrum_ref`
  - `decode_binary_array` → `pub(crate) fn decode_binary_array`
  - `BinaryArrayMeta` → `pub(crate) struct BinaryArrayMeta`

Then in `indexed_mzml.rs`, import them directly:
```rust
use crate::mzml::{decode_binary_array, get_attr, parse_scan_from_id, parse_scan_from_spectrum_ref, BinaryArrayMeta};
```

- [ ] **Step 3: Register module in lib.rs**

Add to `crates/spectrum-io/src/lib.rs` after `pub mod index;`:

```rust
pub mod indexed_mzml;
pub use indexed_mzml::IndexedMzMLReader;
```

- [ ] **Step 4: Write tests**

Append to `indexed_mzml.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small.mzml")
    }

    fn indexed_fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small_indexed.mzml")
    }

    #[test]
    fn open_plain_mzml_builds_scan_index() {
        let reader = IndexedMzMLReader::open(&fixture_path()).unwrap();
        assert_eq!(reader.index().len(), 10);
        assert_eq!(reader.index().source(), IndexSource::BuiltFromScan);
    }

    #[test]
    fn open_indexed_mzml_uses_native_index() {
        let path = indexed_fixture_path();
        if !path.exists() {
            return; // Skip if fixture not generated yet
        }
        let reader = IndexedMzMLReader::open(&path).unwrap();
        assert_eq!(reader.index().len(), 10);
        assert_eq!(reader.index().source(), IndexSource::NativeIndex);
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
        // Verify same content as standard reader
        let standard = MzMLReader.read_spectrum(&fixture_path(), 7).unwrap();
        assert_eq!(spec.mz_array.len(), standard.mz_array.len());
        assert!((spec.retention_time_sec - standard.retention_time_sec).abs() < 0.01);
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
            assert_eq!(idx_spec.intensity_array.len(), std_spec.intensity_array.len());
            assert!((idx_spec.retention_time_sec - std_spec.retention_time_sec).abs() < 0.001);
        }
    }
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p protein-copilot-spectrum-io -- indexed_mzml 2>&1 | tail -20`
Expected: all tests pass

- [ ] **Step 6: Commit**

```bash
git add crates/spectrum-io/src/indexed_mzml.rs crates/spectrum-io/src/mzml.rs crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add IndexedMzMLReader with O(1) scan lookup

Seeks to byte offset from ScanIndex, parses single <spectrum> node.
Implements full SpectrumReader trait (delegates bulk ops to MzMLReader).
All 10 scans verified identical to standard reader output.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 4: IndexedMgfReader

**Files:**
- Create: `crates/spectrum-io/src/indexed_mgf.rs`
- Modify: `crates/spectrum-io/src/lib.rs`

- [ ] **Step 1: Create indexed_mgf.rs**

MGF indexing is simpler — just record the byte offset of each `BEGIN IONS` line.

Create `crates/spectrum-io/src/indexed_mgf.rs`:

```rust
//! Indexed MGF reader — O(1) spectrum lookup via line offset index.
//!
//! Scans the file once to record byte offsets of `BEGIN IONS` lines,
//! then seeks directly for subsequent reads.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use protein_copilot_core::spectrum::{MsLevel, PrecursorInfo, Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;
use crate::index::{IndexSource, ScanIndex};
use crate::mgf::MgfReader;
use crate::reader::SpectrumReader;

/// Indexed MGF reader with O(1) scan lookup.
pub struct IndexedMgfReader {
    index: ScanIndex,
    path: PathBuf,
}

impl IndexedMgfReader {
    /// Opens an MGF file and builds a scan index by recording
    /// `BEGIN IONS` byte offsets.
    pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
        let index = build_mgf_index(path)?;
        tracing::debug!(
            path = %path.display(),
            scans = index.len(),
            "built MGF scan index"
        );
        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }

    /// Returns a reference to the underlying scan index.
    pub fn index(&self) -> &ScanIndex {
        &self.index
    }

    /// Reads a single spectrum by seeking to its byte offset.
    fn read_spectrum_at_offset(&self, scan: u32, offset: u64) -> Result<Spectrum, SpectrumIoError> {
        let file = File::open(&self.path).map_err(|e| SpectrumIoError::IoError {
            path: self.path.clone(),
            source: e,
        })?;
        let mut reader = BufReader::new(file);
        reader.seek(SeekFrom::Start(offset)).map_err(|e| SpectrumIoError::IoError {
            path: self.path.clone(),
            source: e,
        })?;

        // Parse one MGF block: BEGIN IONS ... END IONS
        parse_single_mgf_block(&mut reader, &self.path, scan)
    }
}

impl SpectrumReader for IndexedMgfReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        MgfReader.read_all(path)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        MgfReader.read_summary(path)
    }

    fn read_spectrum(&self, _path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        match self.index.get_offset(scan) {
            Some(offset) => self.read_spectrum_at_offset(scan, offset),
            None => Err(SpectrumIoError::ScanNotFound {
                path: self.path.clone(),
                scan,
            }),
        }
    }

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        MgfReader.for_each_spectrum(path, handler)
    }
}

/// Builds a scan index for an MGF file by recording `BEGIN IONS` byte offsets.
fn build_mgf_index(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    use std::collections::HashMap;

    let file = File::open(path).map_err(|e| {
        if e.kind() == std::io::ErrorKind::NotFound {
            SpectrumIoError::FileNotFound { path: path.to_path_buf() }
        } else {
            SpectrumIoError::IoError { path: path.to_path_buf(), source: e }
        }
    })?;
    let mut reader = BufReader::new(file);
    let mut offsets = HashMap::new();
    let mut fallback_scan: u32 = 0;
    let mut byte_pos: u64 = 0;
    let mut line = String::new();

    // First pass: find BEGIN IONS positions
    let mut begin_positions: Vec<u64> = Vec::new();
    loop {
        let line_start = byte_pos;
        line.clear();
        let n = reader.read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;
        if n == 0 { break; }
        byte_pos += n as u64;
        if line.trim() == "BEGIN IONS" {
            begin_positions.push(line_start);
        }
    }

    // Second pass: for each BEGIN IONS, read ahead to find SCANS= header
    for &pos in &begin_positions {
        fallback_scan += 1;
        reader.seek(SeekFrom::Start(pos))
            .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;

        let mut scan = fallback_scan;
        line.clear();
        // Skip the BEGIN IONS line itself
        reader.read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;

        // Read header lines until we hit a peak line or END IONS
        for _ in 0..20 {
            line.clear();
            let n = reader.read_line(&mut line)
                .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;
            if n == 0 { break; }
            let trimmed = line.trim();
            if trimmed == "END IONS" || trimmed.is_empty() { break; }
            if let Some((key, val)) = trimmed.split_once('=') {
                if key.trim().eq_ignore_ascii_case("SCANS") {
                    if let Ok(s) = val.trim().parse::<u32>() {
                        scan = s;
                    }
                    break;
                }
            }
            // If line looks like a peak (starts with a digit), stop looking for headers
            if trimmed.starts_with(|c: char| c.is_ascii_digit()) { break; }
        }
        offsets.insert(scan, pos);
    }

    Ok(ScanIndex::new(offsets, IndexSource::BuiltFromScan))
}

/// Parses a single MGF block starting from `BEGIN IONS`.
fn parse_single_mgf_block<R: BufRead>(
    reader: &mut R,
    path: &Path,
    expected_scan: u32,
) -> Result<Spectrum, SpectrumIoError> {
    let mut line = String::new();
    let mut pepmass_mz: Option<f64> = None;
    let mut pepmass_int: Option<f64> = None;
    let mut charge: Option<i32> = None;
    let mut rt_sec: Option<f64> = None;
    let mut scan: Option<u32> = None;
    let mut mz_values = Vec::new();
    let mut intensity_values = Vec::new();
    let mut found_begin = false;

    loop {
        line.clear();
        let n = reader.read_line(&mut line)
            .map_err(|e| SpectrumIoError::IoError { path: path.to_path_buf(), source: e })?;
        if n == 0 { break; }
        let trimmed = line.trim();

        if trimmed == "BEGIN IONS" {
            found_begin = true;
            continue;
        }

        if !found_begin { continue; }

        if trimmed == "END IONS" {
            break;
        }

        if let Some((key, val)) = trimmed.split_once('=') {
            match key.trim().to_uppercase().as_str() {
                "PEPMASS" => {
                    let parts: Vec<&str> = val.split_whitespace().collect();
                    pepmass_mz = parts.first().and_then(|v| v.parse().ok());
                    pepmass_int = parts.get(1).and_then(|v| v.parse().ok());
                }
                "CHARGE" => {
                    let s = val.trim();
                    charge = if let Some(n) = s.strip_suffix('+') {
                        n.trim().parse().ok()
                    } else if let Some(n) = s.strip_suffix('-') {
                        n.trim().parse::<i32>().ok().map(|v| -v)
                    } else {
                        s.parse().ok()
                    };
                }
                "RTINSECONDS" => { rt_sec = val.trim().parse().ok(); }
                "SCANS" => { scan = val.trim().parse().ok(); }
                _ => {}
            }
        } else {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(mz), Ok(int)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                    mz_values.push(mz);
                    intensity_values.push(int);
                }
            }
        }
    }

    let mz = pepmass_mz.ok_or_else(|| SpectrumIoError::ParseError {
        path: path.to_path_buf(),
        line: 0,
        detail: format!("scan {expected_scan} missing PEPMASS"),
    })?;

    crate::util::sort_peaks_by_mz(&mut mz_values, &mut intensity_values);

    Spectrum::new(
        scan.unwrap_or(expected_scan),
        MsLevel::MS2,
        rt_sec.unwrap_or(0.0),
        vec![PrecursorInfo {
            mz,
            charge,
            intensity: pepmass_int,
            isolation_window: None,
            source_scan: None,
        }],
        mz_values,
        intensity_values,
    ).map_err(|e| SpectrumIoError::ValidationError {
        scan: expected_scan,
        detail: e.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small.mgf")
    }

    #[test]
    fn open_builds_index() {
        let reader = IndexedMgfReader::open(&fixture_path()).unwrap();
        assert_eq!(reader.index().len(), 10);
    }

    #[test]
    fn read_spectrum_by_index() {
        let reader = IndexedMgfReader::open(&fixture_path()).unwrap();
        let spec = reader.read_spectrum(&fixture_path(), 1).unwrap();
        assert_eq!(spec.scan_number, 1);
    }

    #[test]
    fn indexed_vs_standard_all_scans_match() {
        let indexed = IndexedMgfReader::open(&fixture_path()).unwrap();
        for scan in 1..=10 {
            let idx_spec = indexed.read_spectrum(&fixture_path(), scan).unwrap();
            let std_spec = MgfReader.read_spectrum(&fixture_path(), scan).unwrap();
            assert_eq!(idx_spec.scan_number, std_spec.scan_number);
            assert_eq!(idx_spec.mz_array.len(), std_spec.mz_array.len());
            assert!((idx_spec.retention_time_sec - std_spec.retention_time_sec).abs() < 0.001);
        }
    }

    #[test]
    fn read_spectrum_not_found() {
        let reader = IndexedMgfReader::open(&fixture_path()).unwrap();
        let err = reader.read_spectrum(&fixture_path(), 999).unwrap_err();
        assert!(err.to_string().contains("999"));
    }
}
```

- [ ] **Step 2: Register module in lib.rs**

Add to `crates/spectrum-io/src/lib.rs`:

```rust
pub mod indexed_mgf;
pub use indexed_mgf::IndexedMgfReader;
```

- [ ] **Step 3: Run tests**

Run: `cargo test -p protein-copilot-spectrum-io -- indexed_mgf 2>&1 | tail -15`
Expected: all tests pass

- [ ] **Step 4: Commit**

```bash
git add crates/spectrum-io/src/indexed_mgf.rs crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add IndexedMgfReader with O(1) scan lookup

Scans MGF file once to record BEGIN IONS byte offsets.
Subsequent read_spectrum() calls seek directly to the block.
All 10 scans verified identical to standard MgfReader output.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 5: `create_indexed_reader()` factory + full test suite

**Files:**
- Modify: `crates/spectrum-io/src/lib.rs`

- [ ] **Step 1: Add `create_indexed_reader()` factory**

In `crates/spectrum-io/src/lib.rs`, add after `create_reader()`:

```rust
/// Creates an [`IndexedMzMLReader`] or [`IndexedMgfReader`] for the given file.
///
/// Indexed readers build a scan→offset map on first open, enabling
/// O(1) `read_spectrum()` calls. For operations that need all spectra
/// (read_all, for_each_spectrum), they delegate to the standard reader.
///
/// Prefer this over [`create_reader`] when you'll call `read_spectrum()`
/// multiple times on the same file.
pub fn create_indexed_reader(path: &Path) -> Result<Box<dyn SpectrumReader>, SpectrumIoError> {
    let info = detect_format(path)?;
    match info.format {
        SpectrumFormat::MzML => {
            let reader = IndexedMzMLReader::open(path)?;
            Ok(Box::new(reader))
        }
        SpectrumFormat::Mgf => {
            let reader = IndexedMgfReader::open(path)?;
            Ok(Box::new(reader))
        }
    }
}
```

- [ ] **Step 2: Add integration-style tests**

Add to the `tests` module in `lib.rs`:

```rust
    #[test]
    fn create_indexed_reader_mzml() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small.mzml");
        let reader = create_indexed_reader(&path).unwrap();
        let spec = reader.read_spectrum(&path, 5).unwrap();
        assert_eq!(spec.scan_number, 5);
    }

    #[test]
    fn create_indexed_reader_mgf() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/small.mgf");
        let reader = create_indexed_reader(&path).unwrap();
        let spec = reader.read_spectrum(&path, 3).unwrap();
        assert_eq!(spec.scan_number, 3);
    }
```

- [ ] **Step 3: Run all spectrum-io tests**

Run: `cargo test -p protein-copilot-spectrum-io 2>&1 | tail -15`
Expected: all tests pass (existing + new)

- [ ] **Step 4: Run clippy**

Run: `cargo clippy -p protein-copilot-spectrum-io -- -D warnings 2>&1 | tail -10`
Expected: 0 warnings

- [ ] **Step 5: Commit**

```bash
git add crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add create_indexed_reader() factory function

Public API for creating indexed readers that support O(1) scan lookup.
Auto-detects format (mzML/MGF) and builds appropriate index.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 6: MCP Server reader cache + tool migration

**Files:**
- Modify: `crates/mcp-server/Cargo.toml`
- Modify: `crates/mcp-server/src/tools.rs`

- [ ] **Step 1: Add `lru` dependency**

In the workspace `Cargo.toml`, add to `[workspace.dependencies]`:

```toml
lru = "0.12"
```

In `crates/mcp-server/Cargo.toml`, add:

```toml
lru = { workspace = true }
```

- [ ] **Step 2: Add reader cache to ProteinCopilotServer**

In `crates/mcp-server/src/tools.rs`, add import at top:

```rust
use lru::LruCache;
use std::num::NonZeroUsize;
```

Add a new type alias and constant near the other cache types:

```rust
const MAX_READER_CACHE_SIZE: usize = 20;
type ReaderCache = Arc<Mutex<LruCache<PathBuf, Arc<dyn SpectrumReader>>>>;
```

Add field to `ProteinCopilotServer`:

```rust
pub struct ProteinCopilotServer {
    tool_router: ToolRouter<Self>,
    registry: protein_copilot_search_engine::EngineRegistry,
    run_cache: RunCache,
    dia_cache: Arc<Mutex<OrderedDiaCache>>,
    /// Cached indexed spectrum readers for O(1) scan lookup.
    reader_cache: ReaderCache,
}
```

Update `new()`:

```rust
    pub fn new() -> Self {
        let mut registry = protein_copilot_search_engine::EngineRegistry::new();
        registry.register(Box::new(SimpleSearchEngine::new()));
        Self {
            tool_router: Self::tool_router(),
            registry,
            run_cache: Arc::new(Mutex::new(OrderedRunCache::new())),
            dia_cache: Arc::new(Mutex::new(OrderedDiaCache::new())),
            reader_cache: Arc::new(Mutex::new(LruCache::new(
                NonZeroUsize::new(MAX_READER_CACHE_SIZE).unwrap(),
            ))),
        }
    }
```

Add helper method:

```rust
    /// Get or create an indexed spectrum reader for a file path.
    fn get_or_create_reader(
        &self,
        path: &Path,
    ) -> Result<Arc<dyn SpectrumReader>, ErrorData> {
        let mut cache = self
            .reader_cache
            .lock()
            .map_err(|_| mcp_err(ErrorCode::INTERNAL_ERROR, "reader cache lock failed"))?;
        if let Some(reader) = cache.get(path) {
            return Ok(Arc::clone(reader));
        }
        let reader: Arc<dyn SpectrumReader> = Arc::from(
            protein_copilot_spectrum_io::create_indexed_reader(path)
                .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?,
        );
        cache.put(path.to_path_buf(), Arc::clone(&reader));
        Ok(reader)
    }
```

Note: `Arc<dyn SpectrumReader>` requires `SpectrumReader: Send + Sync`, which it already has (the trait bound includes `Send + Sync`).

- [ ] **Step 3: Migrate `get_spectrum` tool**

Change `get_spectrum` from:

```rust
        let info = protein_copilot_spectrum_io::detect_format(path)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectrum = reader
            .read_spectrum(path, input.scan_number)
```

To:

```rust
        let reader = self.get_or_create_reader(path)?;
        let spectrum = reader
            .read_spectrum(path, input.scan_number)
```

- [ ] **Step 4: Migrate `annotate_spectrum` tool**

Change the spectrum reading section from:

```rust
        let info = protein_copilot_spectrum_io::detect_format(&spectrum_file)
            .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
        let reader = protein_copilot_spectrum_io::create_reader(&info);
        let spectrum = reader
            .read_spectrum(&spectrum_file, input.scan_number)
```

To:

```rust
        let reader = self.get_or_create_reader(&spectrum_file)?;
        let spectrum = reader
            .read_spectrum(&spectrum_file, input.scan_number)
```

- [ ] **Step 5: Build and test**

Run: `cargo build --workspace 2>&1 | tail -10`
Expected: success

Run: `cargo test --workspace 2>&1 | tail -15`
Expected: all tests pass

Run: `cargo clippy --workspace -- -D warnings 2>&1 | tail -10`
Expected: 0 warnings

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml crates/mcp-server/Cargo.toml crates/mcp-server/src/tools.rs
git commit -m "feat(mcp): add reader cache for O(1) spectrum lookup

Adds LRU cache of IndexedReader instances to ProteinCopilotServer.
Migrates get_spectrum and annotate_spectrum to use cached readers.
Repeated reads of the same file now skip index building.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 7: Integration verification

**Files:**
- All crates (verification only)

- [ ] **Step 1: Full workspace build**

Run: `cargo build --workspace 2>&1 | tail -5`
Expected: success

- [ ] **Step 2: Full test suite**

Run: `cargo test --workspace 2>&1 | tail -20`
Expected: all tests pass

- [ ] **Step 3: Clippy clean**

Run: `cargo clippy --workspace -- -D warnings 2>&1 | tail -5`
Expected: 0 warnings

- [ ] **Step 4: Verify total test count**

Run: `cargo test --workspace 2>&1 | grep -oP '\d+ passed' | awk -F' ' '{sum+=$1} END{print sum" total tests"}'`
Expected: 460+ (was 444 + new indexed reader tests)

- [ ] **Step 5: Commit if any integration fixes needed**

```bash
git add -A
git commit -m "fix: integration fixes for scan index cache

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
