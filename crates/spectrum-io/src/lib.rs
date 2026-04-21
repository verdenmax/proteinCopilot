//! # ProteinCopilot Spectrum I/O
//!
//! Library crate for parsing mass spectrometry spectrum files into
//! the shared [`protein_copilot_core::spectrum`] types.
//!
//! Supported formats:
//! - **MGF** (Mascot Generic Format) — text-based, widely used
//! - **mzML** — PSI standard XML format
//!
//! # Usage
//!
//! ```no_run
//! use std::path::Path;
//! use protein_copilot_spectrum_io::{detect_format, create_reader};
//!
//! let path = Path::new("data/sample.mgf");
//! let file_info = detect_format(path).unwrap();
//! let reader = create_reader(&file_info);
//! let summary = reader.read_summary(path).unwrap();
//! println!("Total spectra: {}", summary.total_spectra);
//! ```

pub mod disk_cache;
pub mod error;
pub mod index;
pub mod indexed_mgf;
pub mod indexed_mzml;
pub mod mgf;
pub mod mzml;
pub mod reader;
mod util;

pub use error::SpectrumIoError;
pub use indexed_mgf::IndexedMgfReader;
pub use indexed_mzml::IndexedMzMLReader;
pub use reader::SpectrumReader;

use std::path::Path;

use protein_copilot_core::spectrum::{SpectrumFileInfo, SpectrumFormat};

/// Detects the spectrum file format from the file path.
///
/// Format detection uses file extension:
/// - `.mgf` → [`SpectrumFormat::Mgf`]
/// - `.mzml` → [`SpectrumFormat::MzML`]
///
/// Also verifies the file exists and records its size.
pub fn detect_format(path: &Path) -> Result<SpectrumFileInfo, SpectrumIoError> {
    if !path.exists() {
        return Err(SpectrumIoError::FileNotFound {
            path: path.to_path_buf(),
        });
    }

    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase());

    let format = match ext.as_deref() {
        Some("mgf") => SpectrumFormat::Mgf,
        Some("mzml") => SpectrumFormat::MzML,
        _ => {
            return Err(SpectrumIoError::UnknownFormat {
                path: path.to_path_buf(),
            });
        }
    };

    let metadata = std::fs::metadata(path).map_err(|e| SpectrumIoError::IoError {
        path: path.to_path_buf(),
        source: e,
    })?;

    Ok(SpectrumFileInfo {
        path: path.to_string_lossy().to_string(),
        format,
        file_size_bytes: metadata.len(),
    })
}

/// Creates the appropriate [`SpectrumReader`] for the given file info.
///
/// Returns a boxed trait object that can read spectra from the file.
pub fn create_reader(info: &SpectrumFileInfo) -> Box<dyn SpectrumReader> {
    match info.format {
        SpectrumFormat::Mgf => Box::new(mgf::MgfReader),
        SpectrumFormat::MzML => Box::new(mzml::MzMLReader),
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn create_temp_file(ext: &str) -> tempfile::NamedTempFile {
        let mut builder = tempfile::Builder::new();
        builder.suffix(ext);
        let mut f = builder.tempfile().unwrap();
        writeln!(f, "test content").unwrap();
        f
    }

    #[test]
    fn detect_format_mgf() {
        let f = create_temp_file(".mgf");
        let info = detect_format(f.path()).unwrap();
        assert_eq!(info.format, SpectrumFormat::Mgf);
        assert!(info.file_size_bytes > 0);
    }

    #[test]
    fn detect_format_mzml() {
        let f = create_temp_file(".mzml");
        let info = detect_format(f.path()).unwrap();
        assert_eq!(info.format, SpectrumFormat::MzML);
    }

    #[test]
    fn detect_format_case_insensitive() {
        let f = create_temp_file(".MgF");
        let info = detect_format(f.path()).unwrap();
        assert_eq!(info.format, SpectrumFormat::Mgf);
    }

    #[test]
    fn detect_format_unknown_extension() {
        let f = create_temp_file(".raw");
        let err = detect_format(f.path()).unwrap_err();
        assert!(err.to_string().contains("detect format"));
    }

    #[test]
    fn detect_format_file_not_found() {
        let err = detect_format(Path::new("/nonexistent/file.mgf")).unwrap_err();
        assert!(err.to_string().contains("not found"));
    }

    #[test]
    fn create_reader_returns_correct_type() {
        let mgf_info = SpectrumFileInfo {
            path: "test.mgf".to_string(),
            format: SpectrumFormat::Mgf,
            file_size_bytes: 100,
        };
        let _reader = create_reader(&mgf_info);

        let mzml_info = SpectrumFileInfo {
            path: "test.mzml".to_string(),
            format: SpectrumFormat::MzML,
            file_size_bytes: 100,
        };
        let _reader = create_reader(&mzml_info);
    }

    #[test]
    fn create_indexed_reader_mzml() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mzml");
        let reader = create_indexed_reader(&path).unwrap();
        let spec = reader.read_spectrum(&path, 5).unwrap();
        assert_eq!(spec.scan_number, 5);
    }

    #[test]
    fn create_indexed_reader_mgf() {
        let path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/small.mgf");
        let reader = create_indexed_reader(&path).unwrap();
        let spec = reader.read_spectrum(&path, 3).unwrap();
        assert_eq!(spec.scan_number, 3);
    }
}
