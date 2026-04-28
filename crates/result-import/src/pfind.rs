//! Parser for pFind .spectra result files (skeleton).
//!
//! pFind results include scan numbers directly, so scan matching
//! is not required. Implementation pending sample file availability.

use std::path::Path;

use crate::unimod::UnimodDb;
use crate::{ImportedPsm, ResultImportError, ResultParser};

/// Parser for pFind .spectra result files.
pub struct PFindParser;

impl ResultParser for PFindParser {
    fn parse(
        &self,
        path: &Path,
        _unimod: &UnimodDb,
    ) -> Result<Vec<ImportedPsm>, ResultImportError> {
        let _span = tracing::info_span!("parse_pfind",
            file = %path.display(),
        ).entered();

        if !path.exists() {
            return Err(ResultImportError::FileNotFound {
                path: path.to_path_buf(),
            });
        }

        // TODO: Implement when sample .spectra file is available.
        // pFind results include scan numbers, so scan matching is not needed.
        Err(ResultImportError::Other(
            "pFind .spectra parser not yet implemented — awaiting sample file".to_string(),
        ))
    }
}

/// Detect whether a file is a pFind .spectra file.
pub fn detect(path: &Path) -> bool {
    path.extension().is_some_and(|e| e == "spectra")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_spectra_extension() {
        assert!(detect(Path::new("result.spectra")));
        assert!(!detect(Path::new("result.json")));
    }

    #[test]
    fn parse_returns_not_implemented() {
        let db = UnimodDb::builtin();
        let result = PFindParser.parse(Path::new("/nonexistent.spectra"), &db);
        assert!(result.is_err());
    }
}
