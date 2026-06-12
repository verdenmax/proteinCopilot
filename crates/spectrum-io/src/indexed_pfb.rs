//! Indexed PFB reader: O(1) scan lookup + cached scan metadata via a
//! reused [`crate::index::ScanIndex`].
//!
//! On [`IndexedPfbReader::open`], the footer offset table is read and each
//! record's property header is parsed (skipping peak blobs) to build a
//! `scan → (offset, rt, ms_level, isolation_window)` index. Bulk operations
//! delegate to the streaming [`crate::pfb::PfbReader`].

use std::collections::HashMap;
use std::io::{Seek, SeekFrom};
use std::path::{Path, PathBuf};

use protein_copilot_core::spectrum::{Spectrum, SpectrumSummary};

use crate::error::SpectrumIoError;
use crate::index::{IndexSource, ScanIndex, ScanMeta};
use crate::pfb::{self, PfbReader};
use crate::reader::{Ms2ScanMeta, ScanMetaInfo, SpectrumReader};

/// PFB reader backed by a [`ScanIndex`] for O(1) scan lookup.
///
/// Bulk operations (`read_all`, `read_summary`, `for_each_spectrum`) delegate
/// to the streaming [`PfbReader`]; `read_spectrum`, `list_scan_meta`,
/// `list_ms2_meta`, and `find_by_rt` use the cached index.
pub struct IndexedPfbReader {
    index: ScanIndex,
    path: PathBuf,
}

impl IndexedPfbReader {
    /// Opens a PFB file and builds a scan index from its footer + property headers.
    pub fn open(path: &Path) -> Result<Self, SpectrumIoError> {
        let index = build_pfb_index(path)?;
        Ok(Self {
            index,
            path: path.to_path_buf(),
        })
    }

    /// Returns a reference to the underlying scan index.
    pub fn index(&self) -> &ScanIndex {
        &self.index
    }
}

/// Builds a [`ScanIndex`] by reading footer offsets and each record's property
/// header (peak blobs are skipped — only seeks, no peak reads).
fn build_pfb_index(path: &Path) -> Result<ScanIndex, SpectrumIoError> {
    let mut r = crate::util::open_buffered(path)?;
    let header = pfb::read_header(&mut r, path)?;

    r.seek(SeekFrom::Start(header.addr_list_addr))
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?;
    let mut offsets = Vec::with_capacity(header.scan_num as usize);
    for _ in 0..header.scan_num {
        offsets.push(pfb::read_i64(&mut r, path, "footer offset")? as u64);
    }

    let mut entries: HashMap<u32, ScanMeta> = HashMap::with_capacity(offsets.len());
    for off in offsets {
        r.seek(SeekFrom::Start(off))
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;
        let prop = pfb::read_property_str(&mut r, path)?;
        let toks: Vec<&str> = prop.split('\t').collect();
        let scan = match toks.first().and_then(|t| t.trim().parse::<u32>().ok()) {
            Some(s) => s,
            None => continue,
        };
        let ms_level = toks
            .get(1)
            .and_then(|t| t.trim().parse::<u8>().ok())
            .unwrap_or(0);
        let rt_seconds = toks
            .get(2)
            .and_then(|t| t.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let isolation_window = pfb::isolation_window_from_tokens(&toks);
        if entries
            .insert(
                scan,
                ScanMeta {
                    offset: off,
                    rt_seconds,
                    ms_level,
                    isolation_window,
                },
            )
            .is_some()
        {
            tracing::warn!(scan, path = %path.display(), "duplicate scan number in PFB index; keeping last");
        }
    }

    Ok(ScanIndex::from_meta(entries, IndexSource::BuiltFromScan))
}

impl SpectrumReader for IndexedPfbReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        PfbReader.read_all(path)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        PfbReader.read_summary(path)
    }

    fn read_spectrum(&self, _path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        let offset = self
            .index
            .get_offset(scan)
            .ok_or_else(|| SpectrumIoError::ScanNotFound {
                path: self.path.clone(),
                scan,
            })?;
        let mut r = crate::util::open_buffered(&self.path)?;
        r.seek(SeekFrom::Start(offset))
            .map_err(|e| SpectrumIoError::IoError {
                path: self.path.clone(),
                source: e,
            })?;
        let prop = pfb::read_property_str(&mut r, &self.path)?;
        let (mz, intensity) = pfb::read_peaks(&mut r, &self.path)?;
        pfb::build_spectrum(&prop, mz, intensity, &self.path)
    }

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        PfbReader.for_each_spectrum(path, handler)
    }

    fn list_scan_meta(&self, _path: &Path) -> Result<Vec<ScanMetaInfo>, SpectrumIoError> {
        Ok(self
            .index
            .iter_meta()
            .map(|(&scan, m)| ScanMetaInfo {
                scan_number: scan,
                ms_level: m.ms_level,
                rt_min: m.rt_seconds / 60.0,
                isolation_window: m.isolation_window,
            })
            .collect())
    }

    fn list_ms2_meta(&self, _path: &Path) -> Result<Vec<Ms2ScanMeta>, SpectrumIoError> {
        Ok(self
            .index
            .iter_meta()
            .filter(|(_, m)| m.ms_level == 2)
            .map(|(&scan, m)| Ms2ScanMeta {
                scan_number: scan,
                rt_min: m.rt_seconds / 60.0,
                isolation_window: m.isolation_window,
            })
            .collect())
    }

    fn find_by_rt(
        &self,
        _path: &Path,
        rt_min: f64,
        precursor_mz: f64,
        rt_tolerance_min: f64,
    ) -> Result<Option<(u32, f64)>, SpectrumIoError> {
        Ok(self
            .index
            .find_by_rt(rt_min, precursor_mz, rt_tolerance_min))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pfb::PfbReader;
    use crate::reader::SpectrumReader;

    fn write_pfb(recs: &[(&str, Vec<f64>, Vec<f64>)]) -> (tempfile::TempDir, std::path::PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.pfb");
        let header_size: u64 = 24;
        let mut body: Vec<u8> = Vec::new();
        let mut offsets: Vec<u64> = Vec::new();
        for (prop, mz, inten) in recs {
            offsets.push(header_size + body.len() as u64);
            let pb = prop.as_bytes();
            body.extend_from_slice(&(pb.len() as i32).to_le_bytes());
            body.extend_from_slice(pb);
            body.extend_from_slice(&(mz.len() as i32).to_le_bytes());
            for &m in mz {
                body.extend_from_slice(&m.to_le_bytes());
            }
            for &v in inten {
                body.extend_from_slice(&v.to_le_bytes());
            }
        }
        let addr_list_addr = header_size + body.len() as u64;
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(&0i32.to_le_bytes());
        out.extend_from_slice(&0i32.to_le_bytes());
        out.extend_from_slice(&0i32.to_le_bytes());
        out.extend_from_slice(&(addr_list_addr as i64).to_le_bytes());
        out.extend_from_slice(&(recs.len() as i32).to_le_bytes());
        out.extend_from_slice(&body);
        for &off in &offsets {
            out.extend_from_slice(&(off as i64).to_le_bytes());
        }
        std::fs::write(&path, &out).unwrap();
        (dir, path)
    }

    fn sample_recs() -> Vec<(&'static str, Vec<f64>, Vec<f64>)> {
        vec![
            ("1\t1\t6.0\tFTMS", vec![100.0, 200.0], vec![10.0, 20.0]),
            (
                "2\t2\t12.0\tFTMS\t2\t1000.0\t50\t501.0\tHCD\t1\t2.0\t27.0\t501.25",
                vec![150.0, 250.0, 350.0],
                vec![5.0, 15.0, 25.0],
            ),
            (
                "3\t2\t18.0\tFTMS\t3\t1200.0\t50\t600.0\tHCD\t1\t2.0\t27.0\t600.5",
                vec![120.0, 220.0],
                vec![7.0, 8.0],
            ),
        ]
    }

    #[test]
    fn open_builds_index() {
        let (_d, p) = write_pfb(&sample_recs());
        let reader = IndexedPfbReader::open(&p).unwrap();
        assert_eq!(reader.index().len(), 3);
    }

    #[test]
    fn read_spectrum_matches_standard() {
        let (_d, p) = write_pfb(&sample_recs());
        let indexed = IndexedPfbReader::open(&p).unwrap();
        for scan in [1u32, 2, 3] {
            let a = indexed.read_spectrum(&p, scan).unwrap();
            let b = PfbReader.read_spectrum(&p, scan).unwrap();
            assert_eq!(a.scan_number, b.scan_number);
            assert_eq!(a.num_peaks(), b.num_peaks());
            assert_eq!(a.ms_level, b.ms_level);
        }
    }

    #[test]
    fn list_scan_meta_from_cache() {
        let (_d, p) = write_pfb(&sample_recs());
        let reader = IndexedPfbReader::open(&p).unwrap();
        let mut metas = reader.list_scan_meta(&p).unwrap();
        metas.sort_by_key(|m| m.scan_number);
        assert_eq!(metas.len(), 3);
        assert_eq!(metas[0].ms_level, 1);
        assert_eq!(metas[1].ms_level, 2);
        assert!((metas[0].rt_min - 0.1).abs() < 1e-9);
    }

    #[test]
    fn list_ms2_meta_filters() {
        let (_d, p) = write_pfb(&sample_recs());
        let reader = IndexedPfbReader::open(&p).unwrap();
        let metas = reader.list_ms2_meta(&p).unwrap();
        assert_eq!(metas.len(), 2);
    }

    #[test]
    fn find_by_rt_uses_index() {
        let (_d, p) = write_pfb(&sample_recs());
        let reader = IndexedPfbReader::open(&p).unwrap();
        let hit = reader.find_by_rt(&p, 0.2, 501.25, 0.05).unwrap();
        assert_eq!(hit.map(|(scan, _)| scan), Some(2));
    }

    #[test]
    fn read_spectrum_not_found() {
        let (_d, p) = write_pfb(&sample_recs());
        let reader = IndexedPfbReader::open(&p).unwrap();
        let err = reader.read_spectrum(&p, 999).unwrap_err();
        assert!(matches!(
            err,
            SpectrumIoError::ScanNotFound { scan: 999, .. }
        ));
    }
}
