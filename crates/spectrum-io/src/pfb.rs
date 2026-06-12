//! PFB (pXtract3 / pParse2+ binary) spectrum reader.
//!
//! Little-endian binary layout:
//! - Header (24 bytes): 3×i32 (reserved) + i64 addr_list_addr + i32 scan_num
//! - scan_num records: i32 prop_len + prop_str (UTF-8, '\t'-separated) +
//!   i32 peak_num + f64×peak_num m/z + f64×peak_num intensity
//! - Footer @ addr_list_addr: i64×scan_num record offsets
//!
//! property_str (tab-separated, by position):
//! `[0]`Scan `[1]`MsType(1=MS1,2=MS2) `[2]`RT(seconds) `[3]`InstrumentType;
//! MS2 adds `[4]`Charge `[5]`MH+ `[6]`IonInjectionTime `[7]`ActivationCenter
//! `[8]`ActivationType `[9]`PrecursorScan `[10]`ActivationWindow `[11]`NCE
//! `[12]`monoisotopicMz.

use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use protein_copilot_core::spectrum::{
    IsolationWindow, MsLevel, PrecursorInfo, Spectrum, SpectrumFormat, SpectrumSummary,
};

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Sanity bound to avoid huge allocations from a corrupt `peak_num`.
///
/// 10M peaks is ~100× above any real spectrum (typical scans have far fewer
/// than 100k peaks) and caps the eager `vec![0u8; n*8]` allocation at ~160 MB
/// (two arrays) rather than ~1.6 GB, mitigating denial-of-service from crafted
/// or corrupt files.
const MAX_PEAKS_PER_SCAN: usize = 10_000_000;

/// Sanity bound for a single property string length.
const MAX_PROP_LEN: usize = 100_000_000;

/// Reader for PFB binary spectrum files.
pub struct PfbReader;

/// Decoded PFB header.
pub(crate) struct PfbHeader {
    pub addr_list_addr: u64,
    pub scan_num: u32,
}

fn parse_err(path: &Path, detail: String) -> SpectrumIoError {
    SpectrumIoError::ParseError {
        path: path.to_path_buf(),
        line: 0,
        detail,
    }
}

fn read_buf(
    r: &mut impl Read,
    buf: &mut [u8],
    path: &Path,
    what: &str,
) -> Result<(), SpectrumIoError> {
    r.read_exact(buf)
        .map_err(|e| parse_err(path, format!("short read while reading {what}: {e}")))
}

fn read_i32(r: &mut impl Read, path: &Path, what: &str) -> Result<i32, SpectrumIoError> {
    let mut b = [0u8; 4];
    read_buf(r, &mut b, path, what)?;
    Ok(i32::from_le_bytes(b))
}

pub(crate) fn read_i64(r: &mut impl Read, path: &Path, what: &str) -> Result<i64, SpectrumIoError> {
    let mut b = [0u8; 8];
    read_buf(r, &mut b, path, what)?;
    Ok(i64::from_le_bytes(b))
}

fn bytes_to_f64_vec(buf: &[u8]) -> Vec<f64> {
    buf.chunks_exact(8)
        .map(|c| {
            let mut a = [0u8; 8];
            a.copy_from_slice(c);
            f64::from_le_bytes(a)
        })
        .collect()
}

/// Reads the 24-byte header from the current (start) position.
pub(crate) fn read_header(r: &mut impl Read, path: &Path) -> Result<PfbHeader, SpectrumIoError> {
    let _e1 = read_i32(r, path, "empty_1")?;
    let _e2 = read_i32(r, path, "empty_2")?;
    let _e3 = read_i32(r, path, "empty_3")?;
    let addr = read_i64(r, path, "addr_list_addr")?;
    let scan_num = read_i32(r, path, "scan_num")?;
    if addr < 0 || scan_num < 0 {
        return Err(parse_err(
            path,
            format!("invalid header: addr_list_addr={addr} scan_num={scan_num}"),
        ));
    }
    let addr_list_addr = addr as u64;
    let scan_num = scan_num as u32;
    let file_size = std::fs::metadata(path)
        .map_err(|e| SpectrumIoError::IoError {
            path: path.to_path_buf(),
            source: e,
        })?
        .len();
    // Integrity gate: a well-formed PFB ends with exactly `scan_num` i64 footer
    // offsets, so the footer is the last bytes of the file. This also bounds
    // `scan_num` against the real file size. (Assumes no trailing data appended
    // after the footer.)
    let expected = (scan_num as u64)
        .checked_mul(8)
        .and_then(|footer| addr_list_addr.checked_add(footer));
    if expected != Some(file_size) {
        return Err(parse_err(
            path,
            format!(
                "header invariant violated: addr_list_addr({addr_list_addr}) + scan_num({scan_num})*8 != file_size({file_size})"
            ),
        ));
    }
    Ok(PfbHeader {
        addr_list_addr,
        scan_num,
    })
}

/// Reads an i32 length-prefixed UTF-8 property string (trailing NUL stripped).
pub(crate) fn read_property_str(r: &mut impl Read, path: &Path) -> Result<String, SpectrumIoError> {
    let prop_len = read_i32(r, path, "property_str_len")?;
    if prop_len < 0 {
        return Err(parse_err(
            path,
            format!("negative property_str_len {prop_len}"),
        ));
    }
    if prop_len as usize > MAX_PROP_LEN {
        return Err(parse_err(
            path,
            format!("property_str_len {prop_len} exceeds sanity bound"),
        ));
    }
    let mut buf = vec![0u8; prop_len as usize];
    read_buf(r, &mut buf, path, "property_str")?;
    Ok(String::from_utf8_lossy(&buf)
        .trim_end_matches('\0')
        .to_string())
}

/// Reads `peak_num` (i32) then the two parallel f64 peak arrays.
pub(crate) fn read_peaks(
    r: &mut impl Read,
    path: &Path,
) -> Result<(Vec<f64>, Vec<f64>), SpectrumIoError> {
    let peak_num = read_i32(r, path, "peak_num")?;
    if peak_num < 0 {
        return Err(parse_err(path, format!("negative peak_num {peak_num}")));
    }
    let n = peak_num as usize;
    if n > MAX_PEAKS_PER_SCAN {
        return Err(parse_err(
            path,
            format!("peak_num {n} exceeds sanity bound"),
        ));
    }
    let mut mz_buf = vec![0u8; n * 8];
    read_buf(r, &mut mz_buf, path, "mz array")?;
    let mut in_buf = vec![0u8; n * 8];
    read_buf(r, &mut in_buf, path, "intensity array")?;
    Ok((bytes_to_f64_vec(&mz_buf), bytes_to_f64_vec(&in_buf)))
}

/// Parses the isolation window from a property string's tab-separated tokens.
///
/// Uses `[7]` ActivationCenter and `[10]` ActivationWindow, returning
/// `(target_mz, lower_offset, upper_offset)` with symmetric half-window
/// offsets. Returns `None` when either field is absent or the window is not
/// positive. Shared by `build_spectrum` and the indexed reader's index builder
/// to keep index metadata and spectra coherent.
pub(crate) fn isolation_window_from_tokens(toks: &[&str]) -> Option<(f64, f64, f64)> {
    let center = toks.get(7).and_then(|t| t.trim().parse::<f64>().ok());
    let window = toks.get(10).and_then(|t| t.trim().parse::<f64>().ok());
    match (center, window) {
        (Some(c), Some(w)) if w > 0.0 => Some((c, w / 2.0, w / 2.0)),
        _ => None,
    }
}

/// Builds a validated `Spectrum` from a property string + peak arrays.
pub(crate) fn build_spectrum(
    property_str: &str,
    mut mz: Vec<f64>,
    mut intensity: Vec<f64>,
    path: &Path,
) -> Result<Spectrum, SpectrumIoError> {
    let toks: Vec<&str> = property_str.split('\t').collect();
    let scan = toks
        .first()
        .and_then(|t| t.trim().parse::<u32>().ok())
        .ok_or_else(|| parse_err(path, format!("missing/invalid Scan in '{property_str}'")))?;
    let ms_type = toks
        .get(1)
        .and_then(|t| t.trim().parse::<u8>().ok())
        .ok_or_else(|| parse_err(path, format!("scan {scan}: missing/invalid MsType")))?;
    let ms_level = match ms_type {
        1 => MsLevel::MS1,
        2 => MsLevel::MS2,
        n => MsLevel::Other(n),
    };
    let rt_min = toks
        .get(2)
        .and_then(|t| t.trim().parse::<f64>().ok())
        .unwrap_or(0.0)
        / 60.0;

    let precursors = if ms_level == MsLevel::MS2 {
        let charge = toks
            .get(4)
            .and_then(|t| t.trim().parse::<i32>().ok())
            .filter(|&c| c != 0);
        let activation_center = toks.get(7).and_then(|t| t.trim().parse::<f64>().ok());
        let precursor_scan = toks.get(9).and_then(|t| t.trim().parse::<u32>().ok());
        let mono_mz = toks.get(12).and_then(|t| t.trim().parse::<f64>().ok());
        let mz_val = mono_mz.or(activation_center).unwrap_or(0.0);
        let isolation_window =
            isolation_window_from_tokens(&toks).map(|(target_mz, lower_offset, upper_offset)| {
                IsolationWindow {
                    target_mz,
                    lower_offset,
                    upper_offset,
                }
            });
        vec![PrecursorInfo {
            mz: mz_val,
            charge,
            intensity: None,
            isolation_window,
            source_scan: precursor_scan,
        }]
    } else {
        Vec::new()
    };

    crate::util::sort_peaks_by_mz(&mut mz, &mut intensity);

    Spectrum::new(scan, ms_level, rt_min, precursors, mz, intensity).map_err(|e| {
        SpectrumIoError::ValidationError {
            scan,
            detail: e.to_string(),
        }
    })
}

impl PfbReader {
    /// Streams every record sequentially from offset 24 (header size).
    fn stream<F>(&self, path: &Path, mut handler: F) -> Result<u32, SpectrumIoError>
    where
        F: FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    {
        let mut r = crate::util::open_buffered(path)?;
        let header = read_header(&mut r, path)?;
        let mut count = 0u32;
        for _ in 0..header.scan_num {
            let prop = read_property_str(&mut r, path)?;
            let (mz, intensity) = read_peaks(&mut r, path)?;
            let spec = build_spectrum(&prop, mz, intensity, path)?;
            count += 1;
            if !handler(spec)? {
                break;
            }
        }
        Ok(count)
    }
}

impl SpectrumReader for PfbReader {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError> {
        let mut out = Vec::new();
        self.stream(path, |s| {
            out.push(s);
            Ok(true)
        })?;
        Ok(out)
    }

    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError> {
        let mut acc = crate::util::SummaryAccumulator::new();
        self.stream(path, |s| {
            acc.observe(&s);
            Ok(true)
        })?;
        acc.into_summary(path, SpectrumFormat::Pfb)
    }

    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError> {
        let mut r = crate::util::open_buffered(path)?;
        let header = read_header(&mut r, path)?;
        r.seek(SeekFrom::Start(header.addr_list_addr))
            .map_err(|e| SpectrumIoError::IoError {
                path: path.to_path_buf(),
                source: e,
            })?;
        let mut offsets = Vec::with_capacity(header.scan_num as usize);
        for _ in 0..header.scan_num {
            offsets.push(read_i64(&mut r, path, "footer offset")? as u64);
        }
        for off in offsets {
            r.seek(SeekFrom::Start(off))
                .map_err(|e| SpectrumIoError::IoError {
                    path: path.to_path_buf(),
                    source: e,
                })?;
            let prop = read_property_str(&mut r, path)?;
            let this_scan = prop
                .split('\t')
                .next()
                .and_then(|t| t.trim().parse::<u32>().ok());
            if this_scan == Some(scan) {
                let (mz, intensity) = read_peaks(&mut r, path)?;
                return build_spectrum(&prop, mz, intensity, path);
            }
        }
        Err(SpectrumIoError::ScanNotFound {
            path: path.to_path_buf(),
            scan,
        })
    }

    fn for_each_spectrum(
        &self,
        path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>,
    ) -> Result<u32, SpectrumIoError> {
        self.stream(path, handler)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::SpectrumReader;
    use protein_copilot_core::spectrum::{MsLevel, SpectrumFormat};

    /// Writes a PFB file with the given records; returns (tempdir, path).
    /// Keep the returned TempDir alive for the duration of the test.
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
                vec![350.0, 150.0, 250.0],
                vec![25.0, 5.0, 15.0],
            ),
            (
                "3\t2\t18.0\tFTMS\t3\t1200.0\t50\t600.0\tHCD\t1\t2.0\t27.0\t600.5",
                vec![120.0, 220.0],
                vec![7.0, 8.0],
            ),
        ]
    }

    #[test]
    fn read_all_parses_records() {
        let (_d, p) = write_pfb(&sample_recs());
        let spectra = PfbReader.read_all(&p).unwrap();
        assert_eq!(spectra.len(), 3);
        assert_eq!(spectra[0].scan_number, 1);
        assert_eq!(spectra[0].ms_level, MsLevel::MS1);
        assert_eq!(spectra[1].ms_level, MsLevel::MS2);
        assert_eq!(
            spectra.iter().map(|s| s.num_peaks()).collect::<Vec<_>>(),
            vec![2, 3, 2]
        );
    }

    #[test]
    fn ms1_has_no_precursor_and_rt_seconds_to_minutes() {
        let (_d, p) = write_pfb(&sample_recs());
        let s = &PfbReader.read_all(&p).unwrap()[0];
        assert!(s.precursors.is_empty());
        assert!((s.retention_time_min - 0.1).abs() < 1e-9);
    }

    #[test]
    fn ms2_precursor_mapping() {
        let (_d, p) = write_pfb(&sample_recs());
        let s = &PfbReader.read_all(&p).unwrap()[1];
        assert_eq!(s.precursors.len(), 1);
        let pr = &s.precursors[0];
        assert!((pr.mz - 501.25).abs() < 1e-9);
        assert_eq!(pr.charge, Some(2));
        assert_eq!(pr.source_scan, Some(1));
        let iw = pr.isolation_window.as_ref().unwrap();
        assert!((iw.target_mz - 501.0).abs() < 1e-9);
        assert!((iw.lower_offset - 1.0).abs() < 1e-9);
        assert!((iw.upper_offset - 1.0).abs() < 1e-9);
        for w in s.mz_array.windows(2) {
            assert!(w[0] <= w[1]);
        }
        assert_eq!(s.intensity_array, vec![5.0, 15.0, 25.0]);
    }

    #[test]
    fn read_summary_counts() {
        let (_d, p) = write_pfb(&sample_recs());
        let sum = PfbReader.read_summary(&p).unwrap();
        assert_eq!(sum.total_spectra, 3);
        assert_eq!(sum.ms1_count, 1);
        assert_eq!(sum.ms2_count, 2);
        assert_eq!(sum.format, SpectrumFormat::Pfb);
    }

    #[test]
    fn read_spectrum_by_scan() {
        let (_d, p) = write_pfb(&sample_recs());
        let s = PfbReader.read_spectrum(&p, 2).unwrap();
        assert_eq!(s.scan_number, 2);
        assert_eq!(s.num_peaks(), 3);
        assert_eq!(s.precursors[0].charge, Some(2));
    }

    #[test]
    fn read_spectrum_not_found() {
        let (_d, p) = write_pfb(&sample_recs());
        let err = PfbReader.read_spectrum(&p, 999).unwrap_err();
        assert!(matches!(
            err,
            SpectrumIoError::ScanNotFound { scan: 999, .. }
        ));
    }

    #[test]
    fn truncated_file_errors() {
        let (_d, p) = write_pfb(&sample_recs());
        let bytes = std::fs::read(&p).unwrap();
        std::fs::write(&p, &bytes[..28]).unwrap();
        let err = PfbReader.read_all(&p).unwrap_err();
        assert!(matches!(err, SpectrumIoError::ParseError { .. }));
    }

    #[test]
    fn ms2_with_minimal_fields_falls_back() {
        let recs = vec![(
            "5\t2\t30.0\tFTMS\t2\t900.0\t50\t450.0\tHCD\t1",
            vec![100.0, 200.0],
            vec![1.0, 2.0],
        )];
        let (_d, p) = write_pfb(&recs);
        let s = &PfbReader.read_all(&p).unwrap()[0];
        let pr = &s.precursors[0];
        assert!((pr.mz - 450.0).abs() < 1e-9);
        assert!(pr.isolation_window.is_none());
        assert_eq!(pr.charge, Some(2));
        assert_eq!(pr.source_scan, Some(1));
    }

    #[test]
    fn detect_and_create_reader_for_pfb() {
        let (_d, p) = write_pfb(&sample_recs());
        let info = crate::detect_format(&p).unwrap();
        assert_eq!(info.format, SpectrumFormat::Pfb);
        assert_eq!(crate::create_reader(&info).read_all(&p).unwrap().len(), 3);
        let ireader = crate::create_indexed_reader(&p).unwrap();
        assert_eq!(ireader.read_spectrum(&p, 3).unwrap().scan_number, 3);
    }

    #[test]
    fn spectrum_format_display_pfb() {
        assert_eq!(SpectrumFormat::Pfb.to_string(), "pfb");
    }

    #[test]
    fn ms2_without_precursor_mz_errors() {
        // 7 tab fields: Scan, MsType, RT, Instrument, Charge, MH+, InjTime —
        // no ActivationCenter[7] and no monoisotopicMz[12] → precursor m/z unknown.
        let recs = vec![(
            "7\t2\t10.0\tFTMS\t2\t900.0\t50",
            vec![100.0, 200.0],
            vec![1.0, 2.0],
        )];
        let (_d, p) = write_pfb(&recs);
        let err = PfbReader.read_all(&p).unwrap_err();
        assert!(matches!(
            err,
            SpectrumIoError::ValidationError { scan: 7, .. }
        ));
    }

    #[test]
    fn read_peaks_rejects_peak_num_above_bound() {
        // 11M peaks is below the historical 100M bound but above the hardened
        // 10M ceiling. It must be rejected up-front with a sanity-bound error
        // (no ~160 MB eager allocation, no read attempt).
        let peak_num: i32 = 11_000_000;
        let mut cursor = std::io::Cursor::new(peak_num.to_le_bytes().to_vec());
        let err = read_peaks(&mut cursor, std::path::Path::new("synthetic.pfb")).unwrap_err();
        match err {
            SpectrumIoError::ParseError { detail, .. } => assert!(
                detail.contains("exceeds sanity bound"),
                "expected sanity-bound rejection, got: {detail}"
            ),
            other => panic!("expected ParseError, got {other:?}"),
        }
        // Pin the hardened ceiling so the bound cannot silently regress.
        assert_eq!(MAX_PEAKS_PER_SCAN, 10_000_000);
    }
}
