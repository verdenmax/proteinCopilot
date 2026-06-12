# PFB 二进制谱图格式支持 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 `spectrum-io` 新增只读 PFB（pXtract3/pParse2+ 二进制）谱图读取器（`PfbReader` + `IndexedPfbReader`），接入 `SpectrumFormat::Pfb` 与 detect/create/indexed，并放宽 XIC 守卫使 PFB 可做 XIC。

**Architecture:** PFB 是小端二进制：24B 头部（含 footer 地址与 scan_num）+ 顺序记录（property_str + 双 f64 峰数组）+ footer（记录偏移表）。`PfbReader` 顺序流式实现批量读取、footer 随机读；`IndexedPfbReader` 打开时扫描属性头（跳过峰 blob）构建复用的 `ScanIndex`（scan→偏移 + 元数据），实现 O(1) 随机读与缓存式 scan 元数据/RT 查找。

**Tech Stack:** Rust（spectrum-io / core / xic crates）、`SpectrumReader` trait、`crate::index::{ScanIndex, ScanMeta, IndexSource}`、`crate::util`（open_buffered/sort_peaks_by_mz/SummaryAccumulator）；测试用 `cargo test` + 测试内合成 `.pfb` 字节。

> 设计文档：`docs/superpowers/specs/2026-06-12-pfb-format-support-design.md`
> 格式已用真实样本核实：小端；头部 24B；记录 `prop_len(i32)+prop_str+peak_num(i32)+mz[f64×n]+inten[f64×n]`；footer=`i64×scan_num` 记录偏移；RT=秒（/60→分钟）；强度=f64；MS1 4 字段 / MS2 13 字段。

---

## 质量门（每个 Task 完成前对所改 crate 执行）

```bash
cargo fmt --check -p <crate>
cargo clippy -p <crate> --all-targets          # 不得有源于本次改动的新警告
cargo test -p <crate>
```
> 注意：`crates/xic/src/extract.rs` 与 `crates/mcp-server` 有**预先存在**的死代码/告警；`cargo clippy --workspace -- -D warnings` 会因此失败，与本特性无关。按「所改 crate 无新增告警」验收。

## 文件结构

| 文件 | 职责 | 变更 |
|------|------|------|
| `crates/core/src/spectrum.rs` | `SpectrumFormat::Pfb` + `Display` | 修改 |
| `crates/spectrum-io/src/pfb.rs` | PFB 解析原语 + `PfbReader`（流式 + footer 随机读） | 创建 |
| `crates/spectrum-io/src/indexed_pfb.rs` | `IndexedPfbReader`（ScanIndex 索引） | 创建 |
| `crates/spectrum-io/src/lib.rs` | 模块声明/导出 + detect_format + create_reader + create_indexed_reader | 修改 |
| `crates/xic/src/extract.rs` | 放宽 2 处格式守卫为 `MzML | Pfb` | 修改 |

## 任务总览

- **Task 1**：`SpectrumFormat::Pfb`（core+Display）+ `pfb.rs`（`PfbReader`）+ `lib.rs` 接入（detect_format/create_reader/create_indexed_reader→暂用 PfbReader）。一个内聚提交，使 PFB 经公共 API 完整可读（非索引）。
- **Task 2**：`indexed_pfb.rs`（`IndexedPfbReader`，复用 `ScanIndex::from_meta`）+ 把 `create_indexed_reader` 的 Pfb 臂切到 `IndexedPfbReader` + 导出。
- **Task 3**：放宽 `xic/extract.rs` 2 处守卫支持 PFB，端到端测试（合成 `.pfb` 过守卫并出 XIC）。

---

### Task 1: `SpectrumFormat::Pfb` + `pfb.rs`（`PfbReader`）+ lib 接入

**Files:**
- Modify: `crates/core/src/spectrum.rs`（枚举 + Display）
- Create: `crates/spectrum-io/src/pfb.rs`
- Modify: `crates/spectrum-io/src/lib.rs`（`pub mod pfb;` + detect_format + create_reader + create_indexed_reader）
- Test: `crates/spectrum-io/src/pfb.rs` 内 `#[cfg(test)] mod tests`

> 本任务把「枚举变体 + 读取器 + 接线」合为一个内聚提交：因 `PfbReader::read_summary` 需要 `SpectrumFormat::Pfb`，且新增枚举变体会破坏 `lib.rs` 两处穷尽 `match`，三者必须同时落地才能编译。

- [ ] **Step 1: 写失败测试 + 声明模块**

在 `crates/spectrum-io/src/lib.rs` 模块声明区（`pub mod mzml;` 附近）加：
```rust
pub mod pfb;
```

创建 `crates/spectrum-io/src/pfb.rs`，先放模块文档 + 测试（实现与 `SpectrumFormat::Pfb` 均未定义→编译失败）：
```rust
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
        assert!((s.retention_time_min - 0.1).abs() < 1e-9); // 6 sec / 60
    }

    #[test]
    fn ms2_precursor_mapping() {
        let (_d, p) = write_pfb(&sample_recs());
        let s = &PfbReader.read_all(&p).unwrap()[1];
        assert_eq!(s.precursors.len(), 1);
        let pr = &s.precursors[0];
        assert!((pr.mz - 501.25).abs() < 1e-9); // monoisotopicMz
        assert_eq!(pr.charge, Some(2));
        assert_eq!(pr.source_scan, Some(1));
        let iw = pr.isolation_window.as_ref().unwrap();
        assert!((iw.target_mz - 501.0).abs() < 1e-9);
        assert!((iw.lower_offset - 1.0).abs() < 1e-9); // ActivationWindow/2
        assert!((iw.upper_offset - 1.0).abs() < 1e-9);
        for w in s.mz_array.windows(2) {
            assert!(w[0] <= w[1]); // sorted despite unsorted input
        }
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
        assert!(matches!(err, SpectrumIoError::ScanNotFound { scan: 999, .. }));
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
        assert!((pr.mz - 450.0).abs() < 1e-9); // fallback to ActivationCenter
        assert!(pr.isolation_window.is_none()); // window missing
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
}
```

- [ ] **Step 2: 跑测试确认编译失败**

Run: `cargo test -p protein-copilot-spectrum-io pfb::`
Expected: 编译错误（`cannot find ... PfbReader`、`no variant ... Pfb`）。

- [ ] **Step 3a: 加 `SpectrumFormat::Pfb` + Display（core）**

`crates/core/src/spectrum.rs`，枚举在 `Mgf,` 之后加：
```rust
    /// pXtract3 / pParse2+ binary PFB format.
    Pfb,
```
`Display` impl 在 `SpectrumFormat::Mgf => write!(f, "mgf"),` 之后加：
```rust
            SpectrumFormat::Pfb => write!(f, "pfb"),
```

- [ ] **Step 3b: 写 `pfb.rs` 实现（插入到模块文档之后、测试模块之前）**

```rust
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use protein_copilot_core::spectrum::{
    IsolationWindow, MsLevel, PrecursorInfo, Spectrum, SpectrumFormat, SpectrumSummary,
};

use crate::error::SpectrumIoError;
use crate::reader::SpectrumReader;

/// Sanity bound to avoid huge allocations from a corrupt `peak_num`.
const MAX_PEAKS_PER_SCAN: usize = 100_000_000;

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
    Ok(PfbHeader {
        addr_list_addr: addr as u64,
        scan_num: scan_num as u32,
    })
}

/// Reads an i32 length-prefixed UTF-8 property string (trailing NUL stripped).
pub(crate) fn read_property_str(r: &mut impl Read, path: &Path) -> Result<String, SpectrumIoError> {
    let prop_len = read_i32(r, path, "property_str_len")?;
    if prop_len < 0 {
        return Err(parse_err(path, format!("negative property_str_len {prop_len}")));
    }
    let mut buf = vec![0u8; prop_len as usize];
    read_buf(r, &mut buf, path, "property_str")?;
    Ok(String::from_utf8_lossy(&buf)
        .trim_end_matches('\0')
        .to_string())
}

/// Reads `peak_num` (i32) then the two parallel f64 peak arrays.
pub(crate) fn read_peaks(r: &mut impl Read, path: &Path) -> Result<(Vec<f64>, Vec<f64>), SpectrumIoError> {
    let peak_num = read_i32(r, path, "peak_num")?;
    if peak_num < 0 {
        return Err(parse_err(path, format!("negative peak_num {peak_num}")));
    }
    let n = peak_num as usize;
    if n > MAX_PEAKS_PER_SCAN {
        return Err(parse_err(path, format!("peak_num {n} exceeds sanity bound")));
    }
    let mut mz_buf = vec![0u8; n * 8];
    read_buf(r, &mut mz_buf, path, "mz array")?;
    let mut in_buf = vec![0u8; n * 8];
    read_buf(r, &mut in_buf, path, "intensity array")?;
    Ok((bytes_to_f64_vec(&mz_buf), bytes_to_f64_vec(&in_buf)))
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
        let activation_window = toks.get(10).and_then(|t| t.trim().parse::<f64>().ok());
        let mono_mz = toks.get(12).and_then(|t| t.trim().parse::<f64>().ok());
        let mz_val = mono_mz.or(activation_center).unwrap_or(0.0);
        let isolation_window = match (activation_center, activation_window) {
            (Some(center), Some(win)) if win > 0.0 => Some(IsolationWindow {
                target_mz: center,
                lower_offset: win / 2.0,
                upper_offset: win / 2.0,
            }),
            _ => None,
        };
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
    /// Streams every record sequentially from offset `HEADER_SIZE`.
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
        self.stream(path, |s| handler(s))
    }
}
```

- [ ] **Step 3c: 接入 `lib.rs`**

`crates/spectrum-io/src/lib.rs`：
- `detect_format` 的 `match` 在 `Some("mzml") => ...` 之后加：
```rust
        Some("pfb") => SpectrumFormat::Pfb,
```
- `create_reader` 的 `match info.format` 加臂：
```rust
        SpectrumFormat::Pfb => Box::new(pfb::PfbReader),
```
- `create_indexed_reader` 的 `match info.format` 加臂（**暂用** `PfbReader`，Task 2 切到 `IndexedPfbReader`）：
```rust
        SpectrumFormat::Pfb => Ok(Box::new(pfb::PfbReader)),
```

- [ ] **Step 4: 跑测试确认通过**

Run: `cargo test -p protein-copilot-spectrum-io pfb::`
Expected: PASS（10 个 pfb 测试通过）。
另跑回归：`cargo test -p protein-copilot-core spectrum`（Display 等仍通过）。

- [ ] **Step 5: 质量门 + 提交**

```bash
cargo fmt -p protein-copilot-spectrum-io -p protein-copilot-core
cargo clippy -p protein-copilot-spectrum-io --all-targets
git add crates/core/src/spectrum.rs crates/spectrum-io/src/pfb.rs crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add PfbReader + SpectrumFormat::Pfb for binary PFB

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: `indexed_pfb.rs` — `IndexedPfbReader`

**Files:**
- Create: `crates/spectrum-io/src/indexed_pfb.rs`
- Modify: `crates/spectrum-io/src/lib.rs`（`pub mod indexed_pfb;` + 导出 + 切换 `create_indexed_reader` 的 Pfb 臂）
- Test: 同文件 `#[cfg(test)] mod tests`

- [ ] **Step 1: 声明模块 + 写失败测试**

`crates/spectrum-io/src/lib.rs` 模块声明区加：
```rust
pub mod indexed_pfb;
```
并在导出区（`pub use indexed_mzml::IndexedMzMLReader;` 附近）加：
```rust
pub use indexed_pfb::IndexedPfbReader;
```

创建 `crates/spectrum-io/src/indexed_pfb.rs`，先放文档 + 测试（`IndexedPfbReader` 未定义→编译失败）：
```rust
//! Indexed PFB reader: O(1) scan lookup + cached scan metadata via a
//! reused [`crate::index::ScanIndex`].
//!
//! On [`IndexedPfbReader::open`], the footer offset table is read and each
//! record's property header is parsed (skipping peak blobs) to build a
//! `scan → (offset, rt, ms_level, isolation_window)` index. Bulk operations
//! delegate to the streaming [`crate::pfb::PfbReader`].

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pfb::PfbReader;
    use crate::reader::SpectrumReader;

    /// Writes a PFB file with the given records; returns (tempdir, path).
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
        assert!((metas[0].rt_min - 0.1).abs() < 1e-9); // 6s/60
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
        // scan 2: RT 12s=0.2min, isolation 500..502; precursor 501.25 inside.
        let hit = reader.find_by_rt(&p, 0.2, 501.25, 0.05).unwrap();
        assert_eq!(hit.map(|(scan, _)| scan), Some(2));
    }

    #[test]
    fn read_spectrum_not_found() {
        let (_d, p) = write_pfb(&sample_recs());
        let reader = IndexedPfbReader::open(&p).unwrap();
        let err = reader.read_spectrum(&p, 999).unwrap_err();
        assert!(matches!(err, SpectrumIoError::ScanNotFound { scan: 999, .. }));
    }
}
```

- [ ] **Step 2: 跑测试确认编译失败**

Run: `cargo test -p protein-copilot-spectrum-io indexed_pfb::`
Expected: 编译错误 `cannot find ... IndexedPfbReader`。

- [ ] **Step 3: 写 `indexed_pfb.rs` 实现（插入到文档之后、测试之前）**

```rust
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

/// Builds a [`ScanIndex`] by reading the footer offsets and each record's
/// property header (peak blobs are skipped — only seeks, no peak reads).
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
            None => continue, // skip malformed record in index
        };
        let ms_level = toks
            .get(1)
            .and_then(|t| t.trim().parse::<u8>().ok())
            .unwrap_or(0);
        let rt_seconds = toks
            .get(2)
            .and_then(|t| t.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let center = toks.get(7).and_then(|t| t.trim().parse::<f64>().ok());
        let window = toks.get(10).and_then(|t| t.trim().parse::<f64>().ok());
        let isolation_window = match (center, window) {
            (Some(c), Some(w)) if w > 0.0 => Some((c, w / 2.0, w / 2.0)),
            _ => None,
        };
        entries.insert(
            scan,
            ScanMeta {
                offset: off,
                rt_seconds,
                ms_level,
                isolation_window,
            },
        );
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
        Ok(self.index.find_by_rt(rt_min, precursor_mz, rt_tolerance_min))
    }
}
```

- [ ] **Step 4: 切换 `create_indexed_reader` 的 Pfb 臂**

`crates/spectrum-io/src/lib.rs` 的 `create_indexed_reader` 中，把 Task 1 临时的
```rust
        SpectrumFormat::Pfb => Ok(Box::new(pfb::PfbReader)),
```
改为：
```rust
        SpectrumFormat::Pfb => {
            let reader = IndexedPfbReader::open(path)?;
            Ok(Box::new(reader))
        }
```

- [ ] **Step 5: 跑测试确认通过**

Run: `cargo test -p protein-copilot-spectrum-io pfb:: indexed_pfb::`
Expected: PASS（Task 1 的 10 项 + Task 2 的 6 项均通过；含 Task 1 的 `detect_and_create_reader_for_pfb` 现在走 `IndexedPfbReader`）。

- [ ] **Step 6: 质量门 + 提交**

```bash
cargo fmt -p protein-copilot-spectrum-io
cargo clippy -p protein-copilot-spectrum-io --all-targets
git add crates/spectrum-io/src/indexed_pfb.rs crates/spectrum-io/src/lib.rs
git commit -m "feat(spectrum-io): add IndexedPfbReader (ScanIndex-backed PFB random access)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: PFB 端到端 XIC + 放宽遗留守卫

**Files:**
- Modify: `crates/xic/src/extract.rs`（放宽 2 处遗留守卫 + 端到端测试）
- Test: `crates/xic/src/extract.rs` 内 `#[cfg(test)] mod tests`

> **实现期发现（重要）**：真正被工具调用的是 `extract_xic_unified`（接收 `&dyn SpectrumReader`，**无格式守卫**）。`extract_xic`/`extract_xic_with_raw`（line ~282/~695 的 `!= MzML` 守卫所在）是**预先存在的死函数**（clippy 已标 never-used），不在实时路径上。因此 PFB→XIC 的实时支持由 Task 1/2（`create_indexed_reader → IndexedPfbReader`）达成；本任务的核心是**端到端测试**该实时路径，并顺带放宽遗留守卫以保持一致/防未来复活。

- [ ] **Step 1: 写失败的端到端测试**

在 `crates/xic/src/extract.rs` 的 `#[cfg(test)] mod tests { ... }` 内新增（若该 mod 不存在则创建在文件末尾）：
```rust
    #[test]
    fn pfb_end_to_end_xic_via_indexed_reader() {
        use protein_copilot_core::search_params::{MassTolerance, ToleranceUnit};

        // Synthetic .pfb: MS1 scan 1 + MS2 scans 2/3/4 (same isolation window).
        fn write_xic_pfb() -> (tempfile::TempDir, std::path::PathBuf) {
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("xic.pfb");
            let recs: Vec<(&str, Vec<f64>, Vec<f64>)> = vec![
                ("1\t1\t60.0\tFTMS", vec![400.0, 500.5, 600.0], vec![100.0, 5000.0, 200.0]),
                (
                    "2\t2\t60.5\tFTMS\t2\t1000.0\t50\t500.5\tHCD\t1\t2.0\t27.0\t500.5",
                    vec![150.0, 250.0],
                    vec![10.0, 20.0],
                ),
                (
                    "3\t2\t61.0\tFTMS\t2\t1000.0\t50\t500.5\tHCD\t1\t2.0\t27.0\t500.5",
                    vec![150.0, 250.0, 350.0],
                    vec![30.0, 40.0, 50.0],
                ),
                (
                    "4\t2\t61.5\tFTMS\t2\t1000.0\t50\t500.5\tHCD\t1\t2.0\t27.0\t500.5",
                    vec![150.0, 250.0],
                    vec![15.0, 25.0],
                ),
            ];
            let header_size: u64 = 24;
            let mut body: Vec<u8> = Vec::new();
            let mut offsets: Vec<u64> = Vec::new();
            for (prop, mz, inten) in &recs {
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

        let (_d, p) = write_xic_pfb();
        let reader = protein_copilot_spectrum_io::create_indexed_reader(&p).unwrap();
        let params = ExtractionParams {
            mz_tolerance: MassTolerance {
                value: 20.0,
                unit: ToleranceUnit::Ppm,
            },
            n_cycles: 2,
            top_n_ions: usize::MAX,
            label_type: None,
            intensity_rule: IntensityRule::MaxInWindow,
        };
        let result = extract_xic_unified(
            reader.as_ref(),
            &p,
            3,
            "PEPTIDEK",
            2,
            500.5,
            &[],
            &params,
            20.0,
        )
        .unwrap();
        assert_eq!(result.xic_data.target_scan, 3);
        assert!(
            !result.raw_scans.ms2_scans.is_empty(),
            "expected MS2 scans extracted from the PFB window"
        );
    }
```
> 若 `tests` 模块缺少 `use super::*;`，需补上（`extract_xic_unified`、`ExtractionParams`、`IntensityRule` 由此引入）。

- [ ] **Step 2: 跑测试确认失败/通过状态**

Run: `cargo test -p protein-copilot-xic pfb_end_to_end_xic_via_indexed_reader`
Expected: 该测试**已能通过**（因 Task 1/2 已让 `create_indexed_reader` 支持 PFB，且 `extract_xic_unified` 无格式守卫）。若失败，按报错修正后再继续。
> 这是回归/集成保护测试：它锁定「PFB 经索引读取器走通 XIC 实时路径」这一行为。

- [ ] **Step 3: 放宽两处遗留守卫（一致性 / 防未来复活）**

`crates/xic/src/extract.rs` 中有**两处文本完全相同**的守卫（在 `extract_xic` ~L282 与 `extract_xic_with_raw` ~L695）：
```rust
    if info.format != protein_copilot_core::spectrum::SpectrumFormat::MzML {
```
将**两处**都改为：
```rust
    if !matches!(
        info.format,
        protein_copilot_core::spectrum::SpectrumFormat::MzML
            | protein_copilot_core::spectrum::SpectrumFormat::Pfb
    ) {
```

- [ ] **Step 4: 跑测试 + 质量门**

Run:
```bash
cargo test -p protein-copilot-xic
cargo fmt -p protein-copilot-xic
cargo clippy -p protein-copilot-xic --all-targets
```
Expected: 测试全过；新端到端测试通过；clippy 无源于本次改动的新告警（`extract_xic`/`extract_xic_with_raw` 的 never-used 告警为**预先存在**，与本次无关）。

- [ ] **Step 5: 提交**

```bash
git add crates/xic/src/extract.rs
git commit -m "feat(xic): support PFB in XIC (end-to-end test) + relax legacy format guards

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## 完成标准

- [ ] `cargo test -p protein-copilot-spectrum-io -p protein-copilot-core -p protein-copilot-xic` 全绿。
- [ ] 各所改 crate `cargo fmt --check` 通过、`cargo clippy --all-targets` 无新增告警。
- [ ] `.pfb` 经 `detect_format`→`SpectrumFormat::Pfb`，`create_reader`→`PfbReader`，`create_indexed_reader`→`IndexedPfbReader`；读谱/summary/按 scan 随机读/scan 元数据/RT 查找均可用。
- [ ] PFB 可走通 XIC 实时路径（`extract_xic_unified` 端到端测试通过）。
- [ ] 实现期对真实样本 `.pfb`（`~/share/.../2th/*.pfb`，800MB，不入库）做一次冒烟：`create_indexed_reader` → 断言 scan_num=80096、scan 1=MS1、scan 2=MS2 且 precursor/窗口合理、末张 RT≈120min。
