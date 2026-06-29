# L4 — spectrum-io（crate 源码骨架）

> 上承 [L3-spectrum-io](L3-spectrum-io.md) 的子系统数据流与 [L2](L2-architecture.md) 的分层依赖。L4 不复述数据流，只聚焦本 crate **对外 API 表面与代码骨架**：导出了什么、签名长什么样、关键常量取值、谁在调。所有签名 / 路径 / 常量均核对自 `crates/spectrum-io/src`。

## 1. 用途 / 位置 / 依赖

一句话：把 `.mzML` / `.mgf` / `.pfb` 解析成 `core::spectrum` 类型，对外只暴露三个自由函数 + 一个 `SpectrumReader` trait。

位置 `crates/spectrum-io`，包名 `protein-copilot-spectrum-io`。生产依赖（非 dev）：

```
spectrum-io -> core   : Spectrum / SpectrumSummary / IsolationWindow / PrecursorInfo / MsLevel / SpectrumFormat / SpectrumFileInfo
spectrum-io -> 外部    : quick-xml(XML 流) + base64(解码) + flate2(zlib) + memchr(SIMD 扫描) + thiserror + tracing
```

dev-only：serde_json + tempfile。本 crate 不调 LLM、不做打分 / FDR。六个 reader 结构体虽全部 `pub`，但规范入口是三个工厂函数：测试 / 小文件用 `create_reader`，需多次随机访问用 `create_indexed_reader`，对外统一收敛到 `Box<dyn SpectrumReader>`。

## 2. 对外 API 表面

`lib.rs` 顶层导出（`util` 是唯一私有模块）：

```rust
pub mod disk_cache; pub mod error; pub mod index;          // 逐行 pub mod，此处简写
pub mod indexed_mgf; pub mod indexed_mzml; pub mod indexed_pfb;
pub mod mgf; pub mod mzml; pub mod pfb; pub mod reader;
mod util;                                        // 仅 crate 内可见
pub use error::SpectrumIoError;
pub use indexed_mgf::IndexedMgfReader;
pub use indexed_mzml::IndexedMzMLReader;
pub use indexed_pfb::IndexedPfbReader;
pub use reader::{ScanMetaInfo, SpectrumReader};  // 注意：Ms2ScanMeta 未顶层 re-export，走 reader::Ms2ScanMeta
```

三个入口函数：

```rust
pub fn detect_format(path: &Path) -> Result<SpectrumFileInfo, SpectrumIoError>;
pub fn create_reader(info: &SpectrumFileInfo) -> Box<dyn SpectrumReader>;
pub fn create_indexed_reader(path: &Path) -> Result<Box<dyn SpectrumReader>, SpectrumIoError>;
```

`SpectrumReader: Send + Sync`，4 个必实现 + 3 个带默认实现（索引读取器覆盖为零 I/O）：

```rust
pub trait SpectrumReader: Send + Sync {
    fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError>;
    fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError>;
    fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError>;
    fn for_each_spectrum(&self, path: &Path,
        handler: &mut dyn FnMut(Spectrum) -> Result<bool, SpectrumIoError>)
        -> Result<u32, SpectrumIoError>;
    // 默认实现走慢路径（read_all / for_each）；IndexedMzMLReader / IndexedPfbReader 覆盖：
    fn list_ms2_meta(&self, path: &Path)  -> Result<Vec<Ms2ScanMeta>, SpectrumIoError> { .. }
    fn list_scan_meta(&self, path: &Path) -> Result<Vec<ScanMetaInfo>, SpectrumIoError> { .. }
    fn find_by_rt(&self, path: &Path, rt_min: f64, precursor_mz: f64,
        rt_tolerance_min: f64) -> Result<Option<(u32, f64)>, SpectrumIoError> { .. }
}
```

两个元数据出参结构（`reader.rs`，`rt_min` 单位为分钟）：

```rust
pub struct ScanMetaInfo { pub scan_number: u32, pub ms_level: u8, pub rt_min: f64,
                          pub isolation_window: Option<(f64, f64, f64)> }
pub struct Ms2ScanMeta  { pub scan_number: u32, pub rt_min: f64,
                          pub isolation_window: Option<(f64, f64, f64)> }
```

读取契约（四个必实现方法共享）：所有 `read_summary` 都借 `util::SummaryAccumulator` 流式累计——逐谱 `observe()`、收尾 `into_summary()` 求峰数与隔离窗宽度中位数并 `validate()`，全程不把整文件载入内存（这正是 AI 推参只吃摘要的前提）。`for_each_spectrum` 的 handler 返回 `Ok(false)` 即提前停止，基础版 `read_spectrum` 正借此命中即停。每条谱最终都过 `core::Spectrum::new()` 的 `validate()`（等长、有限、m/z 升序且 >0、强度 >=0、scan >= 1），不满足即包成 `ValidationError { scan, detail }`。

## 3. 各 Reader 实现与关键常量

每种格式一对：流式基础版（单元结构体，无状态）+ 索引版（持 `ScanIndex` + `PathBuf`）。

| 结构体 | 文件 | 一行职责 |
|--------|------|----------|
| `MgfReader` | mgf.rs | BEGIN/END IONS 文本块流式解析 |
| `MzMLReader` | mzml.rs | quick-xml 事件流 + base64/zlib 二进制解码 |
| `PfbReader` | pfb.rs | pParse2+ 小端二进制（24B 头 + 记录 + footer 偏移表） |
| `IndexedMzMLReader { index, path }` | indexed_mzml.rs | `.idx` 缓存或字节扫描建 ScanIndex，seek O(1) |
| `IndexedMgfReader { index, path }` | indexed_mgf.rs | 扫 `BEGIN IONS` 字节位置建索引 |
| `IndexedPfbReader { index, path }` | indexed_pfb.rs | 读 footer 偏移表 + 属性头建索引 |

索引版的 `read_all` / `read_summary` / `for_each_spectrum` 一律委托回对应基础版；只有 `read_spectrum` 与三个 meta 方法走索引快路。三个 indexed 读取器来源不同（mzML 字节扫描、MGF 扫 `BEGIN IONS`、PFB 读 footer 偏移表），但最终都收敛到同一个 `ScanIndex`，因此 `find_by_rt` / `list_scan_meta` 等上层逻辑对三种格式完全一致。

关键常量（防御畸形输入 + 缓存格式）：

| 文件 | 常量 | 值 |
|------|------|----|
| mzml.rs | `MAX_PEAKS_PER_SPECTRUM` | `500_000` |
| pfb.rs | `MAX_PEAKS_PER_SCAN` | `10_000_000` |
| pfb.rs | `MAX_PROP_LEN` | `100_000_000` |
| disk_cache.rs | `MAGIC` | `b"PCIX"` |
| disk_cache.rs | `VERSION` | `2` |
| disk_cache.rs | `HEADER_SIZE` / `ENTRY_SIZE` | `25` / `46` 字节 |
| index.rs | `TAIL_READ_SIZE` / `MAX_INDEX_READ_SIZE` | `4096` / `10*1024*1024` |
| index.rs (byte_scan 内) | `CHUNK_SIZE` / `TAG_MIN_CONTENT` | `256*1024` / `8192` |

`index.rs` / `disk_cache.rs` 的公共面（L3 未逐一列）：

```rust
// index.rs
pub enum IndexSource { NativeIndex, BuiltFromScan }
pub struct ScanMeta { pub offset: u64, pub rt_seconds: f64, pub ms_level: u8,
                      pub isolation_window: Option<(f64, f64, f64)> }
pub struct ScanIndex { /* entries / source / rt_sorted，字段私有 */ }
//   方法：get_offset / get_meta / len / is_empty / source / iter_meta
//        / offsets / scan_numbers / rt_sorted / find_by_rt / new / from_meta
pub fn build_index_from_native_mzml(path: &Path) -> Result<Option<ScanIndex>, SpectrumIoError>;
pub fn build_index_by_scanning(path: &Path)      -> Result<ScanIndex, SpectrumIoError>;
pub fn build_index_by_byte_scan(path: &Path)     -> Result<ScanIndex, SpectrumIoError>;
// disk_cache.rs
pub fn idx_path(mzml_path: &Path) -> PathBuf;                              // 追加 ".idx"
pub fn file_metadata(path: &Path) -> Result<(u64, u64), SpectrumIoError>;  // (size, mtime_secs)
pub fn load_index(mzml_path: &Path, expected_size: u64, expected_mtime: u64)
    -> Result<Option<ScanIndex>, SpectrumIoError>;
pub fn save_index(mzml_path: &Path, index: &ScanIndex, file_size: u64, file_mtime: u64)
    -> Result<(), SpectrumIoError>;
```

`util.rs`（唯一私有模块）是三个基础读取器共用的骨架：

```rust
pub(crate) fn open_buffered(path: &Path) -> Result<BufReader<File>, SpectrumIoError>; // NotFound 与其它 I/O 分流
pub(crate) fn sort_peaks_by_mz(mz: &mut Vec<f64>, intensity: &mut Vec<f64>);          // 等长才排序
pub(crate) struct SummaryAccumulator { /* total / ms1_count / ms2_count / mz_min.. */ }
//   new() -> observe(&Spectrum) 逐谱累计 -> into_summary(path, format) -> SpectrumSummary
```

`sort_peaks_by_mz` 在 `mz.len() != intensity.len()` 时直接原样返回，把错配留给 `Spectrum::new()` 报干净的 `ArrayLengthMismatch`，避免越界 panic 或静默截断峰列。

错误类型集中在 `error.rs`：`SpectrumIoError`（thiserror，11 变体）+ `From<SpectrumIoError> for core::error::CoreError`（变体清单见 L3 第 6 节，此处不复述）。

## 4. 简化源码片段

`lib.rs` 扩展名分派（`detect_format` 大小写不敏感；`create_reader` 按 `format` 选实现）：

```rust
let format = match ext.as_deref() {       // ext = path.extension().to_lowercase()
    Some("mgf")  => SpectrumFormat::Mgf,
    Some("mzml") => SpectrumFormat::MzML,
    Some("pfb")  => SpectrumFormat::Pfb,
    _ => return Err(SpectrumIoError::UnknownFormat { path: path.to_path_buf() }),
};
// create_reader：单元结构体装箱
match info.format {
    SpectrumFormat::Mgf  => Box::new(mgf::MgfReader),
    SpectrumFormat::MzML => Box::new(mzml::MzMLReader),
    SpectrumFormat::Pfb  => Box::new(pfb::PfbReader),
}
// create_indexed_reader：先 detect_format，再 Indexed*::open(path) 装箱
```

`IndexedMzMLReader::open` 两层解析（命中缓存即返回，否则字节扫描后回写）：

```rust
// 层 1：磁盘缓存（PCIX v2），按 size + mtime 判新鲜
if let Ok((size, mtime)) = disk_cache::file_metadata(path) {
    if let Ok(Some(idx)) = disk_cache::load_index(path, size, mtime) {
        return Ok(Self { index: idx, path: path.to_path_buf() });
    }
}
// 层 2：memchr SIMD 字节扫描，一遍提取 offset + RT / ms_level / 隔离窗
let index = build_index_by_byte_scan(path)?;
let _ = disk_cache::save_index(path, &index, size, mtime); // 回写缓存，失败非致命
```

`mzml::decode_binary_array` 骨架（base64 -> 选 zlib -> 选 32/64bit -> 峰数封顶）：

```rust
pub(crate) fn decode_binary_array(b64: &str, meta: &BinaryArrayMeta, path: &Path)
    -> Result<Vec<f64>, SpectrumIoError> {
    let raw = base64::engine::general_purpose::STANDARD.decode(b64.trim())?;  // -> BinaryDecodeError
    let bytes = if meta.is_zlib {                                             // MS:1000574
        let limit = (MAX_PEAKS_PER_SPECTRUM * 8 + 1024) as u64;               // 防 zlib 炸弹
        let mut out = Vec::new();
        ZlibDecoder::new(&raw[..]).take(limit).read_to_end(&mut out)?;
        if out.len() as u64 == limit { return Err(/* 超限 BinaryDecodeError */); }
        out
    } else { raw };
    let values = if meta.is_64bit {                  // MS:1000523；bytes.len() % 8 != 0 即报错
        bytes.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()
    } else {                                         // MS:1000521；按 % 4 校验
        bytes.chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap()) as f64).collect()
    };
    if values.len() > MAX_PEAKS_PER_SPECTRUM { return Err(/* 超限 */); }
    Ok(values)
}
```

`IndexedMzMLReader::read_spectrum` 快路（offset -> seek -> 只解析单节点）：

```rust
let offset = self.index.get_offset(scan).ok_or(SpectrumIoError::ScanNotFound { .. })?;
self.read_spectrum_at_offset(scan, offset)   // BufReader.seek(SeekFrom::Start(offset)) + 解析单个 <spectrum>
```

## 5. 调用链（谁在调）

本 crate 被 search-engine / xic / result-import / entrapment-analysis / mcp-server 等生产 crate 依赖（param-recommend、integration-tests 仅 dev 引用）。

- **mcp-server**（MCP 服务 bin）：`get_or_create_reader(path)` 内调 `create_indexed_reader`，结果存入 `Arc<Mutex<lru::LruCache<PathBuf, Arc<dyn SpectrumReader>>>>`（容量 `NonZeroUsize::new(8)`），按 `canonicalize` 后路径复用；`read_summary` / `read_spectrum` / `find_by_rt` / `list_scan_meta` 均经此缓存读取器。
- **search-engine**：`simple_engine.rs` 与 `adapters/sage/mod.rs` 调 `create_indexed_reader` 后 `read_all` / `for_each_spectrum`（大文件流式，避免一次性载入）。
- **xic**：`extract.rs` 收注入的 `&dyn SpectrumReader`，先 `list_scan_meta` 规划目标 scan，再 `read_spectrum`（O(1) seek）逐个取谱。
- **entrapment-analysis**：`create_indexed_reader` -> `find_by_rt` 定位 MS2 -> `read_spectrum` 取谱做注释。
- **param-recommend**：仅 **dev-dependency**；生产代码消费的是 `core::SpectrumSummary`（本 crate `read_summary` 的产物，由 mcp-server 串接），不直接调 reader。
- **dia-extraction**：不依赖本 crate，经 `core::Spectrum` 流与 `SpectrumSummary.median_isolation_window_da` 间接受益。

跨边界出错时，`From<SpectrumIoError> for CoreError` 把本 crate 错误转成 `core::CoreError`，再由 mcp-server 包成结构化 MCP 错误（码 + 描述 + 建议），全程不 unwrap / expect。

## 6. 测试入口

源码内联 `#[cfg(test)] mod tests`（各 reader、index、disk_cache、util、lib），共 129 个单测；`tests/integration.rs` 14 个端到端（detect -> create -> read，含 corrupt / truncated / no_binary 容错）；lib.rs 文档示例 1 个 doc-test。fixtures 在 `tests/fixtures/`：`small.mgf` / `small.mzml`（各 10 谱）/ `small_indexed.mzml`（带 indexList）/ `small.mzml.idx`（PCIX 缓存）/ `no_binary.mzml` / `corrupt.{mgf,mzml}` / `truncated.mgf`。

```
cargo test -p protein-copilot-spectrum-io --offline
# 129 passed (unit) + 14 passed (integration) + 1 passed (doc-test)
```

只跑单个模块用 `cargo test -p protein-copilot-spectrum-io --offline mzml::`（按模块路径前缀过滤）。

—— 返回 [README](README.md)。
