# L3 — 谱图读取与解析子系统（spectrum-io）

承接 [L2](L2-architecture.md)。本篇深入 `crates/spectrum-io`：它如何把磁盘上的 mzML / MGF / PFB 文件解析成 `core` 的 `Spectrum` 与 `SpectrumSummary`，如何用索引做 O(1) 随机访问，以及如何防御畸形输入。所有签名、常量、CV accession 均核对自源码。

## 1. 职责与位置

一句话：**把异构的质谱文件统一解析为 `core::spectrum` 类型，并提供"全量 / 流式 / 按 scan 随机访问 / 元数据"四种读取姿势。**

```text
   磁盘文件                spectrum-io                          下游
.mzML/.mgf/.pfb  ->  detect_format + *Reader  -->  SpectrumSummary ----> param-recommend
                     (产出 core::spectrum 类型)    Vec<Spectrum>/流式 --> search-engine
                                                   list_scan_meta/find_by_rt -> xic / dia
```

上游是原始文件；spectrum-io 不调任何 LLM、不做打分/FDR。下游 `param-recommend` 吃 `SpectrumSummary` 推参，`search-engine` 吃 `Spectrum` 做酶切打分，`xic` / `dia-extraction` 吃元数据做定向提取。

## 2. 模块边界

每种格式有一对读取器：`mgf.rs` / `mzml.rs` / `pfb.rs` 是流式基础版，`indexed_*` 是带 `ScanIndex` 的随机访问版（批量操作仍委托回基础版）；`util.rs` 与 `error.rs` 横向支撑全部读取器。

| 文件 | 职责 |
|------|------|
| `lib.rs` | 入口：`detect_format` / `create_reader` / `create_indexed_reader` + 模块导出 |
| `reader.rs` | `SpectrumReader` trait（统一接口）+ `Ms2ScanMeta` / `ScanMetaInfo` |
| `mgf.rs` | `MgfReader`：BEGIN/END IONS 文本流式解析 |
| `mzml.rs` | `MzMLReader`：quick-xml 流式解析 + base64/zlib 二进制解码 |
| `indexed_mzml.rs` | `IndexedMzMLReader`：`ScanIndex` 支撑的 seek 随机访问 |
| `indexed_mgf.rs` | `IndexedMgfReader`：BEGIN IONS 字节偏移索引 |
| `indexed_pfb.rs` | `IndexedPfbReader`：footer 偏移表索引 |
| `index.rs` | `ScanIndex` / `ScanMeta` / `IndexSource` + 三种构建路径 + `find_by_rt` |
| `disk_cache.rs` | PCIX v2 `.idx` 磁盘缓存读写（头 25B + 每条 46B） |
| `pfb.rs` | `PfbReader`：pParse2+ 小端二进制解析 |
| `util.rs` | `open_buffered` / `sort_peaks_by_mz` / `SummaryAccumulator` |
| `error.rs` | `SpectrumIoError` + `From<SpectrumIoError> for CoreError` |

## 3. 关键数据结构（core::spectrum）

```rust
pub struct Spectrum {
    pub scan_number: u32,               // 1-based（0 触发 ZeroScanNumber）
    pub ms_level: MsLevel,              // MS1 / MS2 / Other(u8)
    pub retention_time_min: f64,        // 保留时间，内部统一"分钟"
    pub precursors: Vec<PrecursorInfo>, // MS1 空；DDA 通常 1；DIA 0 或 1
    pub mz_array: Vec<f64>,             // m/z，升序、有限、>0
    pub intensity_array: Vec<f64>,      // 强度 counts，有限、>=0，与 mz_array 等长
}

pub struct PrecursorInfo {
    pub mz: f64,
    pub charge: Option<i32>,            // DIA 常为 None
    pub intensity: Option<f64>,
    pub isolation_window: Option<IsolationWindow>,
    pub source_scan: Option<u32>,       // mzML spectrumRef 指向的 MS1 scan
}

pub struct IsolationWindow {            // 对齐 mzML <isolationWindow>
    pub target_mz: f64,                 // 窗口中心 m/z
    pub lower_offset: f64,              // 下沿偏移（m/z, >=0）
    pub upper_offset: f64,              // 上沿偏移（m/z, >=0）
}

pub struct SpectrumSummary {            // 推参的唯一输入
    pub file_path: String,
    pub format: SpectrumFormat,         // MzML / Mgf / Pfb
    pub total_spectra: u64,
    pub ms1_count: u64,
    pub ms2_count: u64,
    pub mz_range: [f64; 2],             // [min, max]
    pub rt_range_min: [f64; 2],         // [min, max]，分钟
    pub precursor_charge_distribution: HashMap<i32, u64>,
    pub median_peaks_per_spectrum: u32,
    pub median_isolation_window_da: Option<f64>, // 用于 DIA 判别
}
```

单位约定：**m/z 无量纲；保留时间内部统一分钟（`retention_time_min` / `rt_range_min`）；scan_number 从 1 起；强度为探测器 counts；隔离窗偏移以 m/z（约等于 Da）计且 ≥0。** 一处易错点：`index.rs` 的 `ScanMeta.rt_seconds` 存的是**秒**，`list_*_meta` 出口处统一 `/ 60.0` 换成分钟（`Ms2ScanMeta.rt_min` / `ScanMetaInfo.rt_min`）。

`Spectrum::new()` 一律走 `validate()`：等长、有限、mz 升序且 >0、强度 ≥0、scan ≥1、隔离窗偏移 ≥0——任一不满足返回 `SpectrumError`，被读取器包成 `ValidationError { scan, detail }`。

## 4. 主流程

### 4.1 格式探测与读取器选择

```text
detect_format(path) -> SpectrumFileInfo:
  if !path.exists() -> FileNotFound
  ext.to_lowercase():                 # 大小写不敏感
    "mgf"  -> Mgf
    "mzml" -> MzML
    "pfb"  -> Pfb
    _      -> UnknownFormat
  返回 SpectrumFileInfo { path, format, file_size_bytes }

create_reader(&info)        -> Box<dyn SpectrumReader>   # 无索引；测试/小文件
create_indexed_reader(path) -> Box<dyn SpectrumReader>   # Indexed*；多次 read_spectrum
```

`SpectrumReader: Send + Sync` 定义四个核心方法 + 三个可被索引读取器覆盖的默认方法：

```rust
fn read_all(&self, path: &Path) -> Result<Vec<Spectrum>, SpectrumIoError>;
fn read_summary(&self, path: &Path) -> Result<SpectrumSummary, SpectrumIoError>;
fn read_spectrum(&self, path: &Path, scan: u32) -> Result<Spectrum, SpectrumIoError>;
fn for_each_spectrum(&self, path, &mut dyn FnMut(Spectrum)->Result<bool,_>) -> Result<u32,_>;
// 默认实现走 read_all/for_each（慢）；IndexedMzMLReader 覆盖为零 I/O / O(logN)：
fn list_ms2_meta(..) ;  fn list_scan_meta(..) ;  fn find_by_rt(path, rt_min, mz, tol) ;
```

无论哪种读取器，`read_summary` 都用 `util::SummaryAccumulator` 流式累计（逐谱 `observe`，最后 `into_summary` 求中位数并 `validate`），不把整个文件一次性载入内存——这正是 AI 推参只需"摘要"而非全量谱图的前提。`for_each_spectrum` 的 handler 返回 `Ok(false)` 即提前停止，基础版 `read_spectrum` 正是借此"命中即停"。

### 4.2 mzML：流式 XML + base64/zlib 解码

`MzMLReader` 用 quick-xml 事件流，按 CV accession 填 `SpectrumBuilder`：`MS:1000511` ms level、`MS:1000016` 扫描起始时间、`MS:1000744/041/042` 前体 m/z/电荷/强度、`MS:1000827/828/829` 隔离窗 target/lower/upper。`</precursor>` 时 `flush_precursor()`，`</spectrum>` 时 `build()`。二进制数组解码（简化）：

```rust
fn decode_binary_array(b64_text, meta, path) -> Result<Vec<f64>, _> {
    let raw = base64::STANDARD.decode(b64_text.trim())?;        // base64
    let bytes = if meta.is_zlib {                               // MS:1000574
        let limit = (MAX_PEAKS_PER_SPECTRUM*8 + 1024) as u64;   // 防 zlib 炸弹
        ZlibDecoder::new(&raw[..]).take(limit).read_to_end(&mut out)?;
        if out.len() as u64 == limit { return Err(BinaryDecodeError) } // 触顶即报错
        out
    } else { raw };
    let values = if meta.is_64bit {                             // MS:1000523
        bytes.chunks_exact(8).map(f64::from_le_bytes)..         // len %8 != 0 -> 报错
    } else {                                                    // MS:1000521
        bytes.chunks_exact(4).map(|c| f32::from_le_bytes(c) as f64)..  // %4
    };
    if values.len() > MAX_PEAKS_PER_SPECTRUM { return Err(BinaryDecodeError) } // 500_000
    Ok(values)
}
```

保留时间单位换算（内部存分钟）：

```rust
rt_min = match unitAccession {           // MS:1000016 的 value
    "UO:0000031" => v,        // 已是分钟
    "UO:0000010" => v / 60.0, // 秒 -> 分钟
    ""           => { warn!("缺 unitAccession，按蛋白组学惯例当分钟"); v }
};
```

`<spectrum id="...scan=123">` 经 `parse_scan_from_id` 取 scan；`spectrumRef="controllerType=0 ... scan=1234"` 经 `parse_scan_from_spectrum_ref` 取 source_scan（取 `scan=` 后的连续数字，回退整串解析）。

### 4.3 MGF：文本块解析

```text
逐行（trim 后）：
  空行 / '#' 开头        -> 跳过
  "BEGIN IONS"          -> 新建块，fallback_scan += 1（缺 SCANS 时兜底）
  "PEPMASS=mz [int]"  "CHARGE=2+|3-"  "RTINSECONDS=s"(/60 -> 分钟)  "SCANS=n"  TITLE(跳过)
  "mz intensity"        -> 追加峰（split_whitespace >= 2 个可解析数值）
  "END IONS"            -> sort_peaks_by_mz -> Spectrum::new(MS2, isolation_window=None)
缺 END IONS 的截断块：收尾时仍尝试建谱（容错）；CHARGE=0 -> None（无物理意义）
```

### 4.4 随机访问与索引

`IndexedMzMLReader::open` 两层解析：(1) 读 `.mzML.idx`（PCIX v2 磁盘缓存，含 RT/ms_level/隔离窗，按 size+mtime 判新鲜）；(2) `build_index_by_byte_scan`（memchr SIMD 扫 `<spectrum ` needle，256KB 分块 + 8192B 元数据窗口、边界处整块延后），成功后回写 `.idx`。`read_spectrum` 用 `ScanIndex::get_offset(scan)` -> `seek(SeekFrom::Start)` -> 只解析该节点（O(1)）；`read_all/summary/for_each` 委托回普通 `MzMLReader`。

```rust
pub struct ScanMeta { offset: u64, rt_seconds: f64, ms_level: u8,
                      isolation_window: Option<(f64,f64,f64)> }   // index.rs
// find_by_rt：在预排序 rt_sorted: Vec<(f64,u32)> 上 partition_point 二分，
//             过滤 ms_level==2 与隔离窗命中，取 |Δrt| 最小者 -> Option<(scan, Δmin)>
```

PFB / MGF 的 indexed 读取器同理：`IndexedPfbReader` 读 footer 的 `i64×scan_num` 偏移表，`IndexedMgfReader` 扫 `BEGIN IONS` 字节位置，二者都建成同一 `ScanIndex`。磁盘缓存 `.idx` 为 PCIX v2 小端格式：25B 头（`PCIX` 魔数 + 版本 + 源文件 size/mtime + 条数）后接每条 46B（scan、offset、rt、ms_level、隔离窗），按 size+mtime 比对判定是否过期。

### 4.5 DDA vs DIA 判别

spectrum-io 在 `read_summary` 里把 MS2 隔离窗宽度（`lower_offset + upper_offset`）的中位数写进 `median_isolation_window_da`；真正分类在 `dia-extraction::detect_acquisition_mode`：

```rust
let widths = MS2 各前体 (lower_offset + upper_offset);
if widths.is_empty() { Unknown }            // 无 MS2 / 无隔离窗
else if median(widths) > threshold_da { DIA }   // 默认 5 Da：DDA 1-3，DIA 10-25
else { DDA }
```

`param-recommend::detect_dia` 则直接看 `median_isolation_window_da > 5.0`，与上述阈值一致。

## 5. 跨 crate 交互

- **依赖（非 dev）**：仅 `core`（`Spectrum / SpectrumSummary / IsolationWindow / PrecursorInfo / MsLevel / SpectrumFormat / SpectrumFileInfo`）+ quick-xml / base64 / flate2 / memchr。
- **`param-recommend`**：消费 `SpectrumSummary`——`detect_dia` 看 `median_isolation_window_da`，`infer_instrument` 按 `mz_range[1]` 与 `median_peaks_per_spectrum` 打分（高分辨 vs 低分辨）——产出 `AiDecision<SearchParams>`。
- **`search-engine`**：`read_all` / `for_each_spectrum` 喂酶切打分；大 DIA 文件用流式避免一次性载入内存。
- **`xic` / `dia-extraction`**：`list_scan_meta` / `find_by_rt` 先定位目标 scan，再只读需要的谱图。
- **`mcp-server`**：`get_or_create_reader` 用 `lru::LruCache`（容量 8）缓存 `create_indexed_reader` 的结果，按 canonical 路径复用 `IndexedMzMLReader`；错误经 `From<SpectrumIoError> for CoreError` 转成结构化 MCP 错误（码 + 描述 + 建议）。

## 6. 错误处理与防御

`SpectrumIoError`（thiserror）变体：`FileNotFound` / `UnknownFormat` / `UnsupportedFormat` / `IoError{path,source}` / `ParseError{path,line,detail}` / `XmlError` / `BinaryDecodeError` / `ValidationError{scan,detail}` / `ScanNotFound{path,scan}` / `IndexParseError` / `DiskCacheError`。**无 blanket `From<io::Error>`**——所有 I/O 错误在 callsite 用 `map_err` 补上 `path` 上下文。

对畸形输入的防御：

- **峰数上限**：mzML 单谱 `MAX_PEAKS_PER_SPECTRUM = 500_000`；PFB 单 scan `MAX_PEAKS_PER_SCAN = 10_000_000`、属性串 `MAX_PROP_LEN = 100_000_000`，超限报错而非 OOM。
- **zlib 炸弹**：解压用 `.take(MAX_PEAKS*8 + 1024)` 封顶，触顶即 `BinaryDecodeError`。
- **字节长度**：64-bit 数组按 `% 8`、32-bit 按 `% 4` 校验，不整除即报错。
- **数组错配**：`sort_peaks_by_mz` 发现 `mz.len() != intensity.len()` 直接原样返回，留给 `Spectrum::new()` 报干净的 `ArrayLengthMismatch`（避免越界 panic 或静默截断）。
- **缺字段容错**：MGF 缺 `SCANS` 用 `fallback_scan`、缺 `END IONS` 仍建谱、`CHARGE=0` 归 None；mzML 缺 RT 默认 0.0 并 `debug!`。
- **磁盘缓存自愈**：`.idx` 缺失 / 版本不符 / size+mtime 失配 -> `load_index` 返回 `Ok(None)`（非致命），回退字节扫描并重建。

—— 返回 [README](README.md)；逐 crate 源码级细节见 L4。

