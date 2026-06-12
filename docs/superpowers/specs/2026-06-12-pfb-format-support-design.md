# PFB 二进制谱图格式支持设计文档

> **日期**: 2026-06-12
> **状态**: 设计确认（格式布局已用真实样本核实；范围经终端确认）
> **方案**: 在 `spectrum-io` 新增 PFB 读取器（`PfbReader` + `IndexedPfbReader`），完整对齐 mzML/mgf 的 `SpectrumReader` 接入方式

---

## 1. 需求概述

新增对 **pXtract3 / pParse2+ 二进制 PFB** 谱图格式的输入支持，使 PFB 文件能像 mzML/mgf 一样被读取，用于读谱、参数推荐、搜索、标注与 XIC。

PFB 自带 footer 偏移索引，天然支持按 scan 随机读取。真实文件较大（样本 ~750–825 MB / 80,096 张谱图），因此提供带索引的读取器以保证 XIC 等按 scan 随机读场景的性能。

### 1.1 格式核实（基于真实样本）

样本：`.../2th/20190830_HF_ZHW_hela_SILAC_DDIA_500_550_2Da_Rep1.pfb`（813,753,733 字节，scan_num=80,096）。已用 Python 解析头部、首尾记录与 footer，确认：

- **字节序：小端（little-endian）**。
- 头部 24 字节：`3×i32`(全 0，预留) + `i64` addr_list_addr + `i32` scan_num。`file_size - addr_list_addr == scan_num*8`（footer 恰好是 scan_num 个 `i64`）。
- footer 偏移单调递增，`addrs[0]==24`（紧接头部），即「记录序→文件偏移」表。
- **RT 单位：秒**。末张 scan 80096 的 RT=7200.175（=120 min 标准梯度），中点 ≈3536 s。映射时 **/60 转分钟**。
- MS1 记录 property_str 4 字段；MS2 记录 13 字段（含 pXtract3 新增的 ActivationWindow / NCE / monoisotopicMz）。
- 谱峰强度为 **8 字节 double**（与表格/读取示例一致；写入示例里的 `float` 为笔误，已与用户确认）。

---

## 2. 设计决策汇总

| 决策项 | 选择 | 理由 |
|--------|------|------|
| 读取器归属 | `spectrum-io` crate 新增 `pfb.rs` / `indexed_pfb.rs` | PFB 是谱图格式，与 mgf/mzml 同层 |
| 范围 | 完整对齐：`PfbReader`（流式）+ `IndexedPfbReader`（footer/属性头索引） | 800MB 文件的 XIC 按 scan 随机读需高效 |
| 字节序 | 小端固定 | 样本核实；pXtract/pParse 在 Windows x86 生成 |
| 强度宽度 | 8 字节 `f64` | 样本核实 + 用户确认 |
| RT 单位 | 秒 → `/60` 存为 `retention_time_min` | 样本核实（末张 7200s=120min） |
| 扩展名 | `.pfb`（v1 仅此） | 用户当前需求；`.pfc` 伴随文件不在 v1 |
| 写出 | 不做（只读） | 输入格式支持 |
| 索引持久化 | 内存索引（打开时构建），不落盘 | YAGNI；mzML 的磁盘缓存是后续优化 |
| precursor m/z | monoisotopicMz（缺则 ActivationCenter） | DIA 下即窗口中心；与现有 XIC 用法一致 |
| 解析容错 | 可选字段缺失/解析失败→`None`，不失败整张谱 | WIFF/旧版可能缺 ActivationWindow/NCE 等 |
| XIC 守卫 | `xic/extract.rs` 放宽为接受 `MzML | Pfb` | PFB 同样有 MS1+MS2+隔离窗口 |
| 错误类型 | 复用 `SpectrumIoError` | 二进制解析→`ParseError`，按 scan→`ScanNotFound` |

---

## 3. 二进制格式规范

```text
Header (24 bytes, little-endian)
  i32  empty_1            (忽略)
  i32  empty_2            (忽略)
  i32  empty_3            (忽略)
  i64  addr_list_addr     footer 首地址
  i32  scan_num           谱图总数 (MS1+MS2)

Loop Body  ×scan_num   (顺序排列；偏移由 footer 记录)
  i32      property_str_len
  u8[len]  property_str   UTF-8, '\t' 分隔, 去尾 '\0'
  i32      peak_num
  f64[peak_num]  all_peak_mz          (升序)
  f64[peak_num]  all_peak_intensity

Footer @ addr_list_addr
  i64[scan_num]  scan_idx_addr        每条记录的起始文件偏移 (记录序)
```

property_str 字段（`\t` 分隔，按位置）：

| idx | 字段 | 类型 | MS 级别 |
|----:|------|------|---------|
| 0 | Scan | int | 通用 |
| 1 | MsType | int (1=MS1, 2=MS2) | 通用 |
| 2 | RT | float（**秒**） | 通用 |
| 3 | InstrumentType | string | 通用 |
| 4 | Charge | int | MS2 |
| 5 | MH+ | float | MS2 |
| 6 | IonInjectionTime | float | MS2 |
| 7 | ActivationCenter | float | MS2 |
| 8 | ActivationType | string | MS2 |
| 9 | PrecursorScan | int | MS2 |
| 10 | ActivationWindow | double | MS2（pXtract3 新增；可能缺失） |
| 11 | CollisionEnergy(NCE) | double | MS2（pXtract3 新增；可能缺失） |
| 12 | monoisotopicMz | double | MS2 |

---

## 4. property_str → `Spectrum` 映射

```text
scan_number          = parse_u32(tok[0])              // 必需
ms_level             = MsType==1 → MS1, ==2 → MS2, 其它 → Other(n)
retention_time_min   = parse_f64(tok[2]) / 60.0       // 秒→分钟；解析失败→0.0
mz_array/intensity_array = 解析得到的峰（按 m/z 升序）

MS1: precursors = []
MS2: precursors = [ PrecursorInfo {
        mz:             monoisotopicMz(tok[12]) 若缺/解析失败 → ActivationCenter(tok[7])
        charge:         Some(Charge) 若 >0；否则 None
        intensity:      None
        isolation_window: 若 ActivationCenter 与 ActivationWindow 均可解析 →
                          Some{ target_mz: ActivationCenter, lower_offset: ActivationWindow/2,
                                upper_offset: ActivationWindow/2 }；否则 None
        source_scan:    Some(PrecursorScan) 若可解析；否则 None
     } ]
```

- 经 `Spectrum::new(...)` 构造并校验；校验失败 → `SpectrumIoError::ValidationError { scan, detail }`。
- 容错：`tok` 长度不足或某可选字段解析失败时，对应字段取 `None`，不让整张谱失败；但 `Scan`/`MsType` 不可解析时该记录视为损坏 → `ParseError`。

---

## 5. 解析与读取器

### 5.1 解析原语（`pfb.rs` 私有）

- `read_header(&mut impl Read) -> (addr_list_addr: i64, scan_num: u32)`。
- `read_record(&mut impl Read) -> PfbRecord { property_str, mz, intensity }`（顺序读：prop_len → prop_str → peak_num → mz → intensity）。
- `record_to_spectrum(record) -> Result<Spectrum, SpectrumIoError>`（§4 映射）。
- 所有多字节整数/浮点用小端读取（`i32::from_le_bytes` 等）；短读/截断 → `ParseError`。

### 5.2 `PfbReader`（无状态，`create_reader` 用）

实现 `SpectrumReader`：
- `read_all` / `read_summary` / `for_each_spectrum`：从偏移 24 起顺序读 scan_num 条记录，逐张构造 `Spectrum`。`read_summary` 复用 `util::SummaryAccumulator`。
- `read_spectrum(scan)`：读头部 + footer 偏移表；按 footer 偏移逐条只读「属性头」（prop_len+prop_str，再 seek 跳过峰区 `peak_num*16`），匹配 `Scan==scan` 后回读该记录的峰，构造 `Spectrum`；未找到 → `ScanNotFound`。

### 5.3 `IndexedPfbReader`（`indexed_pfb.rs`，`create_indexed_reader` 用）

- `open(path)`：一次流式扫描——读头部，再按记录顺序读每条「属性头」（解析 Scan/MsType/RT/ActivationCenter+Window），用 seek 跳过峰区，构建并缓存：
  - `scan_to_offset: HashMap<u32, u64>`
  - `scan_meta: Vec<ScanMetaInfo>`（scan, ms_level, rt_min, isolation_window）
- `read_spectrum(scan)`：从 `scan_to_offset` 取偏移，seek 后读整条记录 → O(1)。
- `list_scan_meta` / `list_ms2_meta`：直接返回缓存。
- `find_by_rt`：在缓存上筛 MS2 + RT 容差 + 隔离窗口包含，取 RT 最近者。
- `read_all` / `read_summary` / `for_each_spectrum`：委托 `PfbReader`（流式）。

> 索引构建只读属性头（~每条数十字节）并 seek 跳过峰 blob，对 80k 谱图为一次性、亚秒级开销，避免读入 ~800MB 峰数据。

---

## 6. 接入与枚举变更

- `core/src/spectrum.rs`：`SpectrumFormat` 增 `Pfb`；`Display` 增 `Pfb => "pfb"`。
- `spectrum-io/src/lib.rs`：
  - `detect_format`：`Some("pfb") => SpectrumFormat::Pfb`。
  - `create_reader`：`SpectrumFormat::Pfb => Box::new(pfb::PfbReader)`（穷尽匹配，必须加臂）。
  - `create_indexed_reader`：`SpectrumFormat::Pfb => Box::new(IndexedPfbReader::open(path)?)`（穷尽匹配，必须加臂）。
  - 新增 `pub mod pfb; pub mod indexed_pfb;` 与 `pub use indexed_pfb::IndexedPfbReader;`。
- `xic/src/extract.rs:282,695`：守卫由 `info.format != MzML` 改为 `!matches!(info.format, MzML | Pfb)`，使 PFB 可做 XIC（PFB 提供 MS1+MS2+隔离窗口）。
- （可选）`error.rs` 中 `UnknownFormat`/`UnsupportedFormat` 的 `supported` 列表追加 `"pfb"`。

> 其余 `SpectrumFormat` 出现点（`param-recommend/rules.rs`、各测试、`into_summary` 调用）为构造或非穷尽匹配，不受影响；仅 `Display`、`create_reader`、`create_indexed_reader` 三处穷尽点必须更新。

---

## 7. 测试计划

### 7.1 单元测试（合成小 `.pfb`，测试内写字节）

测试辅助 `write_min_pfb(scans) -> tempfile`：按 §3 布局写 1 张 MS1 + 2 张 MS2（含 footer），再断言：

- 头部解析：scan_num、addr_list_addr 正确。
- `read_all`：3 张谱；MS1 无 precursor；MS2 的 precursor mz=monoisotopicMz、charge、isolation_window(target=ActivationCenter, ±ActivationWindow/2)、source_scan 正确。
- RT 换算：秒 `/60` → `retention_time_min`。
- 峰：mz/intensity 数量、数值、升序。
- footer 索引：`IndexedPfbReader` 与 `PfbReader` 的 `read_spectrum(scan)` 返回一致；O(1) 取到目标 scan。
- `read_summary`：total/ms1/ms2 计数、format=Pfb。
- 错误：`read_spectrum(不存在)` → `ScanNotFound`；截断文件 → `ParseError`。
- `detect_format(.pfb)` → Pfb；`create_reader`/`create_indexed_reader` 返回 PFB 读取器；`SpectrumFormat::Pfb.to_string()=="pfb"`。
- 容错：MS2 仅 10 字段（缺 ActivationWindow/NCE/monoisotopicMz）时，mz 回退 ActivationCenter、isolation_window=None，不 panic。

### 7.2 真实样本冒烟（实现期，不入库）

实现期对真实样本 `.pfb` 执行：`detect_format` → `create_indexed_reader` → 断言 scan_num=80096、scan 1=MS1、scan 2=MS2 且 precursor/窗口合理、RT 换算（末张≈120min）。样本 800MB 不提交仓库。

---

## 8. 范围边界（v1）

**做**：读取 `.pfb`（流式 + 索引），接入 detect/create/indexed + `SpectrumFormat::Pfb`，放宽 XIC 守卫支持 PFB。

**不做**：写出 PFB；`.pfc` 伴随文件；索引落盘持久化；非小端。

---

## 9. 涉及文件

| 文件 | 变更 |
|------|------|
| `crates/core/src/spectrum.rs` | `SpectrumFormat::Pfb` + `Display` |
| `crates/spectrum-io/src/pfb.rs` | 新增：解析原语 + `PfbReader` |
| `crates/spectrum-io/src/indexed_pfb.rs` | 新增：`IndexedPfbReader` |
| `crates/spectrum-io/src/lib.rs` | 模块声明/导出 + detect_format + create_reader + create_indexed_reader |
| `crates/spectrum-io/src/error.rs` | （可选）supported 列表追加 "pfb" |
| `crates/xic/src/extract.rs` | 放宽 2 处格式守卫为 `MzML | Pfb` |
| `crates/spectrum-io/tests/` 或内联 | PFB 单元测试 + 合成 fixture 辅助 |
