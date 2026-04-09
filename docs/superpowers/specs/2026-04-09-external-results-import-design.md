# External Search Results Import — Design Specification

> **Date**: 2026-04-09
> **Status**: Draft
> **Scope**: New `result-import` lib crate + MCP tool `import_search_results`

## 1. Problem Statement

当前 ProteinCopilot 的谱图标注（`annotate_spectrum`）和 XIC 提取（`extract_xic`）仅支持通过内部搜索产生的 `run_id` 访问 PSM 结果。用户需要导入外部搜索引擎的结果文件（DIA-NN parquet、自定义 JSON、pFind .spectra）来使用同样的标注和可视化能力。

**核心挑战**：外部结果文件（hela.json、DIA-NN report.parquet）不包含 scan number，需要通过 RT + DIA isolation window 匹配到 mzML 中的 MS2 谱图。

## 2. Design Overview

三层架构：**格式解析 → Scan 匹配 → 集成到现有流程**

```
外部结果文件 + mzML 目录
       │
       ▼
  ┌─────────────────────────────┐
  │  Step 1: 格式解析           │
  │  JSON / Parquet / TSV       │
  │  → Vec<ImportedPsm>         │
  │  (Unimod 查表转修饰)        │
  │  (RT 统一为秒)              │
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Step 2: Scan 匹配          │
  │  RT proximity +             │
  │  DIA isolation window       │
  │  → 填充 matched_scan        │
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Step 3: 转换 + 缓存        │
  │  ImportedPsm → core::Psm   │
  │  → SearchResult             │
  │  → run_results_cache        │
  │  返回 run_id                │
  └─────────────┬───────────────┘
                │
                ▼
  annotate_spectrum(run_id, scan)  ✅
  extract_xic(run_id, scan)       ✅
  generate_summary(run_id)        ✅
```

## 3. RT 单位约定

项目内部统一使用**秒（seconds）**作为 RT 单位：

- `core::Spectrum.retention_time_sec: f64` — 谱图 RT，单位秒
- `core::SpectrumSummary.rt_range_sec` — RT 范围，单位秒
- `spectrum-io` mzML 解析器遇到分钟单位（`UO:0000031`）自动 `× 60.0` 转秒

**外部数据转换规则**（在解析器内完成）：

| 格式 | 原始 RT 单位 | 转换 |
|------|-------------|------|
| custom_json (hela.json) | 分钟 | `rt * 60.0` |
| DIA-NN parquet | 分钟 | `RT * 60.0` |
| pFind .spectra | 待确认 | 预留转换接口 |

所有 `ImportedPsm` 的 `rt_sec` 字段在解析完成后已经是秒。后续 scan 匹配、MatchReport 等均使用秒。

## 4. Core Data Structures

### 4.1 ImportedPsm

```rust
/// 从外部文件解析出的 PSM（scan 匹配前）
/// 位于 crates/result-import/src/lib.rs
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImportedPsm {
    /// 肽段氨基酸序列
    pub sequence: String,
    /// 电荷态
    pub charge: i32,
    /// 母离子 m/z
    pub precursor_mz: f64,
    /// 保留时间（秒），解析时已从外部格式转换
    pub rt_sec: f64,
    /// 修饰列表（已通过 Unimod 查表转换为 core::Modification）
    pub modifications: Vec<Modification>,
    /// 打分（DIA-NN 有 Q.Value，custom_json 可能没有）
    pub score: Option<f64>,
    /// 蛋白质 accession 列表
    pub protein_accessions: Vec<String>,
    /// 关联的 raw 文件名（不含扩展名）
    pub raw_name: String,
    /// Scan 匹配后填充
    pub matched_scan: Option<u32>,
    /// 匹配 RT 偏差（秒）
    pub rt_delta_sec: Option<f64>,
}
```

### 4.2 ImportResult

```rust
/// import_search_results 的完整返回
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImportResult {
    /// 生成的 run_id，可用于 annotate_spectrum / extract_xic
    pub run_id: String,
    /// 匹配统计报告
    pub match_report: MatchReport,
    /// 导入的 PSM 数量（匹配成功的）
    pub imported_psm_count: usize,
    /// 唯一肽段数
    pub unique_peptides: usize,
    /// 蛋白质数
    pub protein_count: usize,
    /// 涉及的 raw 文件列表
    pub raw_files: Vec<String>,
}
```

### 4.3 MatchReport

```rust
/// Scan 匹配质量报告
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchReport {
    pub total_psms: usize,
    pub matched: usize,
    pub unmatched: usize,
    pub median_rt_delta_sec: f64,
    pub max_rt_delta_sec: f64,
    /// 每个 raw 文件的匹配统计
    pub per_file: HashMap<String, FileMatchStats>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileMatchStats {
    pub total: usize,
    pub matched: usize,
    pub ms2_count: usize,
}
```

### 4.4 UnimodDb

```rust
/// Unimod 修饰数据库
/// 位于 crates/result-import/src/unimod.rs
pub struct UnimodDb {
    entries: HashMap<u32, UnimodEntry>,
}

pub struct UnimodEntry {
    pub record_id: u32,
    pub title: String,       // e.g. "Oxidation"
    pub mono_mass: f64,      // e.g. 15.994915
    pub residues: Vec<char>, // e.g. ['M', 'W', 'H']
}

impl UnimodDb {
    /// 从 unimod.xml 加载完整数据库
    pub fn from_xml(path: &Path) -> Result<Self>;

    /// 内置精简版（~20 个常见修饰）
    pub fn builtin() -> Self;

    /// 查表：ID → Modification
    pub fn to_modification(&self, record_id: u32, position: usize, sequence: &str)
        -> Result<Modification>;
}
```

## 5. Format Parsers

### 5.1 Parser Trait

```rust
/// 所有格式解析器实现此 trait
/// 位于 crates/result-import/src/lib.rs
pub trait ResultParser: Send + Sync {
    /// 解析文件为 ImportedPsm 列表
    fn parse(&self, path: &Path, unimod: &UnimodDb) -> Result<Vec<ImportedPsm>>;

    /// 检测文件是否为此格式
    fn detect(path: &Path) -> bool where Self: Sized;
}
```

### 5.2 custom_json 解析器

**输入格式**（hela.json 示例）：
```json
{
  "sequence": "AADLLDDVSQK",
  "charge": 2,
  "modify": [[3, 35]],
  "rt": 12.648,
  "precursor_mz": 588.2993,
  "raw_title": "hela_SILAC_DIA_350_1000_Rep1",
  "protein_names": ["sp|P12345|TEST_HUMAN"]
}
```

**解析逻辑**：
- `serde_json` 反序列化整个数组
- `modify: [[pos, unimod_id]]` → 遍历每对，调用 `UnimodDb::to_modification(id, pos, sequence)`
- `rt * 60.0` → `rt_sec`
- `raw_title` → `raw_name`
- `protein_names` → `protein_accessions`
- `score`: 此格式无打分字段，设为 `None`

### 5.3 DIA-NN parquet 解析器

**输入列**（report.parquet 关键列）：
- `Modified.Sequence`: `_AAAC(UniMod:4)DM(UniMod:35)K_` — 带修饰的序列
- `Precursor.Charge`: 电荷态
- `Precursor.Mz`: 母离子 m/z
- `RT`: 保留时间（分钟）
- `Q.Value`: FDR q-value
- `Run`: raw 文件名
- `Protein.Names`: 蛋白质名称

**解析逻辑**：
- 用 `arrow` + `parquet` crate 读取
- 正则 `\(UniMod:(\d+)\)` 提取修饰 ID 和位置
- 去除 `_` 前后缀和 `(UniMod:XX)` 得到纯序列
- `RT * 60.0` → `rt_sec`
- `Q.Value` → `score`
- 可选：按 `filter_qvalue` 过滤（默认 0.01）

### 5.4 pFind .spectra 解析器（预留）

- pFind 结果自带 scan number → 不需要 scan 匹配步骤
- 接口已定义，等拿到样例文件后实现
- `detect()` 检查文件头部特征

## 6. Scan Matching Algorithm

### 6.1 算法流程

```
输入: Vec<ImportedPsm> + mzml_dir: PathBuf + rt_tolerance_sec: f64

1. 按 raw_name 分组 ImportedPsm
2. 对每个 raw_name:
   a. 打开 mzml_dir/raw_name.mzML（利用 reader cache）
   b. 预扫描所有 MS2 → Vec<Ms2Info { scan, rt_sec, isolation_window }>
   c. 按 rt_sec 排序
3. 对每个 ImportedPsm:
   a. 二分查找 |ΔRT| < rt_tolerance_sec 的候选 MS2
   b. 过滤: precursor_mz ∈ [target_mz - lower_offset, target_mz + upper_offset]
   c. 多候选 → 选 |ΔRT| 最小的
   d. 无候选 → matched_scan = None, 记录 unmatched 原因
4. 生成 MatchReport
```

### 6.2 DDA vs DIA

- **DDA**: isolation_window 窄（~2 Da），通常精确一对一匹配
- **DIA**: isolation_window 宽（~25 Da），同一 MS2 可被多个 PSM 匹配到——这是 DIA 的正常行为

### 6.3 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `rt_tolerance_sec` | 30.0 | RT 匹配容差（秒）。DIA cycle time ~3-4s，30s 足够宽容 |

无 m/z tolerance 参数——直接使用 mzML 中 MS2 的 isolation_window 范围。

### 6.4 Ms2Info 结构

```rust
/// 预扫描 mzML 时收集的 MS2 信息
struct Ms2Info {
    scan_number: u32,
    rt_sec: f64,
    /// isolation window: (target_mz, lower_offset, upper_offset)
    isolation_window: (f64, f64, f64),
}
```

## 7. Integration with Existing Tools

### 7.1 转换为标准 SearchResult

匹配后的 `ImportedPsm` 转换为 `core::Psm`：

| ImportedPsm 字段 | → core::Psm 字段 | 来源 |
|-------------------|-------------------|------|
| `matched_scan` | `spectrum_scan` | scan 匹配 |
| `sequence` | `peptide_sequence` | 直接映射 |
| `modifications` | `modifications` | Unimod 查表 |
| `charge` | `charge` | 直接映射 |
| `precursor_mz` | `precursor_mz` | 直接映射 |
| — | `calculated_mz` | 从 sequence + charge + mods 计算 |
| — | `delta_mass_ppm` | `(precursor_mz - calculated_mz) / calculated_mz × 1e6` |
| `score` | `score` | DIA-NN 有 Q.Value；无打分设为 0.0 |
| `protein_accessions` | `protein_accessions` | 直接映射 |
| — | `is_decoy` | `false`（外部结果不含 decoy） |
| — | `q_value` | DIA-NN 的 Q.Value；其他为 None |

Peptide / Protein 聚合复用 `report` crate 现有逻辑。

### 7.2 engine_info

```rust
EngineInfo {
    name: "imported".to_string(),
    version: format!("{}", format_name),  // "custom_json" / "diann_parquet" / "pfind_spectra"
    supported_features: vec![],
}
```

### 7.3 缓存

生成的 `SearchResult` 存入与 `run_search` 共用的 `run_results_cache`，使得：
- `annotate_spectrum(run_id, scan_number)` — 直接可用
- `extract_xic(run_id, scan_number)` — 直接可用
- `generate_summary(run_id)` — 直接可用
- `export_results(run_id)` — 直接可用

## 8. MCP Tool Interface

### 8.1 import_search_results

```json
{
  "name": "import_search_results",
  "description": "Import external search results (DIA-NN, custom JSON, pFind) and match to mzML scans. Returns a run_id for use with annotate_spectrum, extract_xic, and generate_summary.",
  "inputSchema": {
    "type": "object",
    "required": ["result_file", "mzml_dir"],
    "properties": {
      "result_file": {
        "type": "string",
        "description": "Path to external search result file (.json, .parquet, .spectra)"
      },
      "format": {
        "type": "string",
        "enum": ["auto", "custom_json", "diann_parquet", "pfind_spectra"],
        "default": "auto",
        "description": "Result file format. 'auto' detects from extension."
      },
      "mzml_dir": {
        "type": "string",
        "description": "Directory containing mzML files. File association: raw_name + '.mzML'"
      },
      "unimod_path": {
        "type": "string",
        "description": "Path to unimod.xml. If not provided, uses builtin modification database."
      },
      "rt_tolerance_sec": {
        "type": "number",
        "default": 30.0,
        "description": "RT tolerance in seconds for scan matching."
      },
      "filter_qvalue": {
        "type": "number",
        "default": 0.01,
        "description": "Q-value threshold for filtering (DIA-NN). PSMs with Q.Value > threshold are excluded."
      },
      "run_filter": {
        "type": "string",
        "description": "Optional: only import PSMs from this specific run/raw_title."
      }
    }
  }
}
```

### 8.2 返回值

```json
{
  "run_id": "uuid-string",
  "match_report": {
    "total_psms": 53824,
    "matched": 51200,
    "unmatched": 2624,
    "median_rt_delta_sec": 1.2,
    "max_rt_delta_sec": 28.5,
    "per_file": { ... }
  },
  "imported_psm_count": 51200,
  "unique_peptides": 12345,
  "protein_count": 3456,
  "raw_files": ["hela_SILAC_DIA_350_1000_Rep1", ...]
}
```

## 9. Code Structure

```
crates/
  result-import/              ← 新 lib crate
    Cargo.toml
    src/
      lib.rs                  ← ImportedPsm, ImportResult, MatchReport, ResultParser trait, import()
      custom_json.rs          ← CustomJsonParser: hela.json 格式
      diann.rs                ← DiannParser: DIA-NN report.parquet 格式
      pfind.rs                ← PFindParser: pFind .spectra 格式（预留骨架）
      unimod.rs               ← UnimodDb: XML 解析 + 内置精简表 + ID→Modification 转换
      scan_matcher.rs         ← ScanMatcher: RT + isolation_window 匹配算法
      converter.rs            ← ImportedPsm → core::Psm / SearchResult 转换

  mcp-server/
    src/tools.rs              ← 新增 import_search_results tool handler
```

## 10. Dependencies

| Crate | 用途 | 备注 |
|-------|------|------|
| `arrow` + `parquet` | 读取 DIA-NN parquet | 编译较重（~15MB），但是 DIA-NN 格式刚需 |
| `quick-xml` | 解析 unimod.xml | 项目已有（spectrum-io 使用） |
| `regex` | 解析 DIA-NN Modified.Sequence | `(UniMod:\d+)` 提取 |
| `serde_json` | 解析 custom_json | 项目已有 |

## 11. Error Handling

| 场景 | 错误类型 | 处理方式 |
|------|----------|----------|
| 文件不存在 | `FileNotFound` | 返回路径 + 建议检查 |
| 格式不匹配 | `FormatDetectionFailed` | 列出支持的格式 |
| Unimod ID 未知 | `UnknownUnimodId(u32)` | 返回 ID + 建议检查 unimod.xml |
| mzML 找不到 | `MzmlNotFound(raw_name)` | 列出目录中可用的 mzML 文件 |
| 0 匹配 | 不 panic | 返回 MatchReport，unmatched = total |
| Parquet 缺列 | `MissingColumn(name)` | 列出缺失列 + 期望列 |

## 12. Testing Strategy

| 层级 | 测试内容 | 方式 |
|------|----------|------|
| 单元 | UnimodDb 解析 + 查表 | 内置精简 XML fixture |
| 单元 | custom_json 解析（含修饰转换） | 5 条 PSM 的 JSON fixture |
| 单元 | DIA-NN parquet 解析 | 小型 parquet fixture（程序生成） |
| 单元 | scan 匹配算法（DDA/DIA） | 构造已知 RT + isolation_window 的 mock MS2 |
| 单元 | ImportedPsm → Psm 转换 | 验证 calculated_mz、delta_mass_ppm |
| 集成 | import → cache → annotate 全流程 | fixture JSON + 现有 mzML test fixture |

## 13. Scope Boundaries

### 本次实现
- ✅ `result-import` lib crate（完整）
- ✅ `custom_json` 解析器
- ✅ `diann_parquet` 解析器
- ✅ `pfind_spectra` 解析器（预留骨架 + trait 实现）
- ✅ `UnimodDb`（XML 解析 + 内置精简表）
- ✅ Scan 匹配算法（DIA-aware）
- ✅ MCP tool `import_search_results`
- ✅ 与现有 annotate_spectrum / extract_xic / generate_summary 集成

### 不在范围
- ❌ MaxQuant / Spectronaut 等其他格式支持
- ❌ 批量标注 UI（逐个 PSM 调用现有工具）
- ❌ 反库（decoy）生成——外部结果假设已做 FDR 过滤
