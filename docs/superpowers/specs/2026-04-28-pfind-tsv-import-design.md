# pFind TSV Import — Design Spec

> **Date:** 2026-04-28
> **Status:** Approved
> **Scope:** 支持 pFind 结果 TSV 文件的导入，用于谱图标注和 entrapment 分析

---

## 1. 问题描述

pFind 搜索引擎输出的结果经后处理生成标准 TSV 文件（如 `intersection_non_human.tsv`），包含 PSM 级别的鉴定结果。当前系统的 `import_search_results` 和 `classify_entrapment_hits` 两个 MCP Tool 均不支持该格式。需要在两个 crate 中实现 pFind TSV 解析器，使用户可以：

1. 导入 pFind 结果 → `annotate_spectrum` / `extract_xic` 做谱图标注
2. 导入 pFind 结果 → `classify_entrapment_hits` 做 entrapment 分析

## 2. pFind TSV 格式定义

### 2.1 列定义

| # | 列名 | 类型 | 示例 | 说明 |
|---|------|------|------|------|
| 1 | `FileName` | string | `20190830_HF_ZHW_hela_SILAC_DDIA_500_550_2Da_Rep1` | 原始文件名（无扩展名） |
| 2 | `PeptideSequence` | string | `HNDLDDVGK` | 肽段序列（纯氨基酸，无修饰标注） |
| 3 | `Modifications` | string | `10,Carbamidomethyl[C];` | pFind 修饰格式，可为空 |
| 4 | `PepMass` | f64 | `1011.462138` | 肽段中性质量 |
| 5 | `PredRT` | f64 | `15.122` | 预测保留时间（分钟） |
| 6 | `CleavageType` | int | `3` | 酶切类型 |
| 7 | `ProNCTerm` | int | `0` | N/C 端标记 |
| 8 | `Proteins` | string | `sp\|P50475\|SYAC_RAT/` | 蛋白质 accession，尾部带 `/` |
| 9 | `MH+` | f64 | `1012.470123` | [M+H]⁺ 质量 |
| 10 | `Charge` | i32 | `2` | 电荷态 |
| 11 | `ScanNo` | u32 | `9911` | 谱图扫描号（1-based） |
| 12 | `RawScore` | f64 | `15.184` | 原始打分 |
| 13 | `DeltaMassPPM` | f64 | `0.701` | 质量偏差 (ppm) |
| 14 | `DeltaRT(Min)` | f64 | `-0.167` | 保留时间偏差（分钟） |
| 15 | `FinalScore` | f64 | `4.58717e-05` | 最终打分 |
| 16 | `QValue` | f64 | `0` | q-value |

### 2.2 格式识别标志

pFind TSV 通过 header 行中同时存在以下列名识别：
- `PeptideSequence`
- `ScanNo`
- `FileName`

只要 header 中包含这三个列名，即判定为 pFind TSV 格式。

### 2.3 修饰格式解析

pFind 修饰字符串格式为 `pos,Name[Residue];` 的分号分隔列表：

```
10,Carbamidomethyl[C];
1,Carbamidomethyl[C];5,Carbamidomethyl[C];
0,Acetyl[ProteinN-term];
5,Oxidation[M];
```

解析规则：
- 以 `;` 分割（忽略末尾空项）
- 每项以 `,` 分割为 `(position, name_with_residue)`
- `name_with_residue` 格式为 `Name[Residue]`，提取名称和残基
- 位置是 **1-based**（pFind 约定），需转为 0-based 用于内部表示
- 特殊情况：`ProteinN-term` 作为残基时，position 为 N 端修饰

### 2.4 precursor_mz 计算

已知 `MH+`（[M+H]⁺ 质量）和 `Charge`：

```
precursor_mz = (MH+ + (charge - 1) * PROTON_MASS) / charge
```

其中 `PROTON_MASS = 1.007276` Da。

### 2.5 Protein accession 清理

pFind 输出的 Protein 字段尾部带 `/`，需去除：
- `sp|P50475|SYAC_RAT/` → `sp|P50475|SYAC_RAT`

多个蛋白质用 `/` 分隔：
- `sp|P50475|SYAC_RAT/sp|Q12345|TEST_HUMAN/` → `["sp|P50475|SYAC_RAT", "sp|Q12345|TEST_HUMAN"]`

## 3. 架构设计

### 3.1 影响的 Crate

| Crate | 变更 | 说明 |
|-------|------|------|
| `result-import` | 新增 `pfind_tsv.rs` | 实现 `ResultParser` trait，输出 `ImportedPsm` |
| `entrapment-analysis` | 新增 `loader/pfind_tsv.rs` | 输出 `UnifiedPsm` |
| `mcp-server` | 更新 `tools.rs` | 新增 `pfind_tsv` 格式选项，更新 auto 检测 |

### 3.2 result-import crate

新增 `pfind_tsv.rs`：

```rust
pub struct PFindTsvParser;

impl ResultParser for PFindTsvParser {
    fn parse(&self, path: &Path, unimod: &UnimodDb) -> Result<Vec<ImportedPsm>, ResultImportError>;
}

pub fn detect(path: &Path) -> bool; // 读 header 行检测
```

关键设计决策：
- **ScanNo 直接使用**：pFind TSV 已包含 scan number，设置 `matched_scan = Some(scan_no)`，跳过 RT scan matching
- **Score 使用 FinalScore**：作为主打分
- **RawScore 保留**：存入 `ImportedPsm` 的扩展字段或忽略（MVP 先忽略）

`ImportFormat` 枚举新增 `PFindTsv` variant。

### 3.3 entrapment-analysis crate

新增 `loader/pfind_tsv.rs`：

```rust
pub fn load_pfind_tsv(path: &Path) -> Result<Vec<UnifiedPsm>, EntrapmentError>;
```

`ResultFormat` 枚举新增 `PFindTsv` variant。

pFind TSV 列名到 `UnifiedPsm` 字段的映射：
- `PeptideSequence` → `peptide`
- `Charge` → `charge`
- `ScanNo` → `scan_number`
- `FileName` → `spectrum_file`
- `Proteins` → `protein_ids`（去除尾部 `/`）
- `QValue` → `q_value`
- `PredRT` → `retention_time`（分钟）
- `MH+` + `Charge` → `precursor_mz`（计算得到）
- `Modifications` → `modifications`（pFind 格式解析）

### 3.4 格式自动探测

两个 crate 的格式检测均需更新：

**result-import:**
- `.tsv` 扩展名时，读取 header 行检测是否为 pFind TSV
- 匹配 → `ImportFormat::PFindTsv`
- 不匹配 → 当前行为（报错，因为之前不支持 `.tsv`）

**entrapment-analysis:**
- `.tsv` 扩展名时，读取 header 行检测是否为 pFind TSV
- 匹配 → `ResultFormat::PFindTsv`
- 不匹配 → fallback 到 `ResultFormat::GenericTsv`

### 3.5 MCP Server 集成

**import_search_results:**
- 格式选项新增 `pfind_tsv`
- `auto` 检测支持 pFind TSV
- pFind TSV 导入跳过 scan matching（已有 ScanNo），直接构建结果
- 仍需 `mzml_dir` 参数（供下游 `annotate_spectrum` 使用）

**classify_entrapment_hits:**
- 格式选项新增 `pfind_tsv`
- `auto` 检测支持 pFind TSV

### 3.6 修饰解析共享

pFind 修饰格式 `pos,Name[Residue];` 的解析逻辑需在两个 crate 中使用。放在 `result-import` crate 中，entrapment-analysis 通过自己的解析逻辑处理（两个 crate 输出类型不同：`Modification` vs `UnifiedPsm.modifications`）。

## 4. 数据流

### 4.1 import_search_results 路径

```
pFind TSV → PFindTsvParser.parse()
         → Vec<ImportedPsm>  (matched_scan 已填充)
         → skip scan matching
         → build_search_result()
         → SearchResult (cached)
         → annotate_spectrum / extract_xic
```

### 4.2 classify_entrapment_hits 路径

```
pFind TSV → load_pfind_tsv()
         → Vec<UnifiedPsm>
         → EntrapmentAnalyzer.classify_all()
         → classified results
```

## 5. 测试策略

- 将 `output/intersection_non_human.tsv` 的前 10 行作为 fixture 文件
- 单元测试：修饰解析、protein 清理、precursor_mz 计算
- 集成测试：完整文件解析、格式自动检测
- 边界测试：空修饰字段、单电荷、多蛋白质

## 6. 范围限制

- 仅支持当前示例文件的 TSV 列定义
- 不支持 `.spectra` 格式（未来单独实现）
- 不支持列名大小写不敏感匹配（严格匹配 pFind 输出的大小写）
- RawScore 字段在 MVP 中不保留到最终结果
