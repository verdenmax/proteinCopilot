# Spectrum Annotation & Visualization — 设计文档

> **日期**：2026-03-31
> **状态**：待实施
> **范围**：search-engine (annotate 模块), report (visualize 模块), mcp-server (annotate_spectrum tool)

---

## 1. 问题陈述

用户希望对某一张谱图进行单独分析，查看其与候选肽段的碎片离子匹配详情。支持两种场景：
1. 查看已有搜索结果中某个 PSM 的匹配细节
2. 手动指定一条肽段序列，对某张谱图进行匹配并展示结果

## 2. 设计目标

| # | 目标 | 验收标准 |
|---|------|---------|
| G1 | 结构化标注数据 | annotate() 返回每个峰的标注信息（离子类型、delta ppm） |
| G2 | 交互式 HTML 可视化 | 生成自包含 HTML，浏览器打开可悬停看精确数值 |
| G3 | 肽段序列解析图 | HTML 上方显示序列 + b/y 离子覆盖图 |
| G4 | 两种输入模式 | 从已有 PSM 查看 或 手动指定肽段匹配 |
| G5 | 复用现有逻辑 | 复用 matching.rs 的 b/y 离子生成和 tolerance 检查 |

## 3. 架构

### 3.1 数据结构（search-engine/src/annotate.rs）

```rust
/// 单个峰的标注信息
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct AnnotatedPeak {
    pub mz: f64,
    pub intensity: f64,
    pub annotation: Option<IonAnnotation>,
}

/// 离子标注
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct IonAnnotation {
    pub ion_type: IonType,
    pub ion_number: u32,
    pub theoretical_mz: f64,
    pub delta_mz: f64,
    pub delta_ppm: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema)]
pub enum IonType { B, Y }

/// 理论离子及其匹配状态
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct TheoreticalIon {
    pub ion_type: IonType,
    pub number: u32,
    pub theoretical_mz: f64,
    pub matched: bool,
    pub matched_mz: Option<f64>,
    pub delta_ppm: Option<f64>,
}

/// 完整的谱图标注结果
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumAnnotation {
    pub scan_number: u32,
    pub retention_time_sec: f64,
    pub peptide_sequence: String,
    pub charge: i32,
    pub precursor_mz: f64,
    pub theoretical_mz: f64,
    pub delta_mass_ppm: f64,
    pub score: f64,
    pub matched_ions: u32,
    pub total_ions: u32,
    pub protein_accessions: Vec<String>,
    pub peaks: Vec<AnnotatedPeak>,
    pub b_ions: Vec<TheoreticalIon>,
    pub y_ions: Vec<TheoreticalIon>,
    pub modifications: Vec<Modification>,
}
```

### 3.2 annotate() 函数

```rust
pub fn annotate_spectrum(
    spectrum: &Spectrum,
    peptide_sequence: &str,
    charge: i32,
    precursor_tolerance: &MassTolerance,
    fragment_tolerance: &MassTolerance,
    fixed_modifications: &[Modification],
) -> Result<SpectrumAnnotation, SearchEngineError>
```

内部逻辑：
1. 计算肽段理论质量（含固定修饰），算出 theoretical_mz
2. 检查 precursor m/z 是否在 tolerance 内
3. 调用 `generate_b_ions()` / `generate_y_ions()` 生成理论碎片离子
4. 对每个理论离子，在实验 m/z 数组中做 binary search 找匹配
5. 对每个实验峰，记录是否被标注为某个 b/y 离子
6. 构建并返回 `SpectrumAnnotation`

复用 `matching.rs` 中已有的：
- `generate_b_ions()` / `generate_y_ions()` — 理论离子生成
- `within_tolerance()` — 容差判断
- 需要将这些函数从 `pub(crate)` 提升为 `pub`

### 3.3 HTML 可视化渲染（report/src/visualize.rs）

```rust
pub fn render_annotation_html(
    annotation: &SpectrumAnnotation,
    output_path: &Path,
) -> Result<(), ReportError>
```

生成自包含 HTML 文件，结构：

**上半部分 — 肽段序列解析图：**
- 氨基酸序列水平排列
- 上方：b 离子编号（b1, b2, ...），匹配的红色，未匹配灰色
- 下方：y 离子编号（y1, y2, ...），匹配的蓝色，未匹配灰色
- 类似学术论文中的碎片离子覆盖图

**下半部分 — 谱图：**
- X 轴：m/z，Y 轴：相对强度（归一化到 100%）
- 灰色竖线：未匹配实验峰
- 红色竖线 + 标签：匹配的 b 离子（标注 `b3`）
- 蓝色竖线 + 标签：匹配的 y 离子（标注 `y7`）
- JS 交互：悬停显示 m/z、intensity、ion type、delta ppm

**信息面板：**
- Scan, RT, Charge, Peptide, Protein
- Score, Matched ions / Total ions, Δ mass ppm

**实现方式：**
- HTML 模板存放在 `crates/report/templates/annotation.html`
- Rust 通过 `include_str!` 加载模板
- 将 `SpectrumAnnotation` 序列化为 JSON 嵌入 `<script>` 标签
- 模板中的 JS 读取 JSON 数据，使用 Canvas/SVG 渲染谱图
- 零外部依赖（不引入 D3.js 等库），纯原生 JS + SVG

### 3.4 MCP Tool

新增 `annotate_spectrum` tool：

**输入（两种模式合一）：**
```rust
struct AnnotateSpectrumInput {
    // 模式 1：从已有搜索结果
    run_id: Option<String>,

    // 模式 2：手动指定
    file_path: Option<String>,
    peptide_sequence: Option<String>,
    charge: Option<i32>,
    fixed_modifications: Option<Vec<Modification>>,

    // 共用
    scan_number: u32,
    output_path: Option<String>,
    fragment_tolerance: Option<MassTolerance>,
}
```

**逻辑：**
- 若提供 `run_id`：从 RunCache 中取 SearchResult，找到 scan_number 对应的 PSM，用 PSM 的 peptide/charge/mods 做标注
- 若提供 `file_path` + `peptide_sequence` + `charge`：读取谱图，用指定参数做标注
- 默认 fragment_tolerance：0.02 Da
- 默认 output_path：`./annotation_scan{N}.html`

**输出：**
```rust
struct AnnotateResult {
    output_path: String,
    scan_number: u32,
    peptide_sequence: String,
    charge: i32,
    score: f64,
    matched_ions: u32,
    total_ions: u32,
    delta_mass_ppm: f64,
    protein_accessions: Vec<String>,
    message: String,
}
```

### 3.5 Agent 指令更新

在 `proteomics-search.agent.md` 的 tools 列表添加 `annotate_spectrum`。

新增使用场景：
```
谱图标注：
  - 用户说"看一下 scan 1234 的匹配情况"
    → 调用 annotate_spectrum(run_id=xxx, scan_number=1234)
    → 告知用户"标注文件已生成，请在浏览器中打开 xxx.html"
    → 基于 score/matched_ions 给出简短解读

  - 用户说"用 PEPTIDEK 去匹配 scan 100"
    → 调用 annotate_spectrum(file_path=xxx, scan_number=100, peptide_sequence="PEPTIDEK", charge=2)
    → 展示匹配结果
```

## 4. 文件变更

| 操作 | 文件 | 说明 |
|------|------|------|
| 修改 | `search-engine/src/matching.rs` | 将 `generate_b_ions`, `generate_y_ions`, `within_tolerance` 提升为 `pub` |
| 新建 | `search-engine/src/annotate.rs` | SpectrumAnnotation 数据结构 + annotate_spectrum() |
| 修改 | `search-engine/src/lib.rs` | 导出 annotate 模块 |
| 新建 | `report/templates/annotation.html` | 自包含 HTML 模板（JS + SVG 渲染） |
| 新建 | `report/src/visualize.rs` | render_annotation_html() |
| 修改 | `report/src/lib.rs` | 导出 visualize 模块 |
| 修改 | `mcp-server/src/tools.rs` | annotate_spectrum MCP tool |
| 修改 | `.github/agents/proteomics-search.agent.md` | 添加 annotate_spectrum 工具和使用场景 |

## 5. 错误处理

| 场景 | 处理 |
|------|------|
| scan_number 在文件中不存在 | INVALID_PARAMS，"scan N not found in file" |
| run_id 中找不到该 scan 的 PSM | INVALID_PARAMS，"no PSM found for scan N in this search" |
| peptide_sequence 含非标准氨基酸 | INVALID_PARAMS，"sequence contains non-standard residue" |
| output_path 目录不存在 | 自动创建目录，失败则 IO_ERROR |
| 谱图无 precursor 信息 | INVALID_PARAMS，"spectrum has no precursor information" |

## 6. 测试策略

| 测试 | 内容 |
|------|------|
| annotate.rs 单元测试 | 构造 Spectrum + 已知肽段 → 验证 b/y 离子标注正确 |
| annotate.rs 边界测试 | 空谱图、无匹配、全匹配 |
| visualize.rs 测试 | 生成 HTML → 验证文件存在、包含关键标签（SVG/canvas、JSON 数据） |
| MCP tool 测试 | 模式 1 和模式 2 的参数验证 |
| e2e 测试 | 搜索 → 取第一个 PSM → annotate → 验证 HTML 生成 |

## 7. 不包含

- 多谱图批量标注（单张为主，循环调用即可）
- Mirror plot（两谱图对比）
- 交互式谱图编辑
- PDF 导出
