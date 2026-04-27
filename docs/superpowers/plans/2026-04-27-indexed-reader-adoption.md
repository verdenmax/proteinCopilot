# 全局索引读取器采用 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将所有 MCP Server tool 中的 `create_reader()` 替换为 `get_or_create_reader()`（索引 + LRU 缓存），并将 `extract_spectrum_precursors` 从 `read_all()` 改为定向索引读取。同步更新 agent/prompt 文档中过时的性能指导。

**Architecture:** MCP Server 已有 `get_or_create_reader()` 函数（IndexedMzMLReader LRU cache，容量 8）。当前有 7 处调用使用 `create_reader()`（非索引、无缓存），其中 `extract_spectrum_precursors` 加载全部 10 万+ 谱图仅为读取 1 个 MS2 + 若干 MS1。library crate（simple_engine、sage adapter）无法访问 MCP Server 缓存，改用 `create_indexed_reader()`。agent/prompt 文档中 XIC 瓶颈描述已过时（XIC 已优化为 <1s），需更新。

**Tech Stack:** Rust, spectrum-io crate, mcp-server tools, agent.md / prompt.md

---

## Task 1: MCP Server — 替换 create_reader() 为 get_or_create_reader()

**Files:**
- Modify: `crates/mcp-server/src/tools.rs`
  - Line 1354: `read_spectra` tool
  - Line 1397: `recommend_params` tool
  - Line 1741: `run_search` tool (params auto-recommend path)
  - Line 1785: `run_search` tool (DIA guard)
  - Line 2684: `extract_dia_precursors` tool
  - Line 3347: `prepare_search` tool

**变更说明：** 这 6 处都用 `create_reader()`（返回非索引的 MzMLReader/MgfReader），应替换为 `self.get_or_create_reader()` 以复用 LRU 缓存中的 IndexedMzMLReader。注意 `get_or_create_reader()` 接受 `&Path` 并返回 `Result<Arc<dyn SpectrumReader>, ErrorData>`，无需先调用 `detect_format()`。

- [ ] **Step 1: 替换 `read_spectra` (line 1354)**

```rust
// 替换前 (line 1352-1357):
let info = protein_copilot_spectrum_io::detect_format(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let reader = protein_copilot_spectrum_io::create_reader(&info);
let summary = reader
    .read_summary(path)

// 替换后:
let reader = self.get_or_create_reader(path)?;
let summary = reader
    .read_summary(path)
```

- [ ] **Step 2: 替换 `recommend_params` (line 1394-1400)**

```rust
// 替换前:
let path = std::path::Path::new(fp);
let info = protein_copilot_spectrum_io::detect_format(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let reader = protein_copilot_spectrum_io::create_reader(&info);
reader
    .read_summary(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?

// 替换后:
let path = std::path::Path::new(fp);
let reader = self.get_or_create_reader(path)?;
reader
    .read_summary(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?
```

- [ ] **Step 3: 替换 `run_search` params auto-recommend path (line 1739-1744)**

```rust
// 替换前:
let info = protein_copilot_spectrum_io::detect_format(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let reader = protein_copilot_spectrum_io::create_reader(&info);
let summary = reader
    .read_summary(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

// 替换后:
let reader = self.get_or_create_reader(path)?;
let summary = reader
    .read_summary(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
```

- [ ] **Step 4: 替换 `run_search` DIA guard (line 1784-1786)**

```rust
// 替换前:
if let Ok(info) = protein_copilot_spectrum_io::detect_format(first_path) {
    let reader = protein_copilot_spectrum_io::create_reader(&info);
    if let Ok(summary) = reader.read_summary(first_path) {

// 替换后:
if let Ok(reader) = self.get_or_create_reader(first_path) {
    if let Ok(summary) = reader.read_summary(first_path) {
```

注意：此处原代码用 `if let Ok` 容错模式，替换后保持一致。

- [ ] **Step 5: 替换 `extract_dia_precursors` (line 2682-2687)**

```rust
// 替换前:
let info = protein_copilot_spectrum_io::detect_format(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let reader = protein_copilot_spectrum_io::create_reader(&info);
let spectra = reader
    .read_all(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

// 替换后:
let reader = self.get_or_create_reader(path)?;
let spectra = reader
    .read_all(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
```

注意：`extract_dia_precursors` 需要全部谱图（DIA 提取算法要求），因此 `read_all()` 保留。但使用缓存 reader 让后续同文件操作受益。

- [ ] **Step 6: 替换 `prepare_search` (line 3345-3350)**

```rust
// 替换前:
let info = protein_copilot_spectrum_io::detect_format(first_path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let reader = protein_copilot_spectrum_io::create_reader(&info);
let summary = reader
    .read_summary(first_path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

// 替换后:
let reader = self.get_or_create_reader(first_path)?;
let summary = reader
    .read_summary(first_path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
```

- [ ] **Step 7: 构建验证**

```bash
cargo build -p protein-copilot-mcp-server 2>&1 | head -30
```

Expected: 编译成功，无错误。

- [ ] **Step 8: 运行测试**

```bash
cargo test -p protein-copilot-mcp-server 2>&1 | tail -10
```

Expected: 所有测试通过。

- [ ] **Step 9: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "perf: MCP tools 全面采用 get_or_create_reader() 索引缓存

替换 6 处 create_reader() 为 get_or_create_reader()：
- read_spectra: 复用缓存索引
- recommend_params: 同上
- run_search (2处): params 自动推荐 + DIA 安全检测
- extract_dia_precursors: 缓存 reader 供后续操作复用
- prepare_search: 同上

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 2: extract_spectrum_precursors — 从 read_all() 改为定向索引读取

**Files:**
- Modify: `crates/mcp-server/src/tools.rs` (lines 2770-2803, `extract_spectrum_precursors` 函数)

**变更说明：** 当前 `extract_spectrum_precursors` 调用 `create_reader().read_all()` 将整个 mzML 文件加载到内存（100k+ 谱图），仅为找到 1 个目标 MS2 和附近的 MS1。改为使用 `get_or_create_reader()` + `list_scan_meta()` 索引规划 + 定向 `read_spectrum()` 读取目标 MS2 和 RT 附近的 MS1 谱图。

`extract_single_spectrum_precursors(&[Spectrum], scan, extractor)` 接受 `&[Spectrum]`，我们只需传入目标 MS2 + 附近 MS1（而非全部谱图）即可。

- [ ] **Step 1: 实现定向读取逻辑**

替换 `extract_spectrum_precursors` 函数体中的读取逻辑：

```rust
// 替换前 (line 2772-2778):
let path = Path::new(&input.file_path);
let info = protein_copilot_spectrum_io::detect_format(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let reader = protein_copilot_spectrum_io::create_reader(&info);
let spectra = reader
    .read_all(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;

// 替换后:
let path = Path::new(&input.file_path);
let reader = self.get_or_create_reader(path)?;

// 读取目标 MS2
let target_ms2 = reader
    .read_spectrum(path, input.scan_number)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
if target_ms2.ms_level != protein_copilot_core::spectrum::MsLevel::MS2 {
    return Err(mcp_err(
        ErrorCode::INVALID_PARAMS,
        format!("scan {} is not MS2 (ms_level={:?})", input.scan_number, target_ms2.ms_level),
    ));
}
let target_rt = target_ms2.retention_time_sec;

// 用索引查找附近 MS1 scan（±60s RT 范围）
let scan_metas = reader
    .list_scan_meta(path)
    .map_err(|e| mcp_core_err(protein_copilot_core::error::CoreError::from(e)))?;
let rt_min = target_rt / 60.0; // list_scan_meta 返回 rt_min (分钟)
const MS1_RT_WINDOW_MIN: f64 = 1.0; // ±1 分钟
let nearby_ms1_scans: Vec<u32> = scan_metas
    .iter()
    .filter(|m| m.ms_level == 1 && (m.rt_min - rt_min).abs() <= MS1_RT_WINDOW_MIN)
    .map(|m| m.scan_number)
    .collect();

// 定向读取 MS1 谱图
let mut spectra = vec![target_ms2];
for scan_no in &nearby_ms1_scans {
    match reader.read_spectrum(path, *scan_no) {
        Ok(s) => spectra.push(s),
        Err(_) => {} // 跳过读取失败的 scan
    }
}
```

- [ ] **Step 2: 构建验证**

```bash
cargo build -p protein-copilot-mcp-server 2>&1 | head -30
```

Expected: 编译成功。

- [ ] **Step 3: 运行测试**

```bash
cargo test -p protein-copilot-mcp-server 2>&1 | tail -10
```

Expected: 测试通过。

- [ ] **Step 4: Commit**

```bash
git add crates/mcp-server/src/tools.rs
git commit -m "perf: extract_spectrum_precursors 从 read_all 改为索引定向读取

旧版：create_reader().read_all() 加载 100k+ 谱图
新版：get_or_create_reader() + list_scan_meta() 找附近 MS1
     + read_spectrum() 定向读取目标 MS2 + ±1min MS1
     读取量从 100k+ → ~50 个谱图

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 3: search-engine library crate — 采用 create_indexed_reader()

**Files:**
- Modify: `crates/search-engine/src/simple_engine.rs` (lines 294, 324)
- Modify: `crates/search-engine/src/adapters/sage/mod.rs` (line 96)

**变更说明：** Library crate 无法访问 MCP Server 的 `reader_cache`，但可以使用 `create_indexed_reader()` 替代 `create_reader()` 以获得 `.mzML.idx` 磁盘缓存加速。

`simple_engine.rs` 当前对同一文件创建两个 reader（一次 summary，一次 streaming），合并为一个 indexed reader 复用。

- [ ] **Step 1: 合并 simple_engine.rs 的两次文件读取**

```rust
// 替换前 (lines 287-354 — 两个独立 for loop):
// Loop 1: create_reader → read_summary → 累计 ms2_count
// Loop 2: create_reader → for_each_spectrum → 搜索匹配

// 替换后: 合并为单个 loop，每个文件只创建一个 indexed reader
let mut total_ms2_spectra: u64 = 0;
let mut processed_ms2_spectra: u64 = 0;
let mut psms: Vec<Psm> = Vec::new();

// Pre-scan summaries (still need total count for progress)
for file_path in input_files {
    let reader = protein_copilot_spectrum_io::create_indexed_reader(file_path)
        .map_err(|e| {
            let detail = e.to_string();
            diagnostics.fail_stage(&detail);
            diagnostics.set_error(ErrorCategory::InputData, &detail);
            SearchEngineError::IoError { detail }
        })?;
    let summary = reader.read_summary(file_path).map_err(|e| {
        let detail = e.to_string();
        diagnostics.fail_stage(&detail);
        diagnostics.set_error(ErrorCategory::InputData, &detail);
        SearchEngineError::IoError { detail }
    })?;
    total_ms2_spectra += summary.ms2_count;
}

if total_ms2_spectra == 0 {
    diagnostics.fail_stage("No MS2 spectra found in input files");
    diagnostics.set_error(
        ErrorCategory::InputData,
        "No MS2 spectra found in input files",
    );
    return Err(SearchEngineError::NoInputSpectra);
}

report("Reading spectra", 0.15);

for file_path in input_files {
    let reader = protein_copilot_spectrum_io::create_indexed_reader(file_path)
        .map_err(|e| {
            let detail = e.to_string();
            diagnostics.fail_stage(&detail);
            diagnostics.set_error(ErrorCategory::InputData, &detail);
            SearchEngineError::IoError { detail }
        })?;
    // ... handler 和 for_each_spectrum 保持不变 ...
}
```

注意：`for_each_spectrum` 在 IndexedMzMLReader 上委托给 MzMLReader（streaming），所以功能不变。但 reader 本身有 `.mzML.idx` 缓存，且如果后续需要 `read_spectrum()` 等操作，已有索引可用。

- [ ] **Step 2: 替换 sage adapter (line 96)**

```rust
// 替换前:
let reader = protein_copilot_spectrum_io::create_reader(&file_info);

// 替换后:
let reader = protein_copilot_spectrum_io::create_indexed_reader(path)
    .map_err(|e| CoreError::SearchEngineError {
        engine: "Sage".into(),
        detail: format!("Failed to create indexed reader for {}: {}", path.display(), e),
        suggestion: "Check that the input file exists and is a valid mzML/mgf file".into(),
    })?;
```

注意：`create_indexed_reader` 返回 `Result`，需调整错误处理。同时删除原先的 `detect_format()` + `create_reader()` 两步调用。

- [ ] **Step 3: 构建验证**

```bash
cargo build -p protein-copilot-search-engine 2>&1 | head -30
```

- [ ] **Step 4: 运行测试**

```bash
cargo test -p protein-copilot-search-engine 2>&1 | tail -10
```

- [ ] **Step 5: Commit**

```bash
git add crates/search-engine/
git commit -m "perf: search-engine 采用 create_indexed_reader 复用磁盘索引

- simple_engine: create_reader → create_indexed_reader，复用 .mzML.idx 缓存
- sage adapter: 同上

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 4: 更新 agent/prompt 文档 — XIC 性能指导已过时

**Files:**
- Modify: `.github/prompts/spectrum-annotation.prompt.md` (lines 71-102)
- Modify: `.github/agents/proteomics-search.agent.md` (lines 205-218)

**变更说明：** XIC 提取已从"全文件流式扫描 100-150s"优化为"索引规划 + 定向读取 <1s"（`extract_xic_unified`）。文档中的性能瓶颈描述、批量标注耗时估算、超时设置建议均需更新。

- [ ] **Step 1: 更新 spectrum-annotation.prompt.md 性能指导**

替换 lines 71-102 的性能指导部分：

```markdown
## 性能指导

### 内部处理流程

`annotate_spectrum` 对 mzML 文件执行以下步骤：

| 步骤 | 实现 | 耗时 |
|------|------|------|
| 读取目标 scan | `IndexedMzMLReader` — 磁盘索引 `.mzML.idx` O(1) seek | **毫秒级** |
| 离子匹配 | 内存计算 b/y 理论值 vs 实测峰 | **毫秒级** |
| SILAC 重标 scan 查找 | `find_by_rt()` 索引二分查找 | **毫秒级** |
| XIC 提取 | `extract_xic_unified()` — 索引规划 + 定向读取 ~30 scan | **<1s** |

### 关键认知

- **所有步骤均已 O(1)/O(log N) 优化**：`IndexedMzMLReader` 使用磁盘缓存（`.mzML.idx` sidecar），首次打开时 SIMD byte-scan 构建索引并持久化，后续打开毫秒级加载
- **XIC 已完成索引化优化**：`extract_xic_unified()` 通过 `list_scan_meta()` 从内存索引规划读取目标（~30 scan），不再全文件扫描
- **N 个肽段批量标注完全可行**：每个标注 <2s，12 个肽段 <30s（含 MCP 通信开销）
- **reader_cache LRU（容量 8）**：MCP Server 缓存 IndexedMzMLReader 实例，同一文件的连续操作跳过所有索引加载
- **磁盘索引 .mzML.idx**：PCIX v2 二进制格式，46B/entry，记录每个 scan 的 byte_offset / RT / ms_level / isolation_window。首次打开写入，后续秒开

### scan number 获取策略（优先级从高到低）

1. **从搜索结果/导入结果获取（推荐）**：先 `import_search_results` 或 `run_search` → `export_results` 获取 PSM.tsv，一次性拿到所有肽段的 scan number
2. **直接传 scan_number**：如果用户已知 scan number，直接传入。O(1) seek，毫秒级完成
3. **RT 查找（后备）**：`scan_number=0 + retention_time_min` 模式需要索引二分查找，首次打开大文件可能需要数十秒构建索引

### 批量标注建议

- 先用 `export_results` 一次获取所有目标肽段的 scan number、charge、modifications
- 顺序调用 `annotate_spectrum`（服务端 reader LRU 缓存复用索引，每次 <2s）
- 无需特殊超时设置，默认超时即可
```

- [ ] **Step 2: 更新 proteomics-search.agent.md 性能指导**

替换 lines 205-218 的性能指导：

```markdown
### 性能关键：scan number 优先于 RT 查找

**所有标注操作均已索引化优化，单次标注 <2s（含 XIC 提取）。**

1. **有搜索结果时（最快）**：先 `import_search_results` 或 `run_search` 获取 run_id，然后 `export_results` 得到 PSM 列表（含 scan number）。直接用 `scan_number` 调用标注。
2. **手动模式用 scan_number（推荐）**：`IndexedMzMLReader` 有磁盘索引缓存（`.mzML.idx`），**单 scan 查找是 O(1) seek**，毫秒级完成。
3. **RT 查找（后备）**：`scan_number=0 + retention_time_min` 模式需要索引二分查找，仍然很快。

**批量标注**：
- 每次 `annotate_spectrum` 调用 <2s（XIC 通过 `extract_xic_unified()` 索引定向读取 ~30 scan）
- 顺序调用即可，服务端 reader_cache LRU（容量 8）自动复用索引
- 先用 `export_results` 一次性获取所有 scan number，避免逐个通过 RT 查找

### 索引体系

- **内存索引**：`ScanIndex`（`HashMap<scan_number → ScanMeta>`）+ RT 预排序数组
- **磁盘缓存**：`.mzML.idx`（PCIX v2 二进制，46B/entry），记录每个 scan 的 byte_offset / RT / ms_level / isolation_window
- **reader_cache**：MCP Server LRU 缓存（容量 8），缓存 IndexedMzMLReader 实例
- **XIC 规划**：`list_scan_meta()` 从内存索引查询全部 scan 元数据（亚毫秒），筛选后定向 `read_spectrum()` O(1) seek
```

- [ ] **Step 3: Commit**

```bash
git add .github/prompts/spectrum-annotation.prompt.md .github/agents/proteomics-search.agent.md
git commit -m "docs: 更新 agent/prompt XIC 性能指导 — 反映索引优化后 <2s/标注

- spectrum-annotation.prompt.md: XIC 瓶颈描述更新为 <1s
- proteomics-search.agent.md: 批量标注耗时更新，新增索引体系说明
- 移除过时的 ≥600s 超时建议和 25分钟/12肽段估算

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

## Task 5: 更新 coder.agent.md — 新增 spectrum-io 索引读取规范

**Files:**
- Modify: `.github/agents/coder.agent.md`

**变更说明：** coder agent 当前没有关于 spectrum 文件读取的指导。新增一节，指导开发者在编写涉及谱图文件的代码时优先使用索引读取器。

- [ ] **Step 1: 在 coder.agent.md 合适位置添加 spectrum-io 规范**

在文件末尾或"编码规范"部分添加：

```markdown
## spectrum-io 文件读取规范

### 读取器选择（优先级从高到低）

| 场景 | 推荐方式 | 说明 |
|------|---------|------|
| MCP Server tool 中读取谱图 | `self.get_or_create_reader(path)` | LRU 缓存（容量 8）+ IndexedMzMLReader |
| Library crate 中需要 `read_spectrum()` | `create_indexed_reader(path)` | 有 `.mzML.idx` 磁盘缓存 |
| Library crate 中仅需 `for_each_spectrum()` | `create_indexed_reader(path)` | 索引不影响 streaming，但缓存了元数据 |
| 测试代码 | `create_reader(&info)` | 测试用小文件，索引开销不值得 |

### 禁止模式

- ❌ 在 MCP tool 中使用 `create_reader()` — 必须用 `get_or_create_reader()`
- ❌ 使用 `read_all()` 仅为读取单个 scan — 用 `read_spectrum(path, scan_no)` O(1) seek
- ❌ 使用 `for_each_spectrum()` 仅为查询 scan 元数据 — 用 `list_scan_meta()` 或 `list_ms2_meta()` 从内存索引读取

### 索引体系

- `ScanIndex`：内存 HashMap，scan_number → (byte_offset, RT, ms_level, isolation_window)
- `.mzML.idx`：PCIX v2 磁盘缓存（46B/entry），首次打开自动创建
- `reader_cache`：MCP Server LRU 缓存（容量 8），同一文件复用 IndexedMzMLReader
- `list_scan_meta()`：从 ScanIndex 读取全部 scan 元数据（亚毫秒，零 I/O）
```

- [ ] **Step 2: Commit**

```bash
git add .github/agents/coder.agent.md
git commit -m "docs: coder agent 新增 spectrum-io 索引读取规范

指导开发者优先使用 get_or_create_reader() / create_indexed_reader()，
禁止在 MCP tool 中使用 create_reader()，禁止 read_all() 读单 scan。

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
