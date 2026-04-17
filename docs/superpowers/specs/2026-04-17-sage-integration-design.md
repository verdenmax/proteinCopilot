# M2.0 Sage 搜索引擎集成 — 设计文档

> **方案**: 全量库集成（方案 A）  
> **日期**: 2026-04-17  
> **状态**: 已批准，待实施

---

## 1. 目标

将 [Sage](https://github.com/lazear/sage)（Rust 原生蛋白质组学搜索引擎）作为 `sage-core` 库依赖集成到 ProteinCopilot，实现 `SearchEngineAdapter` trait。用户通过 `engine: "Sage"` 参数即可切换引擎，无需安装外部二进制。

## 2. 选择 Sage 的理由

- **Rust 原生**：与 ProteinCopilot 技术栈一致，可作为库直接调用
- **高性能**：rayon 并行，单文件搜索速度与 MSFragger 相当
- **功能完整**：内置 FDR（picked-peptide + KDE PEP）、LFQ、TMT 定量
- **开源 MIT**：无许可证风险
- **API 设计良好**：sage-core 是独立库 crate，适合嵌入式使用

## 3. 架构设计

### 3.1 文件结构

```
crates/search-engine/src/adapters/
├── mod.rs              ← 新增 pub mod sage;
├── pfind.rs            ← 已有（空占位）
└── sage/
    ├── mod.rs          ← SageAdapter struct + SearchEngineAdapter impl
    ├── convert.rs      ← 双向类型转换
    └── config.rs       ← SearchParams → sage Parameters 构建
```

改动的现有文件：
- `crates/search-engine/Cargo.toml` — 新增 `sage-core` git 依赖、`rayon`
- `crates/search-engine/src/adapters/mod.rs` — 加 `pub mod sage;`
- `crates/mcp-server/src/tools.rs` — `run_search` 从硬编码改为 registry 查找
- `crates/core/src/search_result.rs` — `Psm` 加 `extra` 字段
- `crates/core/src/search_params.rs` — `SearchParams` 加 `engine` 字段

### 3.2 SageAdapter 结构体

```rust
pub struct SageAdapter {
    thread_count: usize,
}
```

**无状态设计**：每次 `search()` 内部独立创建 `IndexedDatabase` 和 `Scorer`。与 `SimpleSearchEngine` 保持一致。`IndexedDatabase` 的生命周期限定在 `spawn_blocking` 闭包中，不存在跨请求的引用问题。

### 3.3 搜索流程

```
Phase 1: 读取谱图（spectrum-io，async 上下文）
    input_files → spectrum-io → Vec<Spectrum> → Vec<RawSpectrum>

Phase 2: 构建 DB + 打分（spawn_blocking，rayon 并行）
    FASTA → sage Parameters → IndexedDatabase
    RawSpectrum → ProcessedSpectrum → Scorer::score() → Vec<Feature>
    进度通过 AtomicUsize + tokio interval (500ms) 轮询

Phase 3: FDR（spawn_blocking）
    sage_core::fdr::spectrum_fdr()
    sage_core::fdr::picked_peptide()
    sage_core::fdr::picked_protein()

Phase 4: 结果转换
    Vec<Feature> + IndexedDatabase → Vec<Psm> → SearchResult
```

### 3.4 rayon/tokio 桥接

sage-core 使用 rayon（同步并行），我们的 MCP Server 使用 tokio（异步）。通过 `tokio::task::spawn_blocking` 将整个 sage 搜索推到阻塞线程池：

```rust
let (features, db) = tokio::task::spawn_blocking(move || {
    let db = sage_params.build(fasta_content)?;
    let scorer = Scorer::new(&db, ...);
    let features: Vec<Feature> = raw_spectra.par_iter()
        .flat_map(|spec| {
            let processed = spec.process(...);
            let results = scorer.score(&processed);
            progress_counter.fetch_add(1, Ordering::Relaxed);
            results
        })
        .collect();
    Ok((features, db))
}).await??;
```

进度回调通过 `Arc<AtomicUsize>` + `tokio::time::interval(500ms)` 桥接到 `ProgressCallback`。

## 4. 类型转换层

### 4.1 Spectrum → RawSpectrum

| 我们的字段 | Sage 字段 | 转换 |
|-----------|----------|------|
| `mz_array: Vec<f64>` | `mz: Vec<f32>` | `as f32`（精度损失可接受，sage 内部全部 f32） |
| `intensity_array: Vec<f64>` | `intensity: Vec<f32>` | `as f32` |
| `retention_time_min: f64`（实际存秒） | `scan_start_time: f32`（分钟） | `(sec / 60.0) as f32` |
| `scan_number: u32` | `id: String` | `.to_string()` |
| `precursors[].charge: Option<i32>` | `charge: Option<u8>` | `.map(\|c\| c as u8)`，负值钳位 |
| `precursors[].mz: f64` | `mz: f32` | `as f32` |
| `ms_level: MsLevel` | `ms_level: u8` | `MS1→1, MS2→2` |
| `precursors[].isolation_window` | `isolation_window: Option<Tolerance>` | `Da(lo+hi / 2)` 或 `Ppm(...)` |
| —（无） | `file_id: usize` | 调用方提供 |

**精度损失评估**：m/z 值 1000.0 的 f32 精度约 ±0.00006 Da，远小于典型搜索容差（10 ppm ≈ 0.01 Da）。sage 自身读 mzML 也是 f32，所有内部算法都在 f32 精度下设计。

### 4.2 Modification → ModificationSpecificity

| 我们的 ModPosition | Sage 映射 |
|-------------------|-----------|
| `Anywhere` + residues | 每个 residue 展开为 `Residue(r as u8)` |
| `AnyNTerm` + 无 residue | `PeptideN(None)` |
| `AnyNTerm` + residues | 每个 residue 展开为 `PeptideN(Some(r))` |
| `ProteinNTerm` | `ProteinN(...)` |
| `AnyCTerm` | `PeptideC(...)` |
| `ProteinCTerm` | `ProteinC(...)` |

固定修饰 → `HashMap<ModificationSpecificity, f32>`  
可变修饰 → `HashMap<ModificationSpecificity, Vec<f32>>`

### 4.3 Feature → Psm

| Sage Feature | 我们的 Psm | 转换 |
|-------------|-----------|------|
| `spec_id` | `spectrum_scan` | 解析 scan number（`parse::<u32>()`） |
| `hyperscore` | `score` | 直接赋值 |
| `spectrum_q` | `q_value` | `Some(f as f64)` |
| `charge` | `charge` | `u8 → i32` |
| `expmass` | `precursor_mz` | `(mass + charge * 1.00728) / charge` |
| `calcmass` | `calculated_mz` | 同上公式 |
| `label == -1` | `is_decoy` | `true` |
| `peptide_idx → db` | `peptide_sequence` | 从 `IndexedDatabase` 反查 |
| `peptide_idx → db` | `protein_accessions` | 从 `IndexedDatabase` 反查 |
| `peptide_idx → db` | `modifications` | 从 `IndexedDatabase` 反查并反转换 |
| `delta_mass_ppm` | `delta_mass_ppm` | 直接赋值（或从 expmass/calcmass 计算） |
| `matched_peaks`, `delta_next`, `discriminant_score`, etc. | `extra` | 保存到 `HashMap<String, serde_json::Value>` |

## 5. 现有系统集成改动

### 5.1 SearchParams 新增 `engine` 字段

```rust
pub struct SearchParams {
    // ... 现有字段不变
    /// 指定使用的搜索引擎名称。None 表示使用默认引擎（SimpleSearch）。
    #[serde(default)]
    pub engine: Option<String>,
}
```

### 5.2 Psm 新增 `extra` 字段

```rust
pub struct Psm {
    // ... 现有字段不变
    /// Engine-specific extra fields (e.g., Sage's matched_peaks, delta_next).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub extra: Option<HashMap<String, serde_json::Value>>,
}
```

### 5.3 tools.rs — registry 查找

```rust
// run_search handler 改为：
let engine_name = params.engine.as_deref().unwrap_or("SimpleSearch");
let engine = state.registry.get(engine_name)
    .ok_or(CoreError::SearchEngineError {
        engine: engine_name.to_string(),
        detail: format!("Engine '{}' not registered", engine_name),
        suggestion: format!("Available engines: {:?}",
            state.registry.list_available().iter().map(|e| &e.name).collect::<Vec<_>>()),
    })?;
```

### 5.4 EngineRegistry 初始化

```rust
// McpState::new() 中：
let mut registry = EngineRegistry::new();
registry.register(Box::new(SimpleSearchEngine::new()));
registry.register(Box::new(SageAdapter::default()));
```

## 6. 错误处理

复用现有 `CoreError::SearchEngineError { engine, detail, suggestion }`，不新增错误枚举。

典型错误场景：
- FASTA 文件不存在 → `detail: "FASTA file not found"`, `suggestion: "Check fasta_path"`
- sage-core 构建 DB 失败 → `detail: "Failed to build indexed database"`, `suggestion: "Check FASTA format"`
- 谱图数为 0 → `detail: "No MS2 spectra found"`, `suggestion: "Check input files"`
- charge 值超出 u8 范围 → 钳位并记录 warning（不报错）

## 7. 依赖管理

```toml
# crates/search-engine/Cargo.toml
[dependencies]
sage-core = { git = "https://github.com/lazear/sage.git", rev = "<pinned-commit>" }
rayon = "1"
```

使用 commit hash 锁定版本（当前最新为 v0.15.0-beta.2 附近），避免 breaking change。实施时确认最新稳定 commit。后续 sage 发布新版时手动升级并测试。

## 8. 测试策略

| 层级 | 测试内容 | 方式 |
|------|---------|------|
| 单元测试 | `spectrum_to_raw()` 转换 | 构造 Spectrum → 验证 RawSpectrum 各字段 |
| 单元测试 | `convert_mod()` 修饰映射 | 覆盖 Anywhere/NTerm/CTerm + 多 residue |
| 单元测试 | `feature_to_psm()` 字段映射 | 构造 Feature → 验证 Psm 各字段 |
| 单元测试 | `build_sage_parameters()` 参数构建 | 验证 enzyme/tolerance/modifications 映射 |
| 集成测试 | 小 FASTA + 小 mgf 端到端搜索 | SageAdapter::search() → 验证 PSM 数 > 0 |
| 集成测试 | engine registry 查找 + 切换 | 注册 Simple + Sage，按名查找验证 |
| 集成测试 | FDR q-value 合理性 | 验证 q_value 在 [0, 1]，单调递增 |

sage-core 作为真实依赖参与测试，不 mock——因为它是纯计算库，速度快。

## 9. 不在本次范围内

- LFQ 定量集成（后续 M2.x 任务）
- TMT 定量集成（后续 M2.x 任务）
- Sage 的 RT 预测/校准
- 多引擎结果合并
- DIA 模式特殊处理（Sage 支持 wide-window，但参数映射更复杂）

## 10. 代码量估算

| 文件 | 新增行数 |
|------|---------|
| `adapters/sage/mod.rs` | ~200 行 |
| `adapters/sage/convert.rs` | ~180 行 |
| `adapters/sage/config.rs` | ~80 行 |
| 现有文件改动 | ~40 行 |
| 测试代码 | ~150 行 |
| **合计** | **~650 行** |
