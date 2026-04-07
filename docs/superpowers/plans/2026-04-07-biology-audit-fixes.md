# Biology Audit Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 5 biology audit findings — add missing SearchParams fields, make digest configurable, track protein-terminal peptide positions, and document future work.

**Architecture:** Add 3 new fields to `SearchParams` (with `#[serde(default)]` for backward compatibility), add 2 boolean fields to `DigestedPeptide`, plumb configurable lengths through `digest()`, add documentation comments and future work tracking.

**Tech Stack:** Rust, serde, schemars, cargo test

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/core/src/search_params.rs` | Modify | Add 3 new fields + validation |
| `crates/search-engine/src/digest.rs` | Modify | Add terminal flags, use configurable lengths |
| `crates/search-engine/src/simple_engine.rs` | Modify | Pass peptide length params to `digest()` |
| `crates/param-recommend/src/preset.rs` | Modify | Set new field defaults in all presets |
| `crates/param-recommend/src/rules.rs` | Modify | Set new field defaults in recommendations |
| `crates/dia-extraction/src/detection.rs` | Modify | Add documentation comment |
| `tasks/001-mvp-proteomics-search-platform.md` | Modify | Add future work items |
| `docs/development.md` | Modify | Add biology conventions section |

---

### Task 1: Add `max_variable_modifications`, `min_peptide_length`, `max_peptide_length` to SearchParams

**Files:**
- Modify: `crates/core/src/search_params.rs:198-271`

- [ ] **Step 1: Write tests for new fields and validation**

Add to the existing `mod tests` block in `crates/core/src/search_params.rs`:

```rust
#[test]
fn max_variable_modifications_default() {
    // JSON without max_variable_modifications should deserialize with default 3
    let json = r#"{
        "enzyme": "Trypsin",
        "missed_cleavages": 2,
        "fixed_modifications": [],
        "variable_modifications": [],
        "precursor_tolerance": {"value": 10.0, "unit": "Ppm"},
        "fragment_tolerance": {"value": 0.02, "unit": "Da"},
        "database_path": "test.fasta",
        "decoy_strategy": "Reverse"
    }"#;
    let params: SearchParams = serde_json::from_str(json).unwrap();
    assert_eq!(params.max_variable_modifications, 3);
    assert_eq!(params.min_peptide_length, 7);
    assert_eq!(params.max_peptide_length, 50);
}

#[test]
fn peptide_length_validation_min_gt_max() {
    let mut params = valid_params();
    params.min_peptide_length = 20;
    params.max_peptide_length = 10;
    assert!(params.validate().is_err());
}

#[test]
fn peptide_length_validation_zero_min() {
    let mut params = valid_params();
    params.min_peptide_length = 0;
    assert!(params.validate().is_err());
}

#[test]
fn max_variable_modifications_validation() {
    let mut params = valid_params();
    params.max_variable_modifications = 11;
    assert!(params.validate().is_err());
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-core -- search_params`
Expected: FAIL — fields don't exist yet

- [ ] **Step 3: Add new fields to SearchParams struct**

In `crates/core/src/search_params.rs`, add to the `SearchParams` struct after `acquisition_mode`:

```rust
    /// Maximum number of variable modifications per peptide (default: 3).
    /// Limits combinatorial explosion during variable modification enumeration.
    #[serde(default = "default_max_variable_modifications")]
    pub max_variable_modifications: u32,

    /// Minimum peptide length in residues (default: 7).
    /// Peptides shorter than this are excluded from search results.
    #[serde(default = "default_min_peptide_length")]
    pub min_peptide_length: u32,

    /// Maximum peptide length in residues (default: 50).
    /// Peptides longer than this are excluded from search results.
    #[serde(default = "default_max_peptide_length")]
    pub max_peptide_length: u32,
```

Add the default functions above the struct:

```rust
fn default_max_variable_modifications() -> u32 { 3 }
fn default_min_peptide_length() -> u32 { 7 }
fn default_max_peptide_length() -> u32 { 50 }

/// Maximum allowed value for `max_variable_modifications`.
const MAX_VARIABLE_MODS_LIMIT: u32 = 10;
```

- [ ] **Step 4: Add validation errors**

Add to `SearchParamsError` enum:

```rust
    /// max_variable_modifications exceeds limit.
    #[error("max_variable_modifications must be <= {max}, got {actual}")]
    TooManyVariableMods {
        actual: u32,
        max: u32,
    },

    /// min_peptide_length is zero.
    #[error("min_peptide_length must be >= 1, got 0")]
    ZeroPeptideLength,

    /// min_peptide_length > max_peptide_length.
    #[error("min_peptide_length ({min}) must be <= max_peptide_length ({max})")]
    InvalidPeptideLengthRange {
        min: u32,
        max: u32,
    },
```

- [ ] **Step 5: Add validation logic**

In `SearchParams::validate()`, add after the modification validation block:

```rust
        if self.max_variable_modifications > MAX_VARIABLE_MODS_LIMIT {
            return Err(SearchParamsError::TooManyVariableMods {
                actual: self.max_variable_modifications,
                max: MAX_VARIABLE_MODS_LIMIT,
            });
        }
        if self.min_peptide_length == 0 {
            return Err(SearchParamsError::ZeroPeptideLength);
        }
        if self.min_peptide_length > self.max_peptide_length {
            return Err(SearchParamsError::InvalidPeptideLengthRange {
                min: self.min_peptide_length,
                max: self.max_peptide_length,
            });
        }
```

- [ ] **Step 6: Update `valid_params()` test helper**

Add the 3 new fields to the `valid_params()` function in tests:

```rust
    fn valid_params() -> SearchParams {
        SearchParams {
            // ... existing fields ...
            acquisition_mode: None,
            max_variable_modifications: 3,
            min_peptide_length: 7,
            max_peptide_length: 50,
        }
    }
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `cargo test -p protein-copilot-core -- search_params`
Expected: ALL PASS

- [ ] **Step 8: Fix compilation in downstream crates**

Any code constructing `SearchParams` directly (presets, tests, simple_engine) will now fail to compile. Fix all instances by adding the 3 new fields with their defaults. Search with:

```bash
cargo build --workspace 2>&1 | grep "missing field"
```

Fix each site by adding:
```rust
max_variable_modifications: 3,
min_peptide_length: 7,
max_peptide_length: 50,
```

- [ ] **Step 9: Run full test suite**

Run: `cargo test --workspace`
Expected: ALL PASS

- [ ] **Step 10: Commit**

```bash
git add -A
git commit -m "feat(core): add max_variable_modifications, min/max_peptide_length to SearchParams

Addresses biology audit findings:
- max_variable_modifications (default 3): prevents combinatorial explosion
- min_peptide_length (default 7): standard proteomics minimum
- max_peptide_length (default 50): standard proteomics maximum
- All fields use #[serde(default)] for backward compatibility
- Validation: max_var_mods <= 10, min_len >= 1, min <= max

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Make digest() use configurable peptide lengths + add protein-terminal flags

**Files:**
- Modify: `crates/search-engine/src/digest.rs`
- Modify: `crates/search-engine/src/simple_engine.rs`

- [ ] **Step 1: Write tests for configurable length and terminal flags**

Add to `crates/search-engine/src/digest.rs` tests:

```rust
#[test]
fn digest_respects_custom_min_length() {
    // With min_length=8, "PEPTIDEK" (8 chars) should be included
    let peptides = digest_with_length(
        "PEPTIDEKANSTHERPEPTIDERLASTPART",
        "P001",
        &Enzyme::Trypsin,
        0,
        8,
        50,
    );
    let seqs: Vec<&str> = peptides.iter().map(|p| p.sequence.as_str()).collect();
    assert!(seqs.contains(&"PEPTIDEK")); // exactly 8
    // Short fragments < 8 should be excluded
    for p in &peptides {
        assert!(p.sequence.len() >= 8, "too short: {}", p.sequence);
    }
}

#[test]
fn digest_respects_custom_max_length() {
    let peptides = digest_with_length(
        "PEPTIDEKANSTHERPEPTIDERLASTPART",
        "P001",
        &Enzyme::Trypsin,
        0,
        6,
        10,
    );
    for p in &peptides {
        assert!(p.sequence.len() <= 10, "too long: {}", p.sequence);
    }
}

#[test]
fn digest_marks_protein_nterm() {
    let peptides = digest_with_length(
        "PEPTIDEKANSTHERPEPTIDERLASTPART",
        "P001",
        &Enzyme::Trypsin,
        0,
        6,
        50,
    );
    // First peptide should be protein N-terminal
    assert!(peptides[0].is_protein_nterm, "first peptide should be N-term");
    assert!(!peptides[0].is_protein_cterm);
    // Last peptide should be protein C-terminal
    let last = peptides.last().unwrap();
    assert!(last.is_protein_cterm, "last peptide should be C-term");
    assert!(!last.is_protein_nterm);
}

#[test]
fn digest_single_peptide_is_both_terminal() {
    // A short protein that produces only one peptide
    let peptides = digest_with_length("PEPTIDEK", "P001", &Enzyme::Trypsin, 0, 6, 50);
    assert_eq!(peptides.len(), 1);
    assert!(peptides[0].is_protein_nterm);
    assert!(peptides[0].is_protein_cterm);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test -p protein-copilot-search-engine -- digest`
Expected: FAIL — `digest_with_length` and terminal flag fields don't exist

- [ ] **Step 3: Add terminal flags to DigestedPeptide**

```rust
#[derive(Debug, Clone)]
pub struct DigestedPeptide {
    /// Amino acid sequence.
    pub sequence: String,
    /// Protein accession this peptide came from.
    pub protein_accession: String,
    /// Monoisotopic neutral mass (Da).
    pub neutral_mass: f64,
    /// Whether this peptide starts at the protein N-terminus.
    pub is_protein_nterm: bool,
    /// Whether this peptide ends at the protein C-terminus.
    pub is_protein_cterm: bool,
}
```

- [ ] **Step 4: Add `digest_with_length()` function**

```rust
/// Digests a protein sequence with configurable peptide length range.
///
/// `min_length` and `max_length` control accepted peptide lengths (in residues).
/// Also tracks whether each peptide is at the protein N- or C-terminus.
pub fn digest_with_length(
    sequence: &str,
    protein_accession: &str,
    enzyme: &Enzyme,
    missed_cleavages: u32,
    min_length: u32,
    max_length: u32,
) -> Vec<DigestedPeptide> {
    let cleavage_sites = find_cleavage_sites(sequence, enzyme);
    let fragments = split_at_sites(sequence, &cleavage_sites);
    let num_fragments = fragments.len();

    let mut peptides = Vec::new();

    for mc in 0..=(missed_cleavages as usize) {
        for (i, window) in fragments.windows(mc + 1).enumerate() {
            let combined: String = window.concat();
            let len = combined.len() as u32;
            if len >= min_length && len <= max_length {
                if let Some(mass) = peptide_mass(&combined) {
                    let is_nterm = i == 0;
                    let is_cterm = i + mc + 1 == num_fragments;
                    peptides.push(DigestedPeptide {
                        sequence: combined,
                        protein_accession: protein_accession.to_string(),
                        neutral_mass: mass,
                        is_protein_nterm: is_nterm,
                        is_protein_cterm: is_cterm,
                    });
                }
            }
        }
    }

    peptides
}
```

- [ ] **Step 5: Update original `digest()` to delegate**

```rust
pub fn digest(
    sequence: &str,
    protein_accession: &str,
    enzyme: &Enzyme,
    missed_cleavages: u32,
) -> Vec<DigestedPeptide> {
    digest_with_length(sequence, protein_accession, enzyme, missed_cleavages, 6, 50)
}
```

- [ ] **Step 6: Fix all DigestedPeptide construction sites**

Search for `DigestedPeptide {` in the codebase and add the new fields. Main sites:
- The old `digest()` body is now replaced by delegation, so no manual fix needed there.
- Any test code constructing `DigestedPeptide` directly needs the new fields.

- [ ] **Step 7: Update `simple_engine.rs` to use configurable lengths**

In `SimpleSearchEngine::run_search_inner()`, find where `digest()` is called and change to:

```rust
let peptides = digest_with_length(
    &entry.sequence,
    &entry.accession,
    &params.enzyme,
    params.missed_cleavages,
    params.min_peptide_length,
    params.max_peptide_length,
);
```

Import `digest_with_length` at the top of `simple_engine.rs`.

- [ ] **Step 8: Run tests**

Run: `cargo test --workspace`
Expected: ALL PASS

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "feat(search-engine): configurable peptide length + protein-terminal tracking

- digest_with_length() accepts min/max peptide length parameters
- DigestedPeptide now tracks is_protein_nterm/is_protein_cterm
- SimpleSearchEngine uses SearchParams.min/max_peptide_length
- Original digest() preserved as backward-compatible wrapper (6-50)

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Update presets and param-recommend with new fields

**Files:**
- Modify: `crates/param-recommend/src/preset.rs`
- Modify: `crates/param-recommend/src/rules.rs`

- [ ] **Step 1: Update all preset functions**

In `crates/param-recommend/src/preset.rs`, add the 3 new fields to every `SearchParams` construction in all 5 preset functions (`standard_preset`, `phospho_preset`, `tmt_preset`, `open_search_preset`, `silac_preset`):

```rust
max_variable_modifications: 3,
min_peptide_length: 7,
max_peptide_length: 50,
```

For phospho preset specifically, use `max_variable_modifications: 3` (important for limiting combinatorial explosion with Phospho + Oxidation).

- [ ] **Step 2: Update recommendation rules**

In `crates/param-recommend/src/rules.rs`, find where `SearchParams` is constructed in `recommend()` and add the 3 new fields with defaults.

- [ ] **Step 3: Run tests**

Run: `cargo test --workspace`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(param-recommend): set biology defaults in all presets

- All presets: max_variable_modifications=3, min_peptide_length=7, max_peptide_length=50
- Recommendation rules output includes new fields

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Add DIA documentation comment

**Files:**
- Modify: `crates/dia-extraction/src/detection.rs`

- [ ] **Step 1: Add documentation comment**

In `crates/dia-extraction/src/detection.rs`, update the doc comment on `detect_acquisition_mode`:

```rust
/// Detects whether spectra were acquired in DDA or DIA mode based on
/// the median isolation window width of MS2 spectra.
///
/// **Detection logic:** Total window width = `lower_offset + upper_offset`.
/// This correctly handles both symmetric windows (e.g., ±12.5 Da → 25 Da total)
/// and asymmetric windows (e.g., -10/+15 Da → 25 Da total) as defined in the
/// mzML specification (lower/upper offsets from target_mz).
///
/// **Threshold:** Typical DDA isolation windows are 1–3 Da; DIA windows are 10–25 Da.
/// The default threshold of 5 Da is conservative and correctly separates the two modes.
/// Some optimized DIA methods use narrow 5 Da windows — these would be classified as DDA,
/// which is an acceptable limitation for auto-detection. Users can override via `acquisition_mode`.
///
/// Returns `AcquisitionMode::Unknown` when no MS2 spectra exist or none
/// carry isolation window information.
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p protein-copilot-dia-extraction`
Expected: ALL PASS (no logic change)

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "docs(dia-extraction): document DIA detection threshold and asymmetric windows

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Update tasks and docs with future work

**Files:**
- Modify: `tasks/001-mvp-proteomics-search-platform.md`
- Modify: `docs/development.md`

- [ ] **Step 1: Add biology audit future work to tasks**

Append to `tasks/001-mvp-proteomics-search-platform.md` at the end, before any trailing `---`:

```markdown

---

## Biology Audit: 已修复与未来工作

> **审计日期**: 2026-04-07
> **审计范围**: 全部 7 个 crate 的生物学/化学计算正确性

### 已验证正确 ✅

| 检查项 | 结果 |
|--------|------|
| 20 种氨基酸单同位素残基质量 | 与 NIST 标准一致 |
| PROTON_MASS (1.007276 Da) | 正确 |
| WATER_MASS (18.010565 Da) | 正确 |
| ¹³C-¹²C 质量差 (1.003355 Da) | 正确 |
| b 离子公式 (Σ残基，不加水) | 正确 |
| y 离子公式 (Σ残基 + H₂O) | 正确 |
| PPM 计算 (分母为理论值) | 正确 |
| 修饰质量 (CAM/Ox/Phospho/TMT/SILAC) | 与 UniMod 一致 |
| 酶切规则 (Trypsin K/R not before P 等) | 正确 |
| DIA 隔离窗口解读 | 正确 |
| 同位素峰间距与检测 | 正确 |

### 已修复 🔧

| # | 问题 | 修复方案 | 影响 |
|---|------|----------|------|
| BIO-1 | `SearchParams` 缺少 `max_variable_modifications` | 添加字段，默认 3，上限 10 | 防止可变修饰组合爆炸 |
| BIO-2 | 肽段长度 6-50 硬编码 | 添加 `min/max_peptide_length` 到 SearchParams | 标准默认 7-50，可配置 |
| BIO-3 | `DigestedPeptide` 缺少蛋白端位标记 | 添加 `is_protein_nterm/cterm` | 为 ProteinNTerm 修饰枚举做准备 |
| BIO-4 | DIA 检测缺少非对称窗口文档 | 补充详细文档注释 | 说明 5 Da 阈值的合理性与局限 |

### 未来工作（Phase 2+）

| # | 工作项 | 优先级 | 依赖 | 说明 |
|---|--------|--------|------|------|
| FW-1 | **ProteinNTerm 修饰枚举** | 高 | 可变修饰枚举实现 | 搜索时仅对 `is_protein_nterm=true` 的肽段应用 ProteinNTerm 修饰（如 Acetyl） |
| FW-2 | **可变修饰组合枚举** | 高 | — | 实现 `max_variable_modifications` 限制的组合生成算法 |
| FW-3 | **ppm 碎片离子容差** | 中 | — | 预设改用 ppm（20 ppm）代替 Da（0.02 Da），对高分辨 Orbitrap HCD 更准确 |
| FW-4 | **多电荷碎片离子** | 中 | — | 对高电荷母离子（z≥3），生成 b²⁺/y²⁺ 碎片离子提高匹配率 |
| FW-5 | **a/c/x/z 离子系列** | 低 | — | ETD/ECD 碎裂模式需要 c/z 离子；CID 辅助可加 a 离子 |
| FW-6 | **原生 FDR 计算** | 高 | M2.2 | 实现 target-decoy FDR + q-value 单调化，不再依赖外部引擎 |
| FW-7 | **负离子模式** | 低 | — | 当前仅支持正离子 [M+nH]ⁿ⁺，负模式需要 [M-nH]ⁿ⁻ |
```

- [ ] **Step 2: Add biology conventions section to development.md**

Append to `docs/development.md`:

```markdown

## 质谱学生物学约定

以下常数和公式已经过审计验证（2026-04-07），与 NIST/UniMod 标准一致。

### 质量常数 (`crates/search-engine/src/chemistry.rs`)

| 常数 | 值 (Da) | 来源 |
|------|---------|------|
| PROTON_MASS | 1.007276 | NIST |
| WATER_MASS | 18.010565 | H₂O 单同位素 |
| C13_C12_MASS_DIFF | 1.003355 | ¹³C - ¹²C |

### 碎片离子公式

- **b 离子**: `b_n = Σ(residue_1..n)` — 不含水
- **y 离子**: `y_n = Σ(residue_{n+1}..end) + H₂O` — 含水（C 端保留 OH，N 端保留 H）
- **m/z 转换**: `ion_mz = (ion_mass + charge × PROTON_MASS) / charge`
- **当前限制**: 仅生成单电荷碎片（b¹⁺, y¹⁺）

### PPM 计算

```
delta_ppm = (observed - theoretical) / theoretical × 1e6
```

分母始终使用**理论值**（不是观测值）。

### 修饰应用规则

- **固定修饰**: 自动应用到所有目标残基
- **可变修饰**: 组合枚举，受 `max_variable_modifications` 限制（默认 3）
- **N 端修饰**: `AnyNTerm` 应用于所有肽段 N 端；`ProteinNTerm` 仅应用于蛋白质第一条肽段
- **C 端修饰**: 同理

### 酶切规则

| 酶 | 规则 | 异常 |
|----|------|------|
| Trypsin | K/R 后切 | P 前不切 |
| Trypsin/P | K/R 后切 | 无异常 |
| Lys-C | K 后切 | — |
| Glu-C | D/E 后切 | — |
| Asp-N | D **前**切 | — |
| Chymotrypsin | F/W/Y/L 后切 | — |

### DIA 检测

- **隔离窗口宽度** = `lower_offset + upper_offset`（总宽度，非半宽）
- **DIA 判定阈值**: 中位窗口宽度 > 5 Da
- **局限**: 5 Da 窄窗口 DIA 会被误判为 DDA（可通过 `acquisition_mode` 手动指定）
```

- [ ] **Step 3: Verify docs formatting**

Read the modified files to ensure markdown is correct:

```bash
head -20 tasks/001-mvp-proteomics-search-platform.md
tail -60 docs/development.md
```

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "docs: add biology audit results and future work tracking

- tasks/001: audit results table (verified/fixed/future work)
- docs/development.md: mass constants, ion formulas, enzyme rules reference

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
