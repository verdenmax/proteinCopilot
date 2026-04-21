# Entrapment Analysis v2 — 编辑距离 + 替换类型标注

> **目标**：减少 L4 漏检，将因长度不同而被误判为"真陷阱"的高同源 trap PSM 正确归入 L2/L3。
> **次要目标**：在 L2 内标注替换类型（Q/K、二肽、通用 near-isobaric）作为参考信息。

---

## 1. 问题背景

v1 使用 Hamming distance（等长比较），仅搜索 `by_length[len]` 桶。以下场景被完全遗漏：

| 场景 | 示例 | 当前结果 | 预期结果 |
|------|------|---------|---------|
| 等质量二肽替换 | N↔GG (trap 7AA vs target 8AA) | L4（不等长） | L2 (Δm=0) |
| 编辑距离 1-2 的 indel 同源 | PEPTIDEK vs PEPTIDK (len diff=1) | L4 | L2/L3 |
| 多残基替换质量互消 | 两个替换 Δm 互相抵消 (len diff=2) | L4 | L2 |

**核心改动**：将 Hamming distance 替换为 Levenshtein edit distance，搜索范围从 `len` 扩展到 `len±2`。

---

## 2. 设计决策记录

| 决策 | 选择 | 理由 |
|------|------|------|
| 分级体系 | 保持 L0-L4 不变 | 最小改动，已有生态兼容 |
| Q/K 和二肽处理 | L2 内加 `substitution_type` 字段标注 | 不影响分级逻辑，仅参考 |
| 性能策略 | k-mer 预筛 + Levenshtein | 无损过滤，减少 90%+ 候选 |
| 等质量表范围 | 仅 N↔GG, Q↔AG (2 对) | 二肽↔二肽由 Hamming+Δm=0 自动捕获 |
| charge 信息 | 可用时记录，不影响分级 | 分级基于质量差，不依赖 charge |
| 核心目标 | 减少漏检 | substitution_type 为附加参考 |
| 配置字段重命名 | `delta_mz_threshold_da` → `delta_mass_threshold_da` | 修正命名误导（mz 暗示质荷比） |

---

## 3. 数据结构变更

### 3.1 新增 SubstitutionType（types.rs）

```rust
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SubstitutionType {
    /// L0/L1/L4 — 无氨基酸替换或无匹配
    None,
    /// I↔L 同分异构体（L1 时自动设置）
    LIIsomer,
    /// Q↔K 替换 (Δm ≈ 36.4 mDa)
    QKSubstitution,
    /// 等质量二肽替换 (N↔GG, Q↔AG)
    IsobaricDipeptide {
        single_residue: char,  // 'N' or 'Q'
        dipeptide: String,     // "GG" or "AG"
    },
    /// 其他 |Δm| < threshold 的 near-isobaric 替换
    NearIsobaric,
    /// |Δm| ≥ threshold，可区分的同源物
    Distinguishable,
}
```

### 3.2 ClassifiedPsm 扩展

```rust
pub struct ClassifiedPsm {
    // ... 所有现有字段保持不变 ...

    /// 替换类型标注（v2 新增，参考信息）
    pub substitution_type: SubstitutionType,

    /// 编辑距离（v2 新增，替代语义上的 hamming distance）
    /// 对于等长匹配，edit_distance == hamming_distance
    pub edit_distance: Option<u32>,

    /// 对齐详情（v2 新增）
    /// 格式："Q7→K" (单替换), "ins:G@5" (插入), "del:A@3" (删除)
    pub alignment_detail: Option<String>,
}
```

现有 `mismatches` 和 `delta_mass_da` 字段保留不变（向后兼容）。`mismatches` 对于等长匹配仍等于 Hamming distance；对于不等长匹配设为 `None`（由 `edit_distance` 提供）。

### 3.3 L2 判定逻辑变更

**v1**：`mismatches == 1 && |Δm| < threshold` → L2

**v2**：`edit_distance <= max_edit_distance && |Δm| < threshold` → L2
       `edit_distance <= max_edit_distance && |Δm| >= threshold` → L3

不再限制 "mismatches == 1" 的硬编码，改为通用的 edit_distance ≤ 配置阈值。

---

## 4. TargetDigestIndex 扩展（digest.rs）

### 4.1 新增 k-mer 倒排索引

```rust
pub struct TargetDigestIndex {
    // 现有字段保持不变
    pub by_length: HashMap<usize, Vec<TargetPeptide>>,
    pub exact_set: HashSet<String>,
    pub normalized_set: HashSet<String>,
    pub exact_to_protein: HashMap<String, String>,
    pub normalized_to_original: HashMap<String, (String, String)>,

    // v2 新增
    kmer_index: HashMap<u64, Vec<u32>>,   // kmer_hash → peptide_ids
    all_peptides: Vec<TargetPeptide>,      // id-indexed flat array
    kmer_k: usize,                         // k-mer 长度
}
```

### 4.2 k-mer 构建逻辑

```
k = min_peptide_length / (max_edit_distance + 1)
  = 6 / (2 + 1) = 2

// 对于 k=2，每条长度 12 的肽段提取 11 个 2-mer
// k-mer 用 FxHash 哈希为 u64 以节省内存
```

### 4.3 新增 find_similar() 方法

```rust
impl TargetDigestIndex {
    /// 查找编辑距离 ≤ max_edit_dist 的所有 target 肽段。
    /// 使用 k-mer 预筛 + Levenshtein 计算。
    pub fn find_similar(
        &self,
        query: &str,
        max_edit_dist: u16,
        len_tolerance: usize,  // 默认 2
    ) -> Vec<SimilarityMatch> {
        // 1. 提取 query 的所有 k-mers
        // 2. 在 kmer_index 中查找候选 peptide_ids（取并集）
        // 3. 去重 + 过滤长度范围 [len-tol, len+tol]
        // 4. 对候选计算 Levenshtein edit distance
        // 5. 保留 edit_distance <= max_edit_dist
        // 6. 计算对齐后的 delta_mass
    }
}

pub struct SimilarityMatch {
    pub target_peptide: String,
    pub target_protein: String,
    pub edit_distance: u32,
    pub delta_mass_da: f64,
    pub alignment_detail: String,
    pub substitution_type: SubstitutionType,
}
```

### 4.4 鸽巢定理保证

对于 edit_distance ≤ d 的两条序列，它们必然共享至少一个长度为 k 的公共子串（k = ⌊n/(d+1)⌋）。因此 k-mer 预筛是**数学上无损的**——不会漏掉任何真实匹配。

---

## 5. classify_single() 新流程（similarity.rs）

```
classify_single(psm, group, index, config):
  if group != Trap → return L4 (不变)

  // Phase 1: 精确匹配（不变）
  if index.has_exact(peptide) → L0
  if index.has_normalized(peptide) → L1, substitution_type=LIIsomer

  // Phase 2: 模糊匹配（v2 升级）
  candidates = index.find_similar(peptide, config.max_mismatches, config.len_tolerance)
  // 内部：k-mer 预筛 → Levenshtein → 长度过滤

  if candidates.is_empty() → L4

  best = candidates.min_by(edit_distance, then |delta_mass|)

  // Phase 3: 分级 + 标注
  substitution_type = categorize_substitution(peptide, best)
  if |best.delta_mass| < config.delta_mass_threshold_da → L2
  else → L3
```

### 5.1 categorize_substitution() 函数

```rust
fn categorize_substitution(trap: &str, best: &SimilarityMatch) -> SubstitutionType {
    let len_diff = trap.len() as i32 - best.target_peptide.len() as i32;

    // 1. 等长 + edit_distance=1 → 检查 Q↔K
    if len_diff == 0 && best.edit_distance == 1 {
        if is_qk_pair_at_diff(trap, &best.target_peptide) {
            return SubstitutionType::QKSubstitution;
        }
    }

    // 2. 长度差 1 → 检查等质量二肽
    if len_diff.abs() == 1 {
        if let Some((single, dipeptide)) = check_isobaric_dipeptide(trap, &best.target_peptide) {
            return SubstitutionType::IsobaricDipeptide {
                single_residue: single,
                dipeptide,
            };
        }
    }

    // 3. 通用分类（threshold = config.delta_mass_threshold_da）
    if best.delta_mass_da.abs() < config.delta_mass_threshold_da {
        SubstitutionType::NearIsobaric
    } else {
        SubstitutionType::Distinguishable
    }
}
```

### 5.2 等质量二肽替换表

```rust
const ISOBARIC_DIPEPTIDES: &[(char, &str)] = &[
    ('N', "GG"),  // 114.04293 Da
    ('Q', "AG"),  // 128.05858 Da
];
```

检测逻辑：对齐 trap 和 target 后，查看 indel 位置是否对应表中的替换对。

---

## 6. Levenshtein 实现

### 6.1 标准算法

使用经典 Wagner-Fischer 算法，O(mn) 时间，O(min(m,n)) 空间（单行优化）。

### 6.2 delta_mass 计算（对齐后）

Levenshtein 对齐产生编辑操作序列。对于每个操作：
- **替换 A→B**：Δm += residue_mass(B) - residue_mass(A)
- **插入 B**：Δm += residue_mass(B)
- **删除 A**：Δm -= residue_mass(A)

需要回溯对齐路径（不能只用距离值），所以实现时保留完整 DP 矩阵用于回溯（仅对通过 k-mer 预筛的候选，数量很少）。

### 6.3 alignment_detail 格式

- 替换：`"D0→N"` （位置0，D替换为N）
- 插入：`"ins:G@5"` （位置5插入G）
- 删除：`"del:A@3"` （位置3删除A）
- 多操作用逗号分隔：`"D0→N,ins:G@5"`

---

## 7. 配置变更（config.rs）

```rust
pub struct SimilarityConfig {
    /// 最大编辑距离（替代原 max_mismatches 的语义）
    #[serde(default = "default_max_mismatches")]
    pub max_mismatches: u16,  // 保留字段名以兼容 v1 YAML

    /// 质量差阈值 (Da)，L2 vs L3 的分界线
    /// v1 名称 delta_mz_threshold_da 已废弃但仍接受
    #[serde(
        default = "default_delta_mass_threshold_da",
        alias = "delta_mz_threshold_da"
    )]
    pub delta_mass_threshold_da: f64,

    // 现有字段不变
    pub require_tryptic_ends: bool,
    pub max_missed_cleavages: u32,

    // v2 新增
    /// 长度容差：搜索 len±len_tolerance 范围内的 target 肽段
    #[serde(default = "default_len_tolerance")]
    pub len_tolerance: usize,  // 默认 2

    /// 启用等质量二肽检测
    #[serde(default = "default_true")]
    pub enable_dipeptide_check: bool,  // 默认 true

    /// 启用 Q/K 近等质量标注
    #[serde(default = "default_true")]
    pub enable_qk_detection: bool,  // 默认 true
}
```

兼容性：v1 的 `delta_mz_threshold_da` 通过 `#[serde(alias)]` 继续接受。

---

## 8. 输出变更

### 8.1 TSV 输出（output.rs）

`classified.tsv` 新增 3 列（追加在末尾，不改变现有列顺序）：

| 列名 | 类型 | 说明 |
|------|------|------|
| `substitution_type` | String | None/LIIsomer/QKSubstitution/IsobaricDipeptide/NearIsobaric/Distinguishable |
| `edit_distance` | u32? | 编辑距离（等长时 = mismatches） |
| `alignment_detail` | String? | 对齐详情 |

### 8.2 HTML 报告（report.rs + template）

- 表格新增 `substitution_type` 列
- Delta-mass histogram 的 hover 增加 substitution_type 信息
- 小于 0.1 Da 的 delta_mass 值显示为 mDa（如 "36.4 mDa"）
- 新增 "Substitution Type" 子饼图（仅 L2 内部）

### 8.3 MCP Tool 输出

`ClassifyEntrapmentOutput` 和 `FindSimilarTargetsOutput` 增加对应字段。

---

## 9. 测试策略

### 9.1 单元测试（新增 ~20 个）

**Levenshtein 函数**：
- 等长序列 → 结果等于 Hamming distance
- 插入/删除 → 正确的 edit distance
- 空序列处理
- delta_mass 计算（替换/插入/删除各场景）
- alignment_detail 格式

**k-mer 索引**：
- k-mer 提取正确性
- 预筛后候选集包含所有真实匹配（无损验证）
- 空索引 / 单条目 / 大量条目

**替换类型检测**：
- Q↔K 检测：PEPTQDE vs PEPTKDE → QKSubstitution
- 二肽检测：PEPNDE vs PEPGGDE → IsobaricDipeptide('N', "GG")
- 通用 near-isobaric：DGFLLDGFPR vs NGFLLDGFPR → NearIsobaric
- 区分：大 Δm → Distinguishable

**classify_single 升级**：
- 不等长匹配正确进入 L2/L3（之前为 L4）
- 等长匹配结果与 v1 一致（回归测试）
- substitution_type 正确标注

### 9.2 回归测试

所有现有 56 个测试必须继续通过。v1 的等长匹配行为不变。

### 9.3 集成测试

新增 mini fixture：
- 包含 Q/K 替换对的 FASTA + PSM
- 包含 N↔GG 二肽替换的 FASTA + PSM
- 包含 len±1/±2 的 indel 同源物
- 端到端验证分级结果

---

## 10. 文件变更清单

| 文件 | 变更类型 | 说明 |
|------|---------|------|
| `crates/entrapment-analysis/src/types.rs` | 修改 | 新增 SubstitutionType 枚举，ClassifiedPsm 新增 3 字段 |
| `crates/entrapment-analysis/src/digest.rs` | 修改 | TargetDigestIndex 新增 k-mer 索引 + find_similar() |
| `crates/entrapment-analysis/src/similarity.rs` | 修改 | classify_single() 升级，新增 Levenshtein + categorize_substitution |
| `crates/entrapment-analysis/src/config.rs` | 修改 | SimilarityConfig 新增 3 字段，重命名 delta_mz |
| `crates/entrapment-analysis/src/output.rs` | 修改 | TSV 新增 3 列 |
| `crates/entrapment-analysis/src/report.rs` | 修改 | HTML 表格/图表新增 substitution_type |
| `crates/entrapment-analysis/templates/entrapment_report.html` | 修改 | 表格列 + mDa 显示 + 子饼图 |
| `crates/entrapment-analysis/src/levenshtein.rs` | **新增** | Levenshtein 算法 + 对齐回溯 + delta_mass 计算 |
| `crates/mcp-server/src/tools.rs` | 修改 | 输出结构体新增字段 |
| `crates/entrapment-analysis/Cargo.toml` | 可能修改 | 如需 rustc-hash (FxHashMap) 依赖 |

---

## 11. 性能预期

| 操作 | v1 | v2 |
|------|----|----|
| 索引构建 | O(N) | O(N × avg_len) — k-mer 提取 |
| 单条 PSM 分类 | O(N/buckets) Hamming | O(k-mer 候选数 × len²) Levenshtein |
| 内存 | ~50MB (50万 peptides) | ~80MB (+k-mer 索引) |

对于典型场景（~1000 trap PSMs，~50万 target peptides），预计 v2 总时间 < 10 秒。
