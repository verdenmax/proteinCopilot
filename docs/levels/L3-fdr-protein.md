# L3 — FDR + 蛋白推断子系统

承接 [L2](L2-architecture.md)。本篇深入两个相邻 crate：`fdr`（target-decoy FDR 与 q-value）与 `protein-inference`（肽到蛋白的推断），讲清它们如何把搜索引擎产出的 PSM 收敛成带统计置信度的肽与蛋白结论。

## 1. 职责与位置

子系统夹在搜索与报告之间，纯确定性、不调用 LLM：

- 上游 `search-engine`：产出 `Vec<Psm>`，每条带 `score`、`is_decoy`，可选 `q_value`。SimpleSearch 在 `simple_engine.rs` 中当 `decoy_strategy != None` 时已调用 `calculate_fdr` 给 PSM 打 `q_value`，并 `retain(|p| !p.is_decoy)` 删除诱饵。
- 本子系统：PSM 级 / 肽级 FDR、肽到蛋白映射、parsimony 最小蛋白集、razor 肽分配、蛋白级 FDR、序列覆盖率。
- 下游 `report`：按 1% FDR 汇总与导出（psm/peptide/protein 三级）。

诱饵有两段生命周期：建库阶段 `decoy.rs` 的 `generate_decoys` 产出 `REV_` 蛋白；PSM 级 FDR 打分后 SimpleSearch 会 `retain` 删除诱饵 PSM。蛋白级 picked FDR 仍依赖 `REV_` 蛋白组与 target 配对竞争，若推断输入已无诱饵证据（无 decoy 胜者），则全部 target 记 `q = 0.0`。

数据流（在 MCP 工具 `infer_proteins` 中编排，纯 ASCII）：

```
search-engine
   |  Vec<Psm> (score, q_value, is_decoy)
   v
build_peptide_protein_map  -> run_parsimony -> assign_razor_peptides
   -> calculate_protein_fdr -> calculate_coverage (可选, 需 FASTA)
   |
   v
Vec<ProteinGroup> (q_value, coverage)  -> report (1% FDR 汇总/导出)
```

## 2. 模块边界

`crates/fdr/src/`：

- `calculation.rs` — 核心 `calculate_fdr`：PSM 级 target-decoy + q-value。
- `peptide_fdr.rs` — 肽级 FDR：`extract_unique_peptides` 取每序列最佳分，`calculate_peptide_fdr` 复用内核。
- `protein_fdr.rs` — 蛋白级 picked target-decoy 配对竞争 `calculate_protein_fdr`。
- `decoy.rs` — `generate_decoys`：Reverse / Shuffle 诱饵库生成。
- `error.rs` — `FdrError` 及到 `CoreError` 的转换。

`crates/protein-inference/src/`：

- `mapper.rs` — `build_peptide_protein_map` 建肽与蛋白二部图、I/L 归一、按 `q_value_threshold` 过滤。
- `parsimony.rs` — `run_parsimony`：不可区分分组 -> 子集吸收 -> 贪心覆盖 -> unique/shared 分类。
- `razor.rs` — `assign_razor_peptides`：共享肽按证据归属唯一蛋白组。
- `coverage.rs` — `calculate_coverage`：基于 FASTA 计算序列覆盖率。
- `error.rs` — `InferenceError`。

`fdr/lib.rs` 对外平铺 `calculate_fdr`、`calculate_peptide_fdr`、`calculate_protein_fdr` 及对应结果类型等公共 API；`protein-inference/lib.rs` 仅声明 `mapper/parsimony/razor/coverage/error` 五个模块，由编排层按序组合。

## 3. 关键数据结构

PSM 是输入单元（`core/search_result.rs`），蛋白组是输出单元（`core/protein_group.rs`）：

```rust
pub struct Psm {                  // 截取关键字段
    pub peptide_sequence: String,
    pub score: f64,               // 越大越好
    pub q_value: Option<f64>,     // FDR 未算时为 None
    pub protein_accessions: Vec<String>,
    pub is_decoy: bool,
}

pub struct ProteinGroup {
    pub leader_accession: String,        // 组代表（成员字母序最前）
    pub member_accessions: Vec<String>,  // 不可区分蛋白
    pub peptides: Vec<String>,
    pub unique_peptides: Vec<String>,    // 仅此组拥有
    pub razor_peptides: Vec<String>,     // 分得的共享肽
    pub score: f64,
    pub q_value: Option<f64>,            // 蛋白级 FDR
    pub coverage: Option<f64>,           // 0.0~1.0
    pub is_decoy: bool,
}
```

诱饵标记统一由 `core::util::is_decoy_accession` 判定（前缀 `REV_/SHUF_/DECOY_/REVERSED_`）。FDR 内部用轻量 `ScoredPsm { index, score, is_decoy }` 承载排序，算完再按 `index` 把 `(index, q)` 映射回原对象；三级 FDR 共用这一条打分路径，避免重复实现。

## 4. FDR 伪代码

`calculate_fdr(psms: &[ScoredPsm]) -> Result<Vec<(usize, f64)>, FdrError>` 是三级 FDR 的共用内核：

```
1. 守卫：空 -> NoPsms；任一 score 非有限 -> InvalidScore；无 decoy -> NoDecoyHits
2. 按 score 降序：sorted.sort_by(|a,b| b.score.total_cmp(&a.score))
3. 自上而下累积：每条 decoy++ 或 target++，
   fdr = (decoys / targets).min(1.0)
4. tie 拉平：对每段等分数极大区间 [lo,hi)，
   令该段全部 = raw_fdrs[hi-1]（区间末值，已计入同分所有 target+decoy）
5. 向后单调：q[i] = q[i].min(q[i+1])
6. 输出 (原始 index, q)；统计 q <= 0.01 的条数
```

第 4、5 步是要点：先把同分组拉平到“区间末”的 FDR，消除输入次序影响；再从尾向头取 min，让 q-value 随分数下降单调不减。两步真实代码：

```rust
let tie_fdr = raw_fdrs[hi - 1];
for fdr in &mut raw_fdrs[lo..hi] { *fdr = tie_fdr; }
// ...
for i in (0..q_values.len().saturating_sub(1)).rev() {
    q_values[i] = q_values[i].min(q_values[i + 1]);
}
```

肽级：`extract_unique_peptides` 每序列保留最高分 PSM（decoy 标记取最佳 PSM 的）并按序列排序，再包成 `ScoredPsm` 走同一内核。1% 是统计与默认过滤口径（`q <= 0.01`）；推断入口的过滤阈值 `q_value_threshold: Option<f64>` 可配。注意两处 decoy 语义不同：`mapper` 判定肽 decoy 用“所有蛋白皆 decoy”规则，而 `extract_unique_peptides` 直接取最佳 PSM 的 `is_decoy`。

## 5. 蛋白推断伪代码

```
build_peptide_protein_map(psms, q_value_threshold)
  - 有 q_value 且 > 阈值则丢弃（无 q_value 的保守保留）
  - normalize_il：I -> L 合并等价肽
  - 肽为 decoy 当且仅当其所有蛋白皆 decoy
  - 无任何 target 蛋白 -> NoTargetProteins

run_parsimony(map)
  1. group_indistinguishable：肽集相同的蛋白并为一组（BTreeSet 作键）
  2. remove_subsets：肽集严格子集被并入超集（按肽数降序、成员字母序定向）
  3. greedy_set_cover：循环挑“新覆盖肽最多”的组，直到全覆盖
  4. build_protein_groups：肽只属一个选中组 -> unique，否则 shared
     leader = members[0]（字母序最前）；按 score 降序、leader 升序排序

assign_razor_peptides(&mut groups, map)
  - 共享肽（出现在 >1 组、且未作为某组 unique）归属“证据最强”的组
  - tie-break：unique 数 -> score -> leader 字母序
  - 写回各组 razor_peptides

calculate_protein_fdr(groups)   // picked target-decoy
  - 用 REV_ 前缀把 target 组与 decoy 组按基础 accession 配对
  - 组内竞争：score >= 者胜，平分判 target 胜
  - 未配对 target 自动胜；未配对 decoy 作为 decoy 胜者入池
  - 无 decoy 胜者：全部 target 记 q = 0.0
  - 否则对胜者列表跑 calculate_fdr，再把 q 映射回 target 胜者
  - 按 score 降序、leader 升序输出

calculate_coverage(&mut groups, fasta)   // 可选
  - 取 leader 的 FASTA 序列，I/L 归一后标记被肽覆盖的残基
  - 把 peptides 与 razor_peptides 去重合并，用滑动 find 标记所有出现位置（重叠只计一次）
  - coverage = 覆盖残基 / 序列长度；缺序列或空序列 -> None
```

razor 的“证据最强”即 unique 肽数最多的组，该步仅在组数 > 1 时执行，否则直接返回空映射。

greedy 挑组的确定性 tie-break 核心：

```rust
a_cover.cmp(&b_cover)                          // 新覆盖肽更多者优先
    .then_with(|| a.score.partial_cmp(&b.score).unwrap_or(Ordering::Equal))
    .then_with(|| b.members[0].cmp(&a.members[0])) // 字母序更前者胜
```

## 6. 确定性要点

相等分数与 HashMap 迭代序都可能破坏可复现，本子系统逐处消除：

- FDR 排序用 `total_cmp`（对 f64 全序，含 -0/NaN），等分数 tie 拉平到“区间末值”，与输入次序无关（见 `fdr_order_independence_with_ties` 测试）。
- `extract_unique_peptides` 输出按序列排序；`group_indistinguishable` 用排序的 `BTreeSet` 键、成员排序、组按成员序排。
- `remove_subsets` 按肽数降序 + `members[0]` 升序，保证严格子集总被并入同一超集。
- greedy / razor / 蛋白 FDR 的 tie-break 都以字母序最前的 leader 收尾；蛋白 FDR 输出 `score 降序 then leader_accession 升序`。
- 诱饵 `Shuffle` 用固定种子 `seed_from_u64(42)`，两次生成一致。
- 分数比较分两档：parsimony / razor 内用 `partial_cmp(..).unwrap_or(Equal)`（此处分数已有限），PSM 级 FDR 排序用更严格、对 NaN/-0 也成全序的 `total_cmp`。

## 7. 错误处理

库代码不 `unwrap/expect`，全走 `Result`：

- `FdrError`：`NoPsms` / `NoDecoyHits` / `InvalidScore`；`From<FdrError> for CoreError` 落到 `CoreError::ValidationError { context, detail, suggestion }`，附带可执行建议（如“无 decoy 时改用 Reverse 策略”）。
- `InferenceError`：`NoPsms` / `NoTargetProteins` / `FdrFailed(String)` / `NoFastaProvided` / `InvalidPeptide(String)`。
- 优雅降级：编排层若 `calculate_protein_fdr` 失败，记 `warn` 并返回不带 q-value 的组而非中断；`calculate_coverage` 对缺 FASTA 的蛋白置 `coverage = None` 且只告警。

---

回到 [README](README.md) 选择其它层级或子系统。
