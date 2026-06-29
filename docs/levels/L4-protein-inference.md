# L4 — protein-inference crate

承接 [L3-fdr-protein](L3-fdr-protein.md)，回溯 [L2](L2-architecture.md)。本篇只聚焦 `crates/protein-inference` 一个 crate 的函数级 API 与算法骨架，签名/常量均按源码核验，完整逻辑以源码为准。

## 1. 用途 + 位置 + 依赖

`protein-copilot-protein-inference` 把 PSM 收敛为蛋白组：建肽到蛋白二部图、parsimony 最小蛋白集、razor 共享肽分配、序列覆盖率。纯确定性、无 LLM、无全局可变状态、无 I/O（覆盖率所需 FASTA 由调用方读入后以 `HashMap` 传入）。蛋白级 picked FDR 不在本 crate，委托依赖 `fdr` 的 `calculate_protein_fdr`。

注意 `lib.rs:10-14` 只声明 5 个 `pub mod`（coverage / error / mapper / parsimony / razor），没有单一 `infer()` 编排函数。"infer 入口"实为调用方 `mcp-server` 的 `infer_proteins` 工具，按序拼接本 crate 的 4 个公共函数加 fdr 的 picked FDR。

位置图（纯 ASCII）:

```
search-engine / result-import: Vec<Psm> (peptide_sequence, score, protein_accessions, is_decoy, q_value?)
   |
   v
mapper::build_peptide_protein_map  -> PeptideProteinMap
parsimony::run_parsimony           -> Vec<ProteinGroup>
razor::assign_razor_peptides       -> HashMap<pep, leader> (并写回 razor_peptides)
fdr::calculate_protein_fdr         -> 写回 q_value (picked, 在 fdr crate)
coverage::calculate_coverage       -> 写回 coverage (可选, 需 FASTA)
   |
   v
report / mcp infer_proteins: 按 1% FDR 汇总导出 Vec<ProteinGroup>
```

依赖（`Cargo.toml`）: `protein-copilot-core`（`Psm` / `ProteinGroup` / `is_decoy_accession`）、`protein-copilot-fdr`（picked FDR）、`thiserror`、`serde`、`tracing`。

## 2. 函数级 API

| 函数 / 类型 | 签名（源码原样） | 位置 |
|---|---|---|
| `build_peptide_protein_map` | `fn build_peptide_protein_map(psms: &[Psm], q_value_threshold: Option<f64>) -> Result<PeptideProteinMap, InferenceError>` | mapper.rs:41 |
| `normalize_il` | `fn normalize_il(sequence: &str) -> String` | mapper.rs:31 |
| `PeptideProteinMap` | `struct { peptide_to_proteins: HashMap<String,HashSet<String>>, protein_to_peptides: HashMap<String,HashSet<String>>, peptide_best_score: HashMap<String,f64>, peptide_is_decoy: HashMap<String,bool> }` | mapper.rs:18 |
| `run_parsimony` | `fn run_parsimony(map: &PeptideProteinMap) -> Result<Vec<ProteinGroup>, InferenceError>` | parsimony.rs:37 |
| `assign_razor_peptides` | `fn assign_razor_peptides(groups: &mut [ProteinGroup], map: &PeptideProteinMap) -> HashMap<String, String>` | razor.rs:21 |
| `calculate_coverage` | `fn calculate_coverage(groups: &mut [ProteinGroup], fasta_sequences: &HashMap<String, String>)` | coverage.rs:25 |
| `ProteinGroup` (core) | `struct { leader_accession, leader_description, member_accessions: Vec<String>, peptides, unique_peptides, razor_peptides: Vec<String>, score: f64, q_value: Option<f64>, coverage: Option<f64>, is_decoy: bool }` | core/protein_group.rs:15 |
| `InferenceError` | `enum { NoPsms, NoTargetProteins, FdrFailed(String), NoFastaProvided, InvalidPeptide(String) }` | error.rs:7 |
| `calculate_protein_fdr` (fdr) | `fn calculate_protein_fdr(groups: &[ProteinGroup]) -> Result<ProteinFdrResult, FdrError>` | fdr/protein_fdr.rs:36 |

内部类型 `IndistinguishableGroup { members: Vec<String>, peptides: HashSet<String>, score: f64 }`（parsimony.rs:18，非 pub）承载分组中间态。`build_peptide_protein_map` 有两类早退错误：空输入 `NoPsms`（mapper.rs:48），全诱饵 `NoTargetProteins`（mapper.rs:128-133）；`run_parsimony` 在 `protein_to_peptides` 为空时返回 `NoPsms`。四个公共函数中 `build_peptide_protein_map` 与 `run_parsimony` 从只读输入产出新值（`PeptideProteinMap` / `Vec<ProteinGroup>`），`assign_razor_peptides` 与 `calculate_coverage` 接收 `&mut [ProteinGroup]` 原地写回。`ProteinGroup.score` 即组内肽最高分（`best_peptide_score` 自 `f64::NEG_INFINITY` fold，parsimony.rs:284），`member_accessions` 恒按字典序排列，首位即 `leader_accession`。

## 3. 算法骨架

整链：肽->蛋白映射(I/L 归一) -> parsimony 最小集 -> razor 分配 -> 蛋白级 picked FDR -> 覆盖率。全程确定性，无随机、无输入顺序依赖。

`build_peptide_protein_map`（mapper.rs:41）:

```
1 可选 q_value_threshold 过滤; 无 q_value 的 PSM 一律保留（保守, 不丢未打分项）
2 normalize_il: 'I' -> 'L'，I/L 等价肽并入同一 key（mapper.rs:31）
3 双向填 peptide_to_proteins / protein_to_peptides; peptide_best_score 取最高分
4 peptide_is_decoy: 仅当全部映射蛋白皆诱饵才记 true（mapper.rs:122, is_decoy_accession）
5 校验至少 1 个 target 蛋白，否则 NoTargetProteins（mapper.rs:128）
```

`run_parsimony`（parsimony.rs:37，4 步）:

```
1 group_indistinguishable: 用 BTreeSet<肽> 做 key，肽集相同的蛋白并为一组（:84）
2 remove_subsets: 肽集按数量降序、首成员字典序定序; 严格子集被吸收进超集（:117,142）
3 greedy_set_cover: 反复挑覆盖未解释肽最多的组, tie-break: cover > score > leader（:180）
4 build_protein_groups: 某肽在选中组里只出现 1 次记 unique，否则 shared（:237,257）
```

最终结果 `sort_by(score 降序 then leader_accession 升序)`（parsimony.rs:70）。leader 恒取组内字典序首成员（members 已排序）。`remove_subsets` 先按肽集大小降序、再按首成员字典序定序（parsimony.rs:121-126），确保同一严格子集总被吸收进同一超集，消除 HashMap 迭代顺序的不确定性。

`assign_razor_peptides`（razor.rs:21）: `groups.len() <= 1` 直接空返回（razor.rs:29）; 建 `pep -> 组下标` 表; 取出现于 >1 组、且非任一组 unique 的共享肽; 每条按 `unique_peptides.len() > score > leader(字典序小者胜)` 分给最佳组; 回写各组 `razor_peptides` 并排序。

蛋白级 picked FDR 由 `fdr::calculate_protein_fdr`（fdr/protein_fdr.rs:36）完成: 以去 `REV_` 前缀的 base accession 配对 target/decoy 组, 组内分高者胜（平局 target 胜），胜者列喂回 PSM 级内核得 q-value; 细节见 [L4-fdr](L4-fdr.md)。注意诱饵判定两处不同口径: mapper 用 `is_decoy_accession`（4 种前缀, core/util.rs:38），picked 配对只 strip `REV_`。

`calculate_coverage`（coverage.rs:25）: 对每组 leader 在 FASTA 序列（I/L 归一）上滑窗 `find` 全部出现位置, 标记 covered 位图, `coverage = covered_count / seq_len`（:73）; leader 不在 FASTA 记 `None`，空序列记 `None`。

## 4. 简化源码片段

run_parsimony 贪心循环（parsimony.rs:189-228）:

```rust
while !uncovered.is_empty() {
    let best_idx = groups.iter().enumerate()
        .filter(|(i, _)| !used[*i])
        .max_by(|(_, a), (_, b)| {
            let a_cover = a.peptides.iter().filter(|p| uncovered.contains(p.as_str())).count();
            let b_cover = b.peptides.iter().filter(|p| uncovered.contains(p.as_str())).count();
            a_cover.cmp(&b_cover)                                         // 1) 覆盖最多
                .then_with(|| a.score.partial_cmp(&b.score).unwrap_or(Ordering::Equal)) // 2) 分高
                .then_with(|| b.members[0].cmp(&a.members[0]))            // 3) leader 字典序小者胜
        }).map(|(i, _)| i);
    match best_idx {
        Some(idx) => { used[idx] = true;
            for pep in &groups[idx].peptides { uncovered.remove(pep.as_str()); }
            selected.push(&groups[idx]); }
        None => break,
    }
}
```

razor tie-break（razor.rs:83-100）:

```rust
let best_idx = indices.iter().copied()
    .max_by(|&a, &b| {
        let (_, a_uniq, a_score, ref a_leader) = group_keys[a];
        let (_, b_uniq, b_score, ref b_leader) = group_keys[b];
        a_uniq.cmp(&b_uniq)                                              // 1) unique 肽多者胜
            .then_with(|| a_score.partial_cmp(&b_score).unwrap_or(Ordering::Equal)) // 2) 分高
            .then_with(|| b_leader.cmp(a_leader))                        // 3) leader 字典序小者胜
    })
    .expect("indices is non-empty for shared peptides");
razor_map.insert(pep.clone(), groups[best_idx].leader_accession.clone());
```

coverage 去重（peptides + razor_peptides 合并进 HashSet，避免重复计数; coverage.rs:52-69）:

```rust
let all_peptides: HashSet<&str> = group.peptides.iter()
    .chain(group.razor_peptides.iter()).map(String::as_str).collect();
for peptide in &all_peptides {
    let normalized_pep = normalize_il(peptide);
    let mut start = 0;
    while let Some(pos) = normalized_fasta[start..].find(&normalized_pep) {
        let abs_pos = start + pos;
        for flag in covered.iter_mut().skip(abs_pos).take(normalized_pep.len()) { *flag = true; }
        start = abs_pos + 1;                                            // +1 抓重叠/重复出现
    }
}
```

## 5. 调用链

- `mcp-server/tools.rs:3527` `infer_proteins` 工具编排（本 crate 内无 orchestrator）: `build_peptide_protein_map`（:3556）-> `run_parsimony`（:3563）-> `assign_razor_peptides`（:3568）-> fdr `calculate_protein_fdr`（:3571，失败仅 warn 降级为无 q-value 组）-> 可选 `calculate_coverage`（:3587，需 `fasta_path`）。
- 上游 PSM 来自 search-engine 或 result-import 的 `SearchResult.psms`；razor 写回的 `razor_peptides` 与 coverage 写回的 `coverage` 都落在同一 `Vec<ProteinGroup>` 上原地修改。
- 下游 `report` 按 `q_value <= 0.01` 的 target 组汇总导出（tools.rs:3595 统计 `groups_at_1pct`）。
- I/L 归一在 mapper（建图）与 coverage（比对）两处各做一次，规则一致（`'I' -> 'L'`），保证肽 key 与 FASTA 比对同口径。
- `assign_razor_peptides` 返回的 `HashMap<pep, leader>` 经 `infer_proteins` 原样放入 `InferenceResult.razor_map`（core/protein_group.rs:53），供下游定量按 razor 归属避免重复计数。

## 6. 测试入口

```
cargo test -p protein-copilot-protein-inference --offline
```

57 个测试 = 43 单元（mapper 9 + parsimony 12 + razor 10 + coverage 12）+ 14 集成（tests/integration.rs，全链 mapper -> parsimony -> razor -> picked FDR -> coverage）。覆盖: 空输入 / 全诱饵两类 `InferenceError`、I/L 等价合并、共享肽 razor 三档 tie-break（unique / score / leader）、子集吸收与不可区分分组、贪心最小集确定性、覆盖率重叠 / 多次出现 / 去重 / 缺序列记 None。集成数据集（integration.rs:45-50）含不可区分对（P002 / P004 肽集相同合为一组）、跨蛋白共享肽（ACDEFGHK 由 P001 / P003 共享）与 I/L 等价（KDEFGHIJ vs KDEFGHLJ），覆盖真实推断歧义。

---

回到 [README](README.md) 选择其它层级或子系统。
