# L4 — fdr crate

承接 [L3-fdr-protein](L3-fdr-protein.md)，回溯 [L2](L2-architecture.md)。本篇只聚焦 `crates/fdr` 一个 crate 的函数级 API 与算法骨架，签名、常量、阈值均按源码核验，完整逻辑以源码为准。

## 1. 用途 + 位置 + 依赖

`protein-copilot-fdr` 提供纯确定性的 target-decoy FDR 估计与 q-value 赋值：建库阶段的诱饵生成、PSM/肽/蛋白三级 FDR、单调 q-value。无 LLM、无全局可变状态、无 I/O。

位置（上游 PSM，下游 report / protein-inference）:

```
search-engine: Vec<Psm> (score, is_decoy)
   |
   v
fdr::generate_decoys       -> REV_/SHUF_ 诱饵库（建库阶段）
fdr::calculate_fdr         -> PSM 级 q-value（SimpleSearch 内联调）
fdr::calculate_peptide_fdr -> 肽级 q-value（复用内核）
fdr::calculate_protein_fdr -> 蛋白 picked FDR（infer_proteins 调）
   |
   v
report / protein-inference: 按 1% FDR 汇总导出
```

依赖（`Cargo.toml`）: `protein-copilot-core`（Psm / ProteinGroup / DecoyStrategy）、`thiserror`、`serde`、`tracing`、`rand 0.8 (std)`。模块导出见 `lib.rs:10-22`: `calculation` / `decoy` / `peptide_fdr` / `protein_fdr` / `error`。

## 2. 函数级 API

| 函数 / 类型 | 签名（源码原样） | 位置 |
|---|---|---|
| `calculate_fdr` | `fn calculate_fdr(psms: &[ScoredPsm]) -> Result<Vec<(usize, f64)>, FdrError>` | calculation.rs:30 |
| `ScoredPsm` | `struct { index: usize, score: f64, is_decoy: bool }` | calculation.rs:13 |
| `calculate_peptide_fdr` | `fn calculate_peptide_fdr(peptides: &[PeptideScore]) -> Result<PeptideFdrResult, FdrError>` | peptide_fdr.rs:44 |
| `extract_unique_peptides` | `fn extract_unique_peptides(psms: &[core::search_result::Psm]) -> Vec<PeptideScore>` | peptide_fdr.rs:86 |
| `PeptideScore` | `struct { sequence: String, best_score: f64, is_decoy: bool }` | peptide_fdr.rs:26 |
| `PeptideFdrResult` | `struct { peptide_q_values: HashMap<String,f64>, target_peptides_at_1pct: u64, total_decoy_peptides: u64, total_peptides: u64 }` | peptide_fdr.rs:13 |
| `calculate_protein_fdr` | `fn calculate_protein_fdr(groups: &[ProteinGroup]) -> Result<ProteinFdrResult, FdrError>` | protein_fdr.rs:36 |
| `ProteinFdrResult` | `struct { groups: Vec<ProteinGroup>, target_groups_at_1pct: u64, total_target_groups: u64, total_decoy_groups: u64 }` | protein_fdr.rs:21 |
| `generate_decoys` | `fn generate_decoys(proteins: &[(String,String,String)], strategy: DecoyStrategy) -> Vec<DecoyProtein>` | decoy.rs:25 |
| `DecoyProtein` | `struct { accession: String, description: String, sequence: String }` | decoy.rs:10 |
| `FdrError` | `enum { NoPsms, NoDecoyHits, InvalidScore }` | error.rs:7 |

常量/默认阈值: 1% FDR 硬编码 `0.01`（calculation.rs:125、peptide_fdr.rs:67、protein_fdr.rs:155）; 诱饵前缀 `REV_`/`SHUF_`（decoy.rs:36,55）; picked 配对前缀 `DECOY_PREFIX = "REV_"`（protein_fdr.rs:17）; Shuffle 固定种子 `seed_from_u64(42)`（decoy.rs:45）; 进度日志步长 `5000`（calculation.rs:56）。

库代码不 `unwrap/expect`，三个错误变体经 `From<FdrError> for CoreError`（error.rs:21-35）统一落到 `CoreError::ValidationError { context, detail, suggestion }`，并附可执行建议: `NoPsms` 提示确认搜索是否产出 PSM、`NoDecoyHits` 提示改用 Reverse 策略、`InvalidScore` 提示检查打分。调用方据此优雅降级而非 panic。

## 3. 算法骨架

`calculate_fdr`（PSM/肽/蛋白三级共用内核）:

```
1 空输入 -> NoPsms; 任一 score 非有限 -> InvalidScore; 无 decoy -> NoDecoyHits
2 sort_by(b.score.total_cmp(&a.score))      # 降序, 对 NaN/-0 也成全序
3 逐位累加 targets/decoys, raw = (decoys/targets).min(1.0)
4 tie 拉平: 每段等分数区间 [lo,hi) 全部置为 raw[hi-1]
5 向后单调: q[i] = q[i].min(q[i+1])
6 zip 回原 index -> Vec<(index, q)>, 统计 q <= 0.01 的数量
```

第 4 步是关键: 同分 PSM 必须共享一个 q-value（取该分数段最后一位、计满全部 target+decoy 后的 FDR），消除输入顺序依赖。

`generate_decoys`（种子策略）: `Reverse` 反转序列但保留末位氨基酸（C 端 K/R，使诱饵肽与 target 切性质相近）; `Shuffle` 用固定种子洗牌 `chars[..last]`，两次生成一致; `None` 返回空 vec。

`calculate_protein_fdr`（picked target-decoy 配对）: 以去 `REV_` 前缀的 base accession 把 target 与 decoy 组配对; 组内竞争分高者胜、平局 target 胜; 未配对 target 自动胜出, 未配对 decoy 作为 decoy 胜者入列; 若无 decoy 胜者则全部 target 记 `q = 0.0`, 否则把胜者列喂回 `calculate_fdr`; 输出按 `score 降序 then leader_accession 升序`。

## 4. 简化源码片段

calculate_fdr 的 tie 拉平 + 向后单调（calculation.rs:98-117）:

```rust
let mut lo = 0;
while lo < sorted.len() {
    let mut hi = lo + 1;
    while hi < sorted.len()
        && sorted[hi].score.total_cmp(&sorted[lo].score) == Ordering::Equal {
        hi += 1;                          // 收拢一段等分数 [lo, hi)
    }
    let tie_fdr = raw_fdrs[hi - 1];       // 段末 FDR（已计满 target+decoy）
    for fdr in &mut raw_fdrs[lo..hi] { *fdr = tie_fdr; }
    lo = hi;
}
let mut q_values = raw_fdrs;
for i in (0..q_values.len().saturating_sub(1)).rev() {
    q_values[i] = q_values[i].min(q_values[i + 1]);   // 向后单调拉平
}
```

decoy reverse 保留末位（decoy.rs:69-79）:

```rust
fn reverse_sequence(seq: &str) -> String {
    let chars: Vec<char> = seq.chars().collect();
    if chars.len() <= 1 { return seq.to_string(); }
    let last = chars[chars.len() - 1];
    let mut middle = chars[..chars.len() - 1].to_vec();
    middle.reverse();
    middle.push(last);                    // "PEPTIDEK" -> "EDITPEPK"
    middle.into_iter().collect()
}
```

protein picked 组内配对竞争（protein_fdr.rs:66-86）:

```rust
for (&target_acc, &target_idx) in &target_map {
    if let Some(&decoy_idx) = decoy_map.get(target_acc) {
        paired_decoy_keys.insert(target_acc.to_string());
        if groups[target_idx].score >= groups[decoy_idx].score {
            winners.push((target_idx, false));   // 平局 target 胜
        } else {
            winners.push((decoy_idx, true));
        }
    } else {
        winners.push((target_idx, false));       // 未配对 target 自动胜
    }
}
for (base_acc, &decoy_idx) in &decoy_map {       // 未配对 decoy 入列为 decoy 胜者
    if !paired_decoy_keys.contains(base_acc) { winners.push((decoy_idx, true)); }
}
```

## 5. 调用链

- `search-engine/simple_engine.rs:121-144`: `decoy_strategy != None` 时把 `Vec<Psm>` 映射为 `ScoredPsm` -> `calculate_fdr`（call 在 131）-> 写回 `psms[idx].q_value`，随后 `retain(|p| !p.is_decoy)`（144）删诱饵; 计算失败仅 `warn` 不中断。
- `search-engine/simple_engine.rs:310,700`: 建库阶段 `generate_decoys(&target_tuples, params.decoy_strategy)` 产 `REV_`/`SHUF_` 库再消化。
- `mcp-server/tools.rs:3571`: `infer_proteins` 工具在 parsimony + razor 之后调 `calculate_protein_fdr(&groups)`，失败时降级为不带 q-value 的组（仅 `warn`）。
- `peptide_fdr` 为公共 API（`calculate_peptide_fdr` / `extract_unique_peptides`），供肽级汇总复用同一内核；三级 FDR 都落到同一个 `calculate_fdr`。

## 6. 测试入口

```
cargo test -p protein-copilot-fdr --offline
```

38 个单元测试（calculation 8 + decoy 8 + peptide_fdr 9 + protein_fdr 13），覆盖: 空输入 / NaN / 无诱饵三类 `FdrError`、tie 顺序无关性、q-value 单调与 `[0,1]` 边界、Reverse 保末位、Shuffle 种子可复现、picked 配对各分支（target 胜 / decoy 胜 / 未配对）与确定性排序。

---

回到 [README](README.md) 选择其它层级或子系统。
