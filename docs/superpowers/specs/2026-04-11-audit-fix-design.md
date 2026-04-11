# Audit Fix Design — 全项目边界问题与一致性修复

> **日期**: 2026-04-11
> **触发**: 5-agent 并行审计发现 2 CRITICAL + 9 HIGH + 11 MEDIUM 共 22 项问题
> **策略**: 分三批（A→B→C）逐步提交，每批独立可验证

---

## 问题分类与降级说明

### C1/C2 降级为 HIGH
`apply_fixed_mods()` 不检查 terminal context，但当前所有预设的 ProteinNTerm/CTerm 修饰均为 **variable mod**（varmod.rs 已正确过滤，FW-1 ✅）。只有用户自定义固定修饰为 ProteinNTerm 时才触发。仍需修复以防未来使用。

### H5 (RT 单位) 确认为当前行为正确
行业惯例（OpenMS、Mascot）：缺失 unitAccession 时默认秒。现代工具（msconvert）总是标注单位。加 warning 日志即可。

### H2-H4 (XSS) 风险评估
本地 HTML 文件，数据来自用户搜索结果，非 Web 服务。风险极低但加转义是好习惯。

---

## 批次 A — 防御性验证（8 项）

### A1: apply_fixed_mods 加 terminal context (C1)
**文件**: `search-engine/matching.rs`
**改动**: 函数签名加 `is_protein_nterm: bool, is_protein_cterm: bool`
```rust
fn apply_fixed_mods(sequence: &str, mods: &[Modification],
                     is_protein_nterm: bool, is_protein_cterm: bool) -> f64 {
    // ProteinNTerm: 仅 is_protein_nterm 时 apply
    // ProteinCTerm: 仅 is_protein_cterm 时 apply
    // AnyNTerm/AnyCTerm/Anywhere: 无条件 apply
}
```
**调用点更新**: matching.rs 两处 `apply_fixed_mods(...)` 调用传入 `peptide.is_protein_nterm/cterm`。

### A2: annotate.rs apply_fixed_mod_mass 同步修复 (C2)
**文件**: `search-engine/annotate.rs`
**改动**: 同 A1 逻辑。annotate_spectrum 不传 terminal context（外部调用），默认 false（保守：不应用 ProteinNTerm 固定修饰）。

### A3: charge > 0 验证 (H1)
**文件**: `xic/extract.rs`
**改动**: extract_xic() 和 extract_xic_with_raw() 入口加:
```rust
if charge <= 0 {
    return Err(XicError::InvalidPeptide { detail: format!("charge must be > 0, got {charge}") });
}
```

### A4: scan_number ≥ 1 验证 (H7)
**文件**: `mcp-server/tools.rs`
**改动**: annotate_spectrum handler 入口调用 `validate_scan_number()`（已存在）。

### A5: duplicate scan warning (H6)
**文件**: `spectrum-io/indexed_mgf.rs`
**改动**: insert 前检查 `offsets.contains_key()`，有则 `tracing::warn!`。

### A6: RT 缺失单位加 warning (H5)
**文件**: `spectrum-io/mzml.rs`
**改动**: `unitAccession` 为空时加 `tracing::warn!("MS:1000016 scan start time missing unitAccession; assuming seconds")`。保持默认秒行为不变。

### A7: database_path 验证前移 (H8)
**文件**: `mcp-server/tools.rs`
**改动**: run_search handler 中 `params.validate()` 移到 cache insert 之前。加文件存在检查。

### A8: extract_xic 错误消息补全 (H9→实际是 M8)
**文件**: `mcp-server/tools.rs`
**改动**: manual mode 错误消息列出全部 5 个必需字段。

---

## 批次 B — 安全加固（3 项）

### B1-B3: JSON→HTML `</script>` 转义 (H2-H4)
**文件**: `report/visualize.rs`, `report/unified_visualize.rs`, `report/xic_visualize.rs`
**改动**: 统一使用辅助函数:
```rust
fn escape_json_for_html(json: &str) -> String {
    json.replace("</script>", "<\\/script>")
        .replace("</Script>", "<\\/Script>")
}
```
在 `serde_json::to_string()` 后调用，然后再注入 HTML template。

---

## 批次 C — 边界处理与文档（11 项）

### C-M1: FDR NoDecoyHits 检查
**文件**: `fdr/calculation.rs`
**改动**: calculate_fdr 入口检查 `psms.iter().any(|p| p.is_decoy)`，无 decoy 返回 `FdrError::NoDecoyHits`。

### C-M2: DIA 无隔离窗 warning
**文件**: `dia-extraction/lib.rs`
**改动**: DIA 模式下 MS2 无 isolation_window 时 `tracing::warn!`。

### C-M3: NaN/Inf RT 验证
**文件**: `dia-extraction/correlation.rs`
**改动**: `find_associated_ms1` 跳过 RT 为 NaN/Inf 的 MS1。

### C-M4 + C-M5: K/R 大小写处理
**文件**: `xic/extract.rs`, `xic/heavy.rs`
**改动**: compute_ion_metadata 和 residue_heavy_delta 中 peptide 先 `.to_uppercase()`。

### C-M6: custom_json 位置验证
**文件**: `result-import/custom_json.rs`
**改动**: 验证 position ∈ [1, sequence.len()] 再调用 to_modification()。

### C-M7: run_search 文件检查
**文件**: `mcp-server/tools.rs`
**改动**: input_files 遍历检查 `Path::exists()`，不存在时返回明确错误。

### C-M8: same_isolation_window 容差可配
**文件**: `xic/extract.rs`
**改动**: 将 1.0 mz / 20% width 提取为 ExtractionParams 字段，默认值不变。

### C-M9: mass_delta 符号约定文档
**文件**: `core/search_params.rs`
**改动**: 加 doc comment "positive=addition, negative=loss"。

### C-M10: summary charge=0 过滤
**文件**: `report/summary.rs`
**改动**: 跳过 charge≤0 的 PSM。

### C-M11: target_scan 未找到时报错
**文件**: `xic/extract.rs`
**改动**: target_pos 为 None 时返回 XicError 而非使用全范围。

---

## 不修改项

- **H9 (same_isolation_window 硬编码容差)**: 归入 C-M8 作为 MEDIUM 处理
- 文档型 LOW 问题：随各批次修改一起补充 doc comment

---

## 测试策略

- 每批修复后 `cargo test --workspace` 全通过
- 新增测试覆盖:
  - A1: 测试 ProteinNTerm 固定修饰 + 非末端肽段 → 不应用
  - A3: 测试 charge=0 和 charge=-1 → 返回错误
  - B1: 测试包含 `</script>` 的 JSON → 转义后不破坏 HTML
  - C-M1: 测试全 target PSMs → NoDecoyHits 错误
