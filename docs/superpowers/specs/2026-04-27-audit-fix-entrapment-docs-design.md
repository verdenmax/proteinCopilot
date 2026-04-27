# 审计修复：Entrapment Tool 文档补充 + XIC 性能描述修正

## 问题

全项目审计发现两个问题：

1. **4 个 entrapment MCP tool 对 AI 不可见**：`classify_entrapment_hits`、`analyze_entrapment_stats`、`find_similar_targets`、`annotate_provenance` 已在 MCP Server 中注册（tools.rs L3583-3881），但未出现在任何 agent 或 prompt 文件中，AI 无法发现和调用。

2. **XIC 性能描述过时**：`proteomics-search.agent.md` L307-308 仍描述 "全文件流式扫描… 100-150s"，但 XIC 已通过 `extract_xic_unified()` 优化为索引定向读取 <1s。该描述会误导 AI 避免调用 XIC 相关功能。

## 修改范围

仅修改 `.github/agents/proteomics-search.agent.md`，共两处：

### 修改 1：补充 4 个 entrapment tool

- **frontmatter tools 列表**：在 `diagnose_search` 后追加 4 行
- **可用工具说明章节**：追加 4 个 tool 的功能描述
- **新增工作流章节**：在 XIC Extraction Workflow 之后新增 `### Entrapment Analysis Workflow`

工作流步骤：
1. 准备 YAML config（定义 target/trap 规则）
2. `classify_entrapment_hits(results_file, config_file, target_fasta)` → 分类 trap PSM 为 L0-L4
3. `analyze_entrapment_stats(classified_file)` → 统计分析
4. 可选深入：`find_similar_targets(peptide, target_fasta)` 查同源肽段
5. 可选深入：`annotate_provenance(file_path, scan_number, trap_sequence)` 碎片离子溯源

### 修改 2：修复 XIC 性能描述

- **删除** L307-308 旧文本
- **替换为**：XIC 已索引化优化（`extract_xic_unified()`），单次调用 <1s（含 7.5GB mzML），顺序调用即可，无需特殊超时

## 不修改

- 不删除废弃函数代码（`extract_xic()` / `extract_xic_with_raw()` 保留）
- 不修改其他 agent/prompt 文件
