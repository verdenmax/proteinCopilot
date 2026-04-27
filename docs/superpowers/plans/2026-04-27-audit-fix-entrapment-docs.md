# 审计修复：Entrapment Tool 文档 + XIC 性能描述 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复全项目审计发现的 2 个文档问题：4 个 entrapment MCP tool 对 AI 不可见 + XIC 性能描述过时。

**Architecture:** 仅修改 `.github/agents/proteomics-search.agent.md` 一个文件，分三处编辑：frontmatter tools 列表、可用工具说明章节、XIC 性能描述段落。新增 Entrapment Analysis Workflow 章节。

**Tech Stack:** Markdown (agent 定义文件)

**Spec:** `docs/superpowers/specs/2026-04-27-audit-fix-entrapment-docs-design.md`

---

### Task 1: frontmatter tools 列表追加 4 个 entrapment tool

**Files:**
- Modify: `.github/agents/proteomics-search.agent.md` (lines 3-26, frontmatter tools 列表)

- [ ] **Step 1: 在 `diagnose_search` 后追加 4 行 tool 名称**

在 line 26 (`  - diagnose_search`) 之后、`---` 之前，追加：

```yaml
  - classify_entrapment_hits
  - analyze_entrapment_stats
  - find_similar_targets
  - annotate_provenance
```

修改后 frontmatter tools 列表变为（从 line 22 开始）：

```yaml
  - list_databases
  - download_database
  - get_database_info
  - diagnose_search
  - classify_entrapment_hits
  - analyze_entrapment_stats
  - find_similar_targets
  - annotate_provenance
---
```

- [ ] **Step 2: 验证 frontmatter 格式正确**

Run: `head -35 .github/agents/proteomics-search.agent.md`
Expected: YAML frontmatter 以 `---` 开始和结束，tools 列表包含 27 个 tool（原 23 个 + 4 个新增）。

---

### Task 2: 可用工具说明章节追加 4 个 tool 描述 + 新增 Entrapment Analysis Workflow

**Files:**
- Modify: `.github/agents/proteomics-search.agent.md` (lines 266-320 区域)

- [ ] **Step 1: 在 `## 可用工具说明` 章节尾部（line 283 `import_search_results` 描述之后）追加 4 个 tool 描述**

在 line 283 之后（`with annotate_spectrum, extract_xic, and generate_summary.` 这行之后），追加：

```markdown

- **classify_entrapment_hits**: 对 trap 数据库 PSM 进行同源性分类。读取搜索结果文件，
  应用 YAML config 中的 target/trap 规则，消化 target FASTA，将每个 trap PSM 分类为
  L0（完全匹配）到 L4（无同源目标）。可选 mzml_dir 参数启用碎片离子溯源。
  输出 classified.tsv、entrapment_report.html 等文件。

- **analyze_entrapment_stats**: 对 classify_entrapment_hits 的输出进行统计分析。
  返回 level 分布、蛋白家族聚类、delta-mass 分析。用于解释 entrapment 分类结果。

- **find_similar_targets**: 查找与给定肽段序列相似的 target 肽段。使用编辑距离
  （同长 Hamming、异长 Levenshtein）比对，返回最近匹配及质量差、替换类型标注。
  用于深入调查单个 trap PSM 的同源性。

- **annotate_provenance**: 单谱图碎片离子来源溯源。生成 mirror plot HTML，
  显示哪些峰来自 trap 肽段、target 肽段、两者共有（shared）或未匹配（unassigned）。
  用于可视化验证 trap PSM 的碎片离子归属。
```

- [ ] **Step 2: 在 `### Single Spectrum Inspection`（line 317）之后新增 `### Entrapment Analysis Workflow` 章节**

在 line 320（`3. Optionally use \`get_spectrum\` to see raw peak data`）之后，追加：

```markdown

### Entrapment Analysis Workflow

Entrapment 分析用于评估搜索引擎的假阳性控制质量：向数据库中添加已知不存在的 trap 蛋白，
检查搜索结果中有多少 trap PSM 通过了 FDR 阈值，并分析它们与真实 target 蛋白的同源性。

1. **准备 YAML config**：定义 target/trap 规则（哪些 accession 前缀是 target、哪些是 trap）
2. **分类**：`classify_entrapment_hits(results_file=xxx, config_file=yyy, target_fasta=zzz)`
   - 可选：`mzml_dir` 参数启用碎片离子溯源（判断 trap PSM 是否来自嵌合谱图）
   - 输出 classified.tsv（每个 trap PSM 标注 L0-L4 级别）和 entrapment_report.html
3. **统计分析**：`analyze_entrapment_stats(classified_file=output/entrapment/classified.tsv)`
   - 返回 level 分布、蛋白家族聚类、delta-mass 分析
4. **深入调查**（可选）：
   - `find_similar_targets(peptide="XXXK", target_fasta=zzz)` — 查找 trap 肽段的同源 target
   - `annotate_provenance(file_path=xxx.mzML, scan_number=1234, trap_sequence="XXXK")` — 碎片离子溯源可视化

**分类级别说明**：
- **L0**：trap 肽段在 target 蛋白组中有完全匹配（razor 肽段错误）
- **L1**：存在 ≤2 个氨基酸替换的 target 同源肽段（近同源假阳性）
- **L2**：存在同源但替换 >2 的 target 肽段
- **L3**：有 target 同源蛋白但无肽段级同源
- **L4**：无任何 target 同源（真正的假阳性）
```

---

### Task 3: 修复 XIC 性能描述

**Files:**
- Modify: `.github/agents/proteomics-search.agent.md` (lines 307-308)

- [ ] **Step 1: 替换 L307-308 的过时性能描述**

替换前（line 307-308）：

```markdown
**⚠️ 性能注意**：`extract_xic` 和 `annotate_spectrum` 内的 XIC 提取都需要全文件流式扫描。
大文件（>5GB）每次调用耗时 100-150s。`annotate_spectrum` 已包含 XIC，通常不需要额外调用 `extract_xic`。
```

替换后：

```markdown
**性能**：XIC 提取已索引优化（`extract_xic_unified`），单次调用 <1s（索引就绪后）。
首次打开新 mzML 文件需构建索引（7.5GB 约 10-30s），索引持久化为 `.mzML.idx` 后永久缓存。
`annotate_spectrum` 已包含 XIC，通常不需要额外调用 `extract_xic`。批量标注可顺序调用，无需特殊超时。
```

---

### Task 4: 验证 + Commit

- [ ] **Step 1: 检查文件格式无误**

Run: `head -35 .github/agents/proteomics-search.agent.md && echo "---" && grep -c "^  - " .github/agents/proteomics-search.agent.md`
Expected: frontmatter 正确闭合，tool 计数为 27

- [ ] **Step 2: 检查新增内容在文件中的位置**

Run: `grep -n "classify_entrapment\|analyze_entrapment\|find_similar_targets\|annotate_provenance\|Entrapment Analysis\|extract_xic_unified" .github/agents/proteomics-search.agent.md`
Expected: 6+ 行匹配，覆盖 frontmatter、工具说明、工作流、性能描述

- [ ] **Step 3: Commit**

```bash
git add .github/agents/proteomics-search.agent.md
git commit -m "docs(agent): add 4 entrapment tools + fix stale XIC perf description

- Add classify_entrapment_hits, analyze_entrapment_stats,
  find_similar_targets, annotate_provenance to tools list
- Add tool descriptions and Entrapment Analysis Workflow section
- Fix XIC performance: '全文件流式扫描 100-150s' → '索引优化 <1s'

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
