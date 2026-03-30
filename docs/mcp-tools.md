# MCP Tools 参考

ProteinCopilot MCP Server 提供 8 个工具，通过 JSON-RPC over stdio 暴露给 LLM。

> **当前搜索引擎**：MVP 使用内置的 **SimpleSearchEngine**（基于 b/y 离子匹配的简化搜索），
> 后续将接入 pFind 作为生产级搜索引擎。SimpleSearch 足以验证完整流程，但搜索质量和性能不如专业引擎。

## 启动方式

```bash
cargo run --release -p protein-copilot-mcp-server
```

或通过 `.mcp.json` 配置自动发现。

---

## read_spectra

读取质谱文件，返回统计摘要。

**输入**：
```json
{
  "file_path": "/data/sample.mgf"
}
```

**输出**：`SpectrumSummary`
```json
{
  "file_path": "/data/sample.mgf",
  "format": "Mgf",
  "total_spectra": 10000,
  "ms1_count": 500,
  "ms2_count": 9500,
  "mz_range": [100.0, 2000.0],
  "rt_range_sec": [60.0, 3600.0],
  "precursor_charge_distribution": {"2": 5000, "3": 3000},
  "median_peaks_per_spectrum": 256
}
```

**使用场景**：分析数据特征，为参数推荐提供输入。

---

## get_spectrum

按 scan 号读取单张谱图。

**输入**：
```json
{
  "file_path": "/data/sample.mgf",
  "scan_number": 42
}
```

**输出**：`Spectrum`（含 mz_array、intensity_array、precursors）

**使用场景**：检查特定谱图的详细数据。

---

## recommend_params

基于谱图特征推荐搜索参数。

**输入**：
```json
{
  "summary": { ... },
  "hints": {
    "experiment_type": "phosphorylation",
    "instrument_type": "Orbitrap",
    "enzyme": "Trypsin"
  }
}
```

**输出**：`AiDecision<SearchParams>`
```json
{
  "decision": { "enzyme": "Trypsin", "precursor_tolerance": {"value": 10, "unit": "Ppm"}, ... },
  "confidence": 0.92,
  "explanation": "Based on m/z range...",
  "input_summary": "10000 spectra, ...",
  "alternatives": ["Open search", "TMT labeled"],
  "evidence": ["m/z range: 100-2000", ...]
}
```

**使用场景**：自动推荐参数，LLM 向用户展示推荐理由。

---

## list_presets

列出所有内置搜索参数预设。

**输入**：无

**输出**：`{ "presets": [...] }`（5 个预设：standard、phospho、tmt、silac、open）

**使用场景**：用户选择预设而非自定义参数。

---

## run_search

执行蛋白质数据库搜索。

**输入**：
```json
{
  "params": {
    "database_path": "/data/human.fasta",
    "enzyme": "Trypsin",
    "missed_cleavages": 2,
    "precursor_tolerance": {"value": 10, "unit": "Ppm"},
    "fragment_tolerance": {"value": 0.02, "unit": "Da"},
    "fixed_modifications": [{"name": "Carbamidomethyl", "mass_delta": 57.021464, "residues": ["C"], "position": "Anywhere"}],
    "variable_modifications": [],
    "decoy_strategy": "Reverse"
  },
  "input_files": ["/data/sample.mgf"]
}
```

**输出**：`SearchResult`（含 PSMs、peptides、proteins、summary、metadata）

**使用场景**：执行搜索，必须在用户确认参数后调用。

---

## check_engine

检查搜索引擎状态。

**输入**：无

**输出**：`{ "engine": {...}, "status": "Healthy" }`

**使用场景**：验证引擎可用性。

---

## generate_summary

从搜索结果生成 FDR 过滤后的统计摘要。

**输入**：
```json
{
  "result": { ... }
}
```

**输出**：`SearchResultSummary`（含 psms_at_1pct_fdr、identification_rate 等）

**使用场景**：搜索完成后，生成供 LLM 解读的统计数据。

---

## export_results

将搜索结果导出为文件。

**输入**：
```json
{
  "result": { ... },
  "output_dir": "./output"
}
```

**输出**：导出文件列表（psm.tsv、peptide.tsv、protein.tsv、result.json、run_metadata.json）

**使用场景**：用户请求保存结果。
