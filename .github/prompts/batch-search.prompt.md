---
mode: agent
description: "批量搜索 — 多文件质谱数据的批量搜索策略"
---

# 批量搜索

处理包含多个质谱文件的实验，如多样品比较、技术重复、分级分离等场景。

## 使用方式

`run_search` 的 `input_files` 参数接受文件列表：

```
run_search(
  params={...},
  input_files=["sample1.mzML", "sample2.mzML", "sample3.mzML"]
)
```

或通过 `prepare_search`：

```
prepare_search(
  input_files=["sample1.mzML", "sample2.mzML", "sample3.mzML"],
  organism="human"
)
```

## 流程

### 1. 数据概览
- 调用 `read_spectra` 检查**第一个文件**的数据特征
- 确认所有文件来自相同实验条件（相同仪器、相同采集模式）

### 2. 参数推荐
- `prepare_search` 或 `recommend_params` 基于第一个文件推荐参数
- 同一批次的文件通常使用**相同的搜索参数**

### 3. 执行搜索
- 将所有文件一次性传给 `run_search`
- 搜索引擎内部合并处理所有谱图
- 返回单个 `run_id`，统一管理

### 4. 蛋白推断（推荐）
- 多文件搜索后，调用 `infer_proteins(run_id=xxx, fasta_path=xxx)` 聚合蛋白结果
- 跨文件的 PSM 被合并后再做蛋白推断，提高覆盖率

### 5. 结果导出
- `generate_summary(run_id)` 显示合并统计
- `export_results(run_id)` 导出包含所有文件结果的 TSV/JSON

## 分组策略

| 场景 | 策略 |
|------|------|
| 同一样品的技术重复 | 合并到一个 run_search 调用 |
| 分级分离（fractionation） | 合并到一个 run_search 调用 |
| 不同实验条件 | 分别搜索，各自独立 run_id |
| 不同物种 | 必须分别搜索（不同数据库） |
| DDA + DIA 混合 | 必须分别处理（不同工作流） |

## 性能预期

- **SimpleSearch**：~100 谱图/秒（单线程）
- **Sage**：~1000-5000 谱图/秒（多线程，取决于 CPU 核数）
- 10 个文件 × 10,000 谱图/文件 = 100,000 谱图 → Sage 约 20-100 秒

## 注意事项
- 所有文件必须存在且可读（搜索前同步校验）
- 参数推荐基于第一个文件的谱图特征
- 建议同一批次文件使用相同的仪器参数
- 跨文件的蛋白推断需要 `infer_proteins`，不会自动执行

## 适用场景
- 多样品比较蛋白质组学实验
- 分级分离样品（SCX、高 pH RP 等）
- 技术或生物学重复实验
- 时间序列实验（不同时间点的样品）
