---
mode: agent
description: "蛋白推断工作流 — 从 PSM 搜索结果到蛋白质水平鉴定，支持 Parsimony 最小集和 Razor 肽段"
---

# 蛋白推断

将肽段-谱图匹配（PSM）结果聚合到蛋白质水平，确定样品中存在哪些蛋白质。

## 输入要求
- 已完成搜索的 run_id（从 `run_search` 获得，status = "Completed"）
- FASTA 数据库路径（与搜索使用的相同）

## 流程

1. 确认搜索已完成：`get_search_status(run_id)` → status = "Completed"
2. 调用 `infer_proteins(run_id=xxx, fasta_path=xxx)`
   - 可选参数：
     - `fdr_threshold`：FDR 阈值，默认 0.01（1%）
     - `min_peptides`：最少肽段数，默认 1
3. 解读推断结果：
   - **蛋白组数**：通过 FDR 阈值的蛋白组数量
   - **Unique peptides**：仅属于一个蛋白的肽段，是蛋白鉴定的最强证据
   - **Shared peptides**：被多个蛋白共享的肽段
   - **Razor peptides**：共享肽段中，归属到证据最多蛋白的那一份
   - **序列覆盖率**：匹配肽段覆盖蛋白序列的百分比

## 算法说明

### Parsimony（最小蛋白集）
- 找到能解释所有鉴定肽段的**最少蛋白质数量**
- 一个肽段可能匹配多个蛋白（同源蛋白、亚型）
- Parsimony 消除冗余：如果蛋白 A 的所有肽段都被蛋白 B 包含，则 A 是冗余的

### Razor 肽段分配
- 共享肽段归属到拥有最多 unique peptides 的蛋白质
- 每个共享肽段只计算一次（不重复计数）
- 这是 MaxQuant / Proteome Discoverer 使用的标准方法

### 多级 FDR
- **PSM FDR 1%** → **肽段 FDR 1%** → **蛋白 FDR 1%**
- 每一级独立过滤，逐级收紧
- 发表级别分析通常使用更严格的蛋白 FDR（如 0.1%）

## 结果质量评估

| 指标 | 参考范围（HeLa 标准样品） |
|------|--------------------------|
| 蛋白组数 | 3,000 - 6,000+ |
| Unique peptides/蛋白 中位数 | ≥ 2 |
| 序列覆盖率中位数 | 15-30% |
| 1-peptide 蛋白占比 | < 30%（过高提示数据深度不足） |

## 适用场景
- 标准蛋白质组学实验的蛋白水平报告
- 比较不同样品的蛋白鉴定差异
- 验证目标蛋白是否被鉴定到
