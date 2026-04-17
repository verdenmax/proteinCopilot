---
mode: agent
description: "Sage 搜索引擎 — 使用 sage-core 库进行高性能蛋白质数据库搜索"
---

# Sage 搜索引擎

使用 sage-core（v0.15.0）进行生产级蛋白质数据库搜索。Sage 是一个高性能开源搜索引擎，
内置于 ProteinCopilot 作为库调用（非子进程）。

## 使用方式

在搜索参数中指定引擎：

**方式 1：通过 prepare_search**
```
prepare_search(input_files=[...], organism="human", engine="Sage")
```

**方式 2：通过 run_search**
```
run_search(params={...engine: "Sage"...}, input_files=[...])
```

## Sage vs SimpleSearch 对比

| 特性 | Sage | SimpleSearch |
|------|------|--------------|
| 打分算法 | Hyperscore + LDA rescoring | 简化打分 |
| 并行化 | rayon 多线程 | 单线程 |
| FDR 计算 | 三级（spectrum/peptide/protein） | 基础 target-decoy |
| 适用场景 | 生产分析 | 快速测试、小规模数据 |
| 数据规模 | 大规模（万级谱图） | 小规模（千级谱图） |

## Sage 搜索流程

1. **参数转换**：ProteinCopilot SearchParams → Sage Parameters
   - 消化酶映射（Trypsin → KR|P 正则规则）
   - 修饰映射（名称 → 质量偏移）
   - 质量容差映射（Da/ppm → SageTolerance）
2. **数据库构建**：FASTA → IndexedDatabase（内存中构建 + 索引）
3. **打分**：rayon 并行谱图匹配，计算 hyperscore
4. **LDA Rescoring**：线性判别分析重打分，综合多维特征
5. **FDR 计算**：
   - Spectrum-level q-value（基于 discriminant_score）
   - Peptide-level q-value（最佳 PSM 代表）
   - Protein-level q-value（picked-protein 方法）

## 结果字段说明

搜索结果的 `extra` 字段包含 Sage 特有信息：

| 字段 | 说明 |
|------|------|
| `hyperscore` | Sage 原始打分（越高越好） |
| `discriminant_score` | LDA 重打分后的综合分数（用于排序） |
| `spectrum_q` | 谱图水平 q-value |
| `peptide_q` | 肽段水平 q-value |
| `protein_q` | 蛋白水平 q-value |
| `delta_hyperscore` | 最佳与次佳 PSM 的 hyperscore 差值 |
| `matched_intensity_pct` | 匹配碎片离子强度占总强度的百分比 |
| `poisson` | 随机匹配概率（Poisson 模型） |

## Sage 特有参数（内部默认值）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| 肽段质量范围 | 500-5000 Da | 过滤过短/过长肽段 |
| min_ion_index | 2 | 跳过前 2 个碎片离子（通常噪声高） |
| max_variable_mods | 2 | 每条肽段最多 2 个可变修饰 |
| chimera | false | 不启用嵌合谱图处理 |

这些参数在 SageAdapter 中硬编码为合理默认值，无需用户调整。

## 引擎健康检查

调用 `check_engine(engine="Sage")` 返回引擎名称、版本、健康状态，以及所有已注册引擎列表。

## 适用场景
- 生产级蛋白质组学搜索
- 大规模数据集（>10,000 谱图）
- 需要多级 FDR 控制的发表级分析
- 需要高性能并行处理时
