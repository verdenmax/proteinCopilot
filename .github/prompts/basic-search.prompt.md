---
mode: agent
description: "基础蛋白质质谱搜索流程 — 从谱图文件到搜索结果的 step-by-step 引导"
---

# 基础搜索流程

执行一次标准的蛋白质质谱数据库搜索。

## 输入要求
- 谱图文件路径（.mgf 或 .mzML 格式）
- FASTA 蛋白数据库路径（如 UniProt Human）— 可选，支持自动下载内置数据库

## 流程

1. 请提供谱图文件路径和 FASTA 数据库路径
2. 我会先读取谱图文件，分析数据特征
3. **数据库选择**（如未指定 FASTA 路径）：
   - 从用户消息中检测物种意图：
     - "人"/"human"/"人类" → `human_swissprot`
     - "小鼠"/"mouse" → `mouse_swissprot`
     - "大肠杆菌"/"E.coli" → `ecoli_swissprot`
     - "酵母"/"yeast" → `yeast_swissprot`
     - "拟南芥"/"Arabidopsis" → `arabidopsis_swissprot`
   - 调用 `list_databases` 检查本地缓存状态
   - 未缓存则建议下载："检测到您需要搜索 [物种] 蛋白质组，需要下载 UniProt Swiss-Prot 数据库。是否下载？"
   - 已缓存则直接使用本地路径作为 `database_path`
   - 提示 cRAP 污染物数据库的可用性
4. 基于数据特征推荐搜索参数（酶、容差、修饰）
5. 向你展示推荐参数和理由，等待确认
6. 确认后执行搜索
7. 生成结果摘要并解读

## 输出
- 搜索结果统计：PSM/肽段/蛋白质鉴定数量
- 鉴定率和质量指标
- 自然语言结果解读
- 可选：导出 TSV/JSON 结果文件

## 适用场景
- 标准蛋白质组学搜索
- 首次分析新数据集
- 不确定该用什么参数时
