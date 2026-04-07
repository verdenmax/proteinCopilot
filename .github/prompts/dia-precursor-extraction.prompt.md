---
mode: agent
description: "DIA 母离子提取 — 从 DIA 数据中提取候选母离子，支持全文件批量提取和单谱图分析"
---

# DIA 母离子提取

从 DIA（数据非依赖采集）质谱数据中提取候选母离子。DIA 谱图的隔离窗口较宽（通常 >5 Da），
包含多个共碎裂母离子，需要通过 MS1 同位素模式分析来还原候选母离子。

## 使用场景

### 场景 1：全文件 DIA 提取（搜索前准备）

用户提供 mzML 文件，需要提取所有 MS2 谱图的候选母离子用于后续数据库搜索。

**流程：**
1. 调用 `read_spectra` 获取文件摘要，确认数据特征
2. 调用 `recommend_params` 确认是否为 DIA 数据（会提示 DIA 检测结果）
3. 调用 `extract_dia_precursors` 进行批量母离子提取
   - 自动检测 DIA/DDA 模式
   - 可配置 output_mode: "pseudo"（每个母离子一张谱图）或 "multi"（保留多母离子）
   - 可配置电荷范围 (min_charge, max_charge)
4. 检查提取统计：提取率、电荷分布
5. 使用返回的 `dia_run_id` 传递给 `run_search` 执行搜索

### 场景 2：单谱图母离子分析

用户指定某个 scan number，查看该谱图的母离子提取详情。

**流程：**
1. 调用 `extract_spectrum_precursors` 指定文件和 scan number
2. 查看结果：
   - 使用了哪个 MS1 谱图（ms1_scan_used）
   - 关联方法（correlation_method）：source_scan / scan_order / rt_nearest
   - 提取到的候选母离子列表（m/z, charge, intensity）
   - 隔离窗口信息
3. 可以调用 `get_spectrum` 查看原始 MS1 和 MS2 谱图数据
4. 可以对提取到的候选进行手动标注验证

## 关键参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| min_charge | 2 | 最小电荷态 |
| max_charge | 5 | 最大电荷态 |
| output_mode | "pseudo" | "pseudo"=拆分, "multi"=保留多母离子 |
| acquisition_mode | auto | 可手动覆盖为 "DDA" 或 "DIA" |

## 判断标准

- **DIA 检测阈值**：中位隔离窗口宽度 > 5 Da
- **同位素模式匹配**：
  - 同位素间距 = 1.00335 / charge
  - 同位素容差 = 0.01 Da
  - 至少 2 个同位素峰
  - 强度递减趋势
- **MS1 关联策略**（三级回退）：
  1. spectrumRef 直接引用
  2. 扫描顺序（向前最近 MS1）
  3. 保留时间最接近

## 注意事项

- DDA 数据会被自动检测并跳过提取（直接返回原始谱图）
- 提取结果缓存在内存中，通过 `dia_run_id` 关联
- 建议先用 `read_spectra` 确认文件内容，再决定提取策略
- 如果提取率偏低，可以调整 min_charge/max_charge 范围
