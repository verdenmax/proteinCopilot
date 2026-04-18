---
mode: agent
description: "搜索诊断 — 分析搜索失败原因、评估结果质量、提供参数调优建议和重试方案"
---

# 搜索诊断

分析搜索运行的失败原因或结果质量异常，提供修复建议。

## 使用时机

- 搜索失败时（`get_search_status` 返回 status = "Failed..."）
- 搜索成功但结果异常时（鉴定率低、蛋白数量少等）
- 用户主动要求分析搜索质量

## 诊断流程

### 1. 获取诊断数据
- 调用 `diagnose_search(run_id=xxx)` 获取 `SearchDiagnostics`
- 前提：搜索已结束（`get_search_status` 的 `has_diagnostics = true`）

### 2. 搜索失败诊断

按以下顺序解读：

1. **error_category** → 确定大方向
   - `InputData`：文件损坏、格式错误、无 MS2 谱图
   - `Parameters`：容差不合理、酶不匹配
   - `Database`：物种不匹配、FASTA 格式错误
   - `Engine`：引擎内部错误、资源不足
2. **failure_stage** → 定位失败发生在哪个阶段
3. **stages[]** → 展示搜索"走了多远"（哪些阶段完成了）
4. **suggestions[]** → 按 priority 排序展示修复方案
5. 如果 suggestions 中有 **param_changes** → 直接向用户展示修改后的参数

### 3. 结果质量诊断

搜索成功但可能有质量问题：

1. **anomalies[]** → 列出检测到的异常
2. 对每个异常解释：
   - `LowIdentificationRate`：PSM 鉴定率过低，最常见原因是物种/酶不匹配
   - `NoDecoyHits`：FDR 无法计算，需检查数据库 decoy 序列
   - `HighFdr`：通过 FDR 的 PSM 太少，结果统计不可靠
   - `NarrowTolerance`：容差过窄可能排除正确匹配
   - `WideTolerance`：容差过宽导致假阳性多
   - `DatabaseMismatch`：物种不匹配嫌疑
   - `SlowSearch`：匹配阶段瓶颈
   - `LowSpectraQuality`：碎片谱质量不足
3. **suggestions[]** → 展示优化建议

### 4. 搜索正常

如果无异常：
- 简要展示各阶段耗时
- 确认结果质量正常
- 建议下一步：`infer_proteins` 或 `export_results`

## 重试搜索

如果 suggestions 包含 param_changes：
1. 向用户展示原始参数和建议调整
2. **等待用户确认**后才能使用新参数调用 `run_search`
3. 保留原始 run_id 供结果对比

## 领域参考值

| 指标 | 正常范围 | 异常阈值 |
|------|---------|---------|
| DDA 鉴定率（HeLa） | 15-40% | < 10% |
| PSM @ 1% FDR | > 1000 | < 50 |
| 搜索耗时（Sage） | 1-5 min | > 10 min |
| 搜索耗时（SimpleSearch） | 5-30 min | > 60 min |
| 前体容差 | 5-20 ppm | < 5 或 > 50 ppm |
