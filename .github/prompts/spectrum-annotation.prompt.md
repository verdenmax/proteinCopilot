---
mode: agent
description: "谱图标注与可视化 — 对单张谱图进行 b/y 碎片离子匹配，支持 SILAC 重标 Mirror Plot，生成交互式 HTML 标注图"
---

# 谱图标注

对单张质谱谱图进行肽段碎片离子匹配分析，生成包含 b/y 离子标注的交互式 HTML 可视化文件。
支持 SILAC 重标实验的 Mirror Plot（轻标朝上、重标朝下）和双窗口 XIC。

## 输入要求

**模式 1：从已有搜索结果标注**
- 搜索的 run_id（从 `run_search` 获得）
- 要标注的 scan number

**模式 2：手动指定肽段标注**
- 谱图文件路径（.mgf 或 .mzML）
- scan number
- 肽段序列（如 "PEPTIDEK"）
- 电荷态（如 2）

## 四种标注场景（决策矩阵）

根据**采集模式**（DDA/DIA）和**是否重标**（SILAC），标注分为 4 种情况：

| 场景 | 采集模式 | 重标 | `label_type` | 输出 |
|------|----------|------|-------------|------|
| ① | DDA | 无 | 不传 | 标准谱图 + 单通道 XIC |
| ② | DDA | SILAC | **必须传** | **Mirror Plot** + 双轨 XIC |
| ③ | DIA | 无 | 不传 | 标准谱图 + DIA cycle XIC |
| ④ | DIA | SILAC | **必须传** | **Mirror Plot** + 双窗口 XIC |

### 关键区别
- **场景 ②（DDA+SILAC）**：重标母离子被仪器作为独立前体选择，在**不同的 DDA scan** 中。工具自动在附近 scan 中按前体 m/z 匹配找到重标 scan。
- **场景 ④（DIA+SILAC）**：重标母离子 m/z 不同，落在**不同的 DIA 隔离窗口**。工具自动在附近 scan 中按隔离窗口包含关系找到重标 scan。
- DIA 自动检测：中位隔离窗口宽度 > 5 Da 判定为 DIA。仪器级别 DDA 窗口通常 < 2 Th，DIA 窗口通常 10-25 Da。

## SILAC / 重标检测（关键）

**必须在调用 `annotate_spectrum` 之前判断是否为 SILAC 实验。**
以下任一条件满足时，必须传递 `label_type` 参数：

1. **用户明确提到**：SILAC、重标、heavy label、轻重标、K+8、R+10
2. **搜索参数中包含 SILAC**：之前的 `run_search` 使用了 SILAC 相关修饰
3. **文件名线索**：文件名包含 `SILAC`、`Heavy`、`H/L` 等关键词
4. **用户要求 mirror plot** 或 "轻重标对比"

**SILAC 数据不传 `label_type` = 结果错误（不是可选项）：**
- **XIC 错误**：DIA 模式下重标母离子 m/z 与轻标不同，落在不同的隔离窗口。不传 `label_type`，XIC 只从轻标窗口取数据，重标通道的 MS2 色谱**完全丢失**，显示的色谱图是**错误的**
- **谱图标注不完整**：缺少重标匹配信息，无法进行轻重标对比验证
- **MS1 XIC 缺失**：不传 `label_type` 只提取轻标母离子 MS1 曲线，重标母离子曲线丢失

**有重标 = 必须使用 Mirror Plot：**
- SILAC 实验的标注**必须**生成 Mirror Plot（轻标朝上 + 重标朝下）
- 不存在"SILAC 数据但只看轻标"的正常用法 — 那是数据丢失

**标准 SILAC label_type 参数：**
```json
{
  "label_type": {
    "Silac": {
      "heavy_k_delta": 8.014199,
      "heavy_r_delta": 10.008269
    }
  }
}
```
> K+8 = ¹³C₆¹⁵N₂-Lysine, R+10 = ¹³C₆¹⁵N₄-Arginine（最常用的 SILAC 标记）

## 流程

1. 确认用户要标注的谱图（scan number）和标注模式
2. **判断是否为 SILAC 实验**（见上方检测规则），如不确定则询问用户
3. 如果是模式 1，从搜索结果中查找该 scan 对应的 PSM
4. 如果是模式 2，使用用户指定的肽段和电荷态
5. 调用 `annotate_spectrum`：
   - 非 SILAC：不传 `label_type`，生成标准谱图
   - **SILAC：必须传 `label_type`**，工具会自动：
     - 计算重标母离子 m/z（轻标理论值 + SILAC 质量偏移 / 电荷）
     - 在附近扫描中找到包含重标母离子的 DIA 窗口 MS2 谱图
     - 对重标谱图进行独立的 b/y 离子匹配
     - 生成 **Mirror Plot**（轻标蓝色朝上 + 重标橙色朝下）
     - XIC 从两个不同的 DIA 窗口提取轻标和重标色谱（实线=轻标，虚线=重标）
6. 生成交互式 HTML 文件：
   - 非 SILAC：标准谱图（b/y 离子标注 + 覆盖图 + XIC）
   - SILAC：**Mirror Plot + 双轨 XIC**（这是唯一正确的标注方式，不传 label_type 的结果是错误的）
7. 向用户报告匹配结果：
   - 轻标匹配分数
   - SILAC 时额外报告：重标匹配分数、重标 scan 号、轻重标 RT 差异
8. 告知用户在浏览器中打开 HTML 文件查看详细标注

## 输出
- 交互式 HTML 标注文件（自包含，无外部依赖）
- 匹配统计摘要（分数、离子匹配数、质量偏差）
- SILAC 时：轻标和重标双份统计
- 自然语言解读（匹配质量评估）

## 结果解读指导

- **score > 0.5**：匹配质量较好，大部分碎片离子被检测到
- **score 0.2-0.5**：部分匹配，可能存在噪声或修饰未考虑
- **score < 0.2**：匹配较差，可能是错误鉴定或参数不合适
- **|Δ ppm| < 5**：前体离子质量匹配精确
- **|Δ ppm| > 10**：质量偏差较大，检查仪器校准或序列是否正确

### SILAC 特有解读
- **轻重标 score 接近**：SILAC 标记成功，两个通道均有良好碎片化
- **重标 score 明显低于轻标**：可能重标富集不足或重标母离子信号弱
- **Mirror Plot 中 y 离子全匹配但 b 离子少**：正常，C 端含 K/R 的 y 离子在 SILAC 中更易检测
- **"No DIA window found for heavy precursor"**：重标母离子 m/z 超出 DIA 窗口范围，无法提取重标谱图

## 适用场景
- 验证搜索结果中某个 PSM 的可靠性
- 手动检查感兴趣的谱图与候选肽段的匹配
- **SILAC 实验**：验证轻重标肽段的共洗脱和碎片化质量
- 教学演示碎片离子匹配原理
- 论文图片准备（HTML 可截图，Mirror Plot 是 SILAC 论文常用图形）
