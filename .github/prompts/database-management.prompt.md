---
mode: agent
description: "FASTA 数据库管理 — 查看、下载和管理蛋白质序列数据库"
---

# 数据库管理

管理蛋白质序列数据库（FASTA 格式），支持内置数据库自动下载和本地缓存。

## 内置数据库

| ID | 物种 | 来源 | 关键词（自动匹配） |
|----|------|------|---------------------|
| human_swissprot | 人 (Homo sapiens) | UniProt Swiss-Prot | human, 人, 人类, homo sapiens, 9606 |
| mouse_swissprot | 小鼠 (Mus musculus) | UniProt Swiss-Prot | mouse, 小鼠, mus musculus, 10090 |
| ecoli_swissprot | 大肠杆菌 (E. coli) | UniProt Swiss-Prot | ecoli, e.coli, 大肠杆菌, escherichia |
| yeast_swissprot | 酵母 (S. cerevisiae) | UniProt Swiss-Prot | yeast, 酵母, saccharomyces |
| arabidopsis_swissprot | 拟南芥 (A. thaliana) | UniProt Swiss-Prot | arabidopsis, 拟南芥 |
| crap | 污染物 | cRAP | contaminant, 污染, crap |

## 流程

### 查看数据库状态
1. 调用 `list_databases()` 查看所有数据库
2. 每个数据库显示：Available（可下载）或 Downloaded（已缓存，含文件大小和蛋白数量）

### 下载数据库
1. 调用 `download_database(database_id="human_swissprot")`
2. 从 UniProt 通过 HTTPS 下载，自动解析 FASTA 统计蛋白数量
3. 缓存到本地目录，下次使用无需重新下载
4. 返回本地文件路径，可直接用作搜索的 `database_path`
5. 如需更新：`download_database(database_id="human_swissprot", force=true)`

### 查看数据库详情
- 调用 `get_database_info(database_id="human_swissprot")` 查看：
  - 蛋白序列数量
  - 文件大小
  - SHA256 校验和
  - 下载时间
  - 前 5 个蛋白 accession（验证正确性）

### 自动数据库解析（推荐）
- 使用 `prepare_search(organism="human")` 时自动处理：
  1. 检查本地缓存
  2. 未缓存则自动下载
  3. 填充到搜索参数的 database_path
- 支持中英文物种名和 NCBI Taxonomy ID

## cRAP 污染物数据库

- **Common Repository of Adventitious Proteins**
- 包含实验室常见污染物：人角蛋白、胰蛋白酶自切产物、BSA 等
- 建议在搜索时将 cRAP 序列合并到物种数据库中
- 高占比的 cRAP 鉴定可能提示样品制备问题

## 适用场景
- 首次使用时下载所需数据库
- 检查已缓存数据库是否需要更新
- 使用自定义 FASTA 数据库时，直接提供 database_path 即可
