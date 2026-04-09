# mzML Scan Index & Spectrum Cache — Design Spec

> **日期**: 2026-04-09
> **状态**: Approved
> **分支**: `feat/scan-index-cache`

## 1. 问题

当前 `spectrum-io` 的 `read_spectrum(path, scan)` 每次调用都重新打开文件、从头线性扫描到目标 scan。没有任何索引或缓存。对于 1GB mzML（~50k scans），读取 scan=25000 需要 2-5 秒。标注 10 个 PSM 需要 20-50 秒。

## 2. 方案

### 2.1 两层缓存

| 层 | 缓存对象 | 生命周期 | 作用 |
|---|---|---|---|
| L1: ScanIndex | `HashMap<u32, u64>` scan→byte_offset | 文件级 | O(1) 定位 |
| L2: Reader Cache | `LruCache<PathBuf, Arc<IndexedMzMLReader>>` | MCP server 进程级 | 避免重复建索引 |

### 2.2 索引来源

1. **原生 indexList** (`<indexedmzML>` 格式): 读取文件末尾 `<indexListOffset>` → seek 到 `<indexList>` → 解析所有 offset。耗时 <5ms。
2. **全扫构建** (无索引 `<mzML>` 格式): 首次打开时线性扫描，记录每个 `<spectrum>` 标签的字节偏移。一次性成本。

### 2.3 新增类型

```rust
pub struct ScanIndex {
    offsets: HashMap<u32, u64>,
    source: IndexSource,
    total_spectra: u32,
}

pub enum IndexSource {
    NativeIndex,    // 来自 mzML <indexList>
    BuiltFromScan,  // 首次全扫描构建
}

pub struct IndexedMzMLReader {
    index: ScanIndex,
    path: PathBuf,
}
```

### 2.4 设计决策

- **索引构建在 spectrum-io crate 内部** — 封装性好，所有调用者自动受益
- **IndexedMzMLReader 是新类型** — 不改动现有无状态 MzMLReader
- **实现 SpectrumReader trait** — 完全兼容现有 API
- **MGF 也做行偏移索引** — IndexedMgfReader
- **MCP Server 层 LRU 缓存** — reader_cache 字段缓存 IndexedReader 实例
- **不做文件变更检测** — MVP 不需要，文件在分析期间不会被修改

## 3. 性能预期

| 指标 | Before | After | 改善 |
|---|---|---|---|
| read_spectrum (1GB, scan=25000) | 2-5s | <10ms | 200-500× |
| 标注 10 个 PSM | 20-50s | <100ms | 200-500× |
| 首次打开（有原生索引） | — | <5ms | 读末尾几KB |
| 首次打开（无索引） | — | 同 read_all 时间 | 一次性成本 |
| 内存占用 | 0 | ~2MB/文件 | 索引 HashMap |

## 4. 影响范围

### 修改
- `crates/spectrum-io/src/mzml.rs` — 新增索引解析和 seek 读取
- `crates/spectrum-io/src/mgf.rs` — 新增行偏移索引
- `crates/spectrum-io/src/lib.rs` — 导出新类型 + create_indexed_reader()
- `crates/spectrum-io/src/reader.rs` — 可选: trait 新增方法
- `crates/mcp-server/src/tools.rs` — 添加 reader_cache, 迁移工具

### 新增
- `crates/spectrum-io/src/index.rs` — ScanIndex, IndexSource 类型
- 测试 fixture: indexed mzML 文件

### 不修改
- `crates/core/` — 无变化
- `crates/xic/` — 后续可选迁移，本次不动
- 现有 SpectrumReader trait 签名 — 保持向后兼容
