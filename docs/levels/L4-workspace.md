# L4 — workspace 工程约定（收尾篇）

返回索引 [README](README.md)，架构全局见 [L2](L2-architecture.md)。本篇不讲单个 crate，而是把 15 个 crate 串起来的 workspace 层约定：构建、依赖、测试、规范、扩展。所有事实核对自根 `Cargo.toml`、`.clippy.toml`、`rustfmt.toml`、`.github/`，标注来源，不臆造。

## 1. workspace 总览

根 `Cargo.toml`（行 1-9）只做两件事：声明成员、固定公共元数据。

```
[workspace]
resolver = "2"                 # 行 2，新依赖解析器
members  = ["crates/*"]        # 行 3，glob 自动纳入 crates/ 下全部目录

[workspace.package]            # 行 5-9，单一真相源
version      = "0.1.0"
edition      = "2021"
rust-version = "1.85"
license      = "MIT"
```

`resolver = "2"` 启用 Rust 2021 的特性解析，避免 dev/build 依赖的 feature 泄漏进正式构建；`members` 用 glob 而非逐个列举，新增目录即自动成为成员，省去维护成本。`cargo metadata --offline` 实测 15 个成员，命名一律 `protein-copilot-*`。其中 core 与 fasta-db 不依赖任何内部 crate，是依赖图的两个根；mcp-server 居于顶端，组装全部库。角色一句话：

| crate | 类型 | 角色 |
|-------|------|------|
| core | lib | 共享领域类型 + SearchEngineAdapter trait，依赖图根 |
| fasta-db | lib | FASTA 库注册/下载/缓存，另一依赖图根 |
| spectrum-io | lib | mzML/mgf 解析 + ScanIndex 随机访问 |
| param-recommend | lib | 确定性搜索参数推荐 |
| search-engine | lib | 引擎调度 + SimpleSearch + Sage/pFind adapter |
| fdr | lib | target-decoy FDR + decoy 库生成 |
| protein-inference | lib | parsimony 蛋白推断 + 蛋白级 FDR |
| xic | lib | XIC 色谱提取（DDA/DIA + SILAC）+ Plotly |
| dia-extraction | lib | DIA 前体离子提取（同位素模式） |
| result-import | lib | 外部结果导入（DIA-NN/pFind/JSON） |
| report | lib | 统计摘要 + TSV/JSON 导出 |
| entrapment-analysis | lib | entrapment 命中分类 + 同源分析 |
| entrapment-cli | bin | entrapment 分类命令行（bin `entrapment`） |
| mcp-server | bin | MCP Server 组装（bin `protein-copilot-mcp`，27 工具） |
| integration-tests | lib* | 跨 crate 集成测试 harness（`test_helpers`，publish=false） |

合计 12 纯库 + 2 bin（mcp-server、entrapment-cli）+ 1 测试 harness。`edition = "2021"` 统一语言版次，`rust-version = "1.85"` 声明最低工具链（MSRV），低于此版本的 cargo 会直接拒绝构建。

## 2. 共享依赖

版本只在根定义两处，子 crate 引用而不重复写版本号：

- `[workspace.package]` -> 元数据（上一节）；子 crate 写 `version.workspace = true` 继承。
- `[workspace.dependencies]`（行 11-72）-> 依赖版本；子 crate 写 `dep = { workspace = true }`。

外部依赖按域分组（版本来自根，全 workspace 唯一）：

| 域 | 依赖 |
|----|------|
| 序列化 | serde 1 (+derive), serde_json 1, schemars 1 (+uuid1,chrono04) |
| 错误 | thiserror 1, anyhow 1 |
| 标识/时间 | uuid 1.10 (+v4,serde), chrono 0.4 (+serde) |
| 异步 | tokio 1 (full), async-trait 0.1 |
| 可观测 | tracing 0.1, tracing-subscriber 0.3 (+env-filter,json) |
| 解析 | quick-xml 0.37, base64 0.22, flate2 1, memchr 2 |
| 网络/哈希 | reqwest 0.12 (rustls-tls+stream, default off), sha2 0.10 |
| 配置/CLI | serde_yaml 0.9, csv 1, clap 4 (+derive) |

内部 path 依赖 13 个（散落于行 40-62，与 reqwest/sha2 等外部依赖交错）：core、spectrum-io、param-recommend、search-engine、report、fdr、dia-extraction、xic、result-import、protein-inference、fasta-db、entrapment-analysis、entrapment-cli。子 crate 同样以 `{ workspace = true }` 引用。

这样做的好处是版本只有一处：升级 serde 只改根里一行，15 个 crate 同步生效，杜绝 workspace 内的版本漂移；`Cargo.lock` 进一步把整棵依赖树（326 包）锁死，保证换机器也能复现构建。

特例（crate-local，不进 workspace）：search-engine 的 `sage-core`（git pin rev cd712d4）与 `rayon`、fdr 的 `rand`、各 crate 的 `tempfile`（dev）。即"通用版本上提，专用版本就地"，避免把只有一个 crate 用的依赖塞进公共清单。

## 3. 构建/测试/质量命令

本环境无网络，所有 cargo 命令必须带 `--offline`（联网取 registry/git 会直接挂；Cargo.lock 已锁 326 包，离线即可命中缓存）。

```
cargo build  --workspace --offline
cargo test   --workspace --offline
cargo clippy --workspace --offline
cargo fmt    --all                          # 格式化无需网络
cargo test -p protein-copilot-fdr --offline # 单 crate
```

- `rustfmt.toml`：`max_width = 100`、`use_field_init_shorthand = true`。
- `.clippy.toml`：`avoid-breaking-exported-api = false`（允许 clippy 提示破坏导出 API 的改法）。
- warnings 当 error：写在 `copilot-instructions.md` §6（"CI 中视为错误"）与 `rust.instructions.md`。注意仓库当前无 `.github/workflows`、根也无 `[workspace.lints]`，故这是提交前人工跑 clippy 守门的约定，而非编译期 `deny` 强制。

clippy 与 fmt 是合并前的硬门槛：任一报错都视为未完成。声称功能完成前，必须贴出 `cargo test --workspace --offline` 全绿与 `cargo clippy --workspace --offline` 零警告的实际输出作为证据。

## 4. 工程规范要点（copilot-instructions.md）

下面几条贯穿全部 15 个 crate，是 review 的硬性检查项，新代码逐条对照：

- 确定性 / LLM 分层（§2.1）：Rust 只做确定性计算，LLM 只做推理解释；禁止 Rust 内调 LLM、禁止把 FDR/打分交给 LLM。
- 库代码无 `unwrap()` / `expect()`（§5.1）；用 `Result` + `?`。
- thiserror per crate（§5.2）：每个 crate 定义自己的 `XxxError` 枚举。
- serde everywhere（§2.5）：所有结果类型 `Serialize + Deserialize`。
- tracing（§2.8）：搜索/FDR/推参等关键操作必须有结构化日志。
- 确定可复现（§2.5）：每次运行生成唯一 `run_id`，落 manifest。
- Adapter 抽象（§2.6）：搜索引擎经 `SearchEngineAdapter` trait 接入，pFind/Sage/Comet 隔离在各自 adapter。
- MCP I/O 结构化（§2.3、§6）：工具输入输出必须 JSON Schema 可描述，禁自由文本；禁全局可变状态，用依赖注入。

## 5. 如何新增一个 crate

1. 建 `crates/<name>/src/`，写 `Cargo.toml`，`name = "protein-copilot-<name>"`，四项元数据写 `*.workspace = true`。`members = ["crates/*"]` 会自动纳入，无需改根 members。
2. 若该 crate 要被别的 crate 依赖：在根 `[workspace.dependencies]` 加一行 `protein-copilot-<name> = { path = "crates/<name>" }`。
3. `lib.rs` 导出清晰公共 API；新增 `error.rs`，用 thiserror 定 `<Name>Error`。
4. 依赖一律 `{ workspace = true }`；仅本 crate 用的小依赖才 crate-local。
5. TDD：先写失败测试，再 `cargo test -p protein-copilot-<name> --offline`、`cargo clippy -p protein-copilot-<name> --offline` 零警告。

由于成员是 glob，第 1 步建好目录后该 crate 即纳入 workspace 构建；第 2 步只在它需要被别的 crate 复用时才登记。

```
[package]
name = "protein-copilot-<name>"
version.workspace = true
edition.workspace = true
rust-version.workspace = true
license.workspace = true

[dependencies]
protein-copilot-core = { workspace = true }  # 内部依赖
thiserror = { workspace = true }             # 外部统一版本
```

## 6. 目录导航

想深入某一块，按下面的目录定位文档与配置（代码看 crates/，规范看 .github/，设计与计划看 docs/）：

```
docs/
  levels/        L1-L4 本套分层文档（README.md 为索引）
  superpowers/
    specs/       设计文档 YYYY-MM-DD-<topic>-design.md
    plans/       实施计划
  architecture/  crate-dependencies.svg + workflow-*.svg + gen_*.py
.github/
  copilot-instructions.md   工程总规范（9 节）
  instructions/             rust.instructions.md
  agents/                   7 个领域 Agent（architect/coder/...）
  prompts/                  Skill / Prompt 模板
crates/          15 个成员（见第 1 节表）
examples/        hela-mix-2da-entrapment.yaml, test_iso_error.rs
tests/fixtures/  small_test.fasta/.mgf, pfind_sample.tsv
```

至此 L1-L4 收尾。返回目录 [README](README.md)。
