# L4 — entrapment-cli crate

承接 [L3-entrapment](L3-entrapment.md)，回溯 [L2](L2-architecture.md)。本篇只聚焦 `crates/entrapment-cli` 一个 crate：它是 entrapment 分析的独立命令行入口，与 MCP Server 平行。包名 `protein-copilot-entrapment-cli`，二进制名 `entrapment`（Cargo.toml:9-11）。flag 名、类型、默认值均按 `src/main.rs` 核验，不臆造。

## 1. 用途 + 位置 + 依赖

`entrapment-cli` 不含任何分析算法，全部逻辑委托给 `entrapment-analysis`。它只做三件事：用 clap 解析参数、初始化 tracing、把子命令分派到 `run_analyze` / `run_report` / `run_inspect`（main.rs:121-143）。它与 `mcp-server` 互不依赖——同一套库既能被 LLM 经 MCP 调用，也能被人在 shell 里直接跑。`main` 先用 `tracing_subscriber::fmt()` 装日志（`RUST_LOG` 缺省取 `info`，main.rs:112-117）再 `Cli::parse()`；参数非法时由 clap 自动打印 usage 并以非零码退出，本 crate 不手写校验。

依赖（Cargo.toml:13-20）很薄：`protein-copilot-entrapment-analysis`（全部业务）、`clap`（derive 解析）、`serde_json`（config 快照序列化）、`tracing` + `tracing-subscriber`（结构化日志，`RUST_LOG` 可调）、`chrono`（运行时间戳）。

三条子命令职责：
- `analyze`：跑完整管线，读结果 + FASTA + config，写四份产物。
- `report`：从已分级 TSV 重生成 HTML——当前为占位，直接返回 "not yet implemented"（main.rs:302-307）。
- `inspect`：单肽探针——内部造一条 dummy `UnifiedPsm`、强制以 `Trap` 分组喂 `classify_single` 以触发相似度比对，打印该肽与目标库的最佳匹配肽/蛋白、mismatch 数与 delta-mass（main.rs:310-385）。

## 2. CLI 参数表

clap derive 默认把字段名转 kebab-case 长 flag，`#[arg(short, long)]` 取首字母为短 flag；每条子命令短 flag 各自独立。

`analyze`（main.rs:44-63）:

| flag | 类型 | 默认 | 含义 |
|---|---|---|---|
| `-r, --results` | String | 必填 | 搜索结果 .parquet/.tsv |
| `-c, --config` | String | 必填 | YAML 配置路径 |
| `-t, --target-fasta` | String | 必填 | 目标 FASTA 库 |
| `-f, --format` | Option<FormatArg> | None（按扩展名探测） | `diann-parquet` 或 `generic-tsv` |
| `-o, --out` | String | `output/entrapment` | 输出目录 |
| `--mzml-dir` | Option<PathBuf> | None | 给定则触发碎片溯源 + 多靶溯源 |

`report`（main.rs:65-72）: `-l, --classified <String>`（短 flag 显式设为 `l`，main.rs:67）、`-o, --out <Option<String>>`（默认 None；注释称将落在同目录 `entrapment_report.html`，但 run_report 未实现）。

`inspect`（main.rs:74-84）: `-p, --peptide <String>` 必填、`-t, --target-fasta <String>` 必填、`-c, --config <Option<String>>`（默认 None，回退 `SimilarityConfig::default()`）。

`FormatArg`（main.rs:88）是 `ValueEnum`，两变体 `DiannParquet`/`GenericTsv`，CLI 取值 `diann-parquet`/`generic-tsv`，`to_result_format` 映射到库的 `ResultFormat`（main.rs:98）。省略 `--format` 时由 `ResultFormat::from_path` 按扩展名定夺：`.parquet` 走 DiannParquet，`.tsv`/`.txt` 先做 pFind 探测、否则 GenericTsv，其它扩展名直接报 `LoaderError`（loader/mod.rs:37）。

## 3. 主流程 run_analyze

`run_analyze`（main.rs:156）是九步线性管线；错误统一 `Box<dyn Error>` 上抛，由 `main` 打印 `Error:` 并 `process::exit(1)`（main.rs:145-148）：

```
1 from_yaml      读 YAML -> EntrapmentConfig             config.rs:309
2 format         --format 或 from_path 按扩展名探测       loader/mod.rs:37
3 load_psms      结果文件 -> Vec<UnifiedPsm>             loader/mod.rs:69
4 Analyzer::new  config + FASTA -> 消化索引              lib.rs:53
5 classify_all   逐 PSM 打 L0-L4 -> Vec<ClassifiedPsm>   lib.rs:79
5b/5c 可选       --mzml-dir 触发单靶 + 多靶溯源(失败仅告警) lib.rs:192/454
6 create_dir_all 建输出目录
7 写产物         classified.tsv + razor_errors.tsv + run_metadata.json
8 render_report  自包含 HTML                             report.rs:130
9 stdout         打印 Summary(target/trap/L0-L4 计数)
```

第 7 步的 `run_metadata.json` 是可复现快照（`EntrapmentRunMetadata`，output.rs:25）：工具版本（`CARGO_PKG_VERSION`）、RFC3339 时间戳、输入与 FASTA 的路径 + sha256、config 的 serde_json 快照、PSM 与各级计数；其中两份 sha256 由 `file_sha256` 算出（output.rs:55），连同 config 快照构成可复现指纹，便于在 CI 或论文附录里比对同一次运行。`--mzml-dir` 缺省时 5b/5c 整段跳过，CLI 退化为"纯分级"；溯源即便报错也只 `eprintln` 告警、不中断管线（main.rs:200-203, 240）。5c 还会先用 `Tagger` 给每条 PSM 重打 group（target/trap/ambiguous）再喂多靶溯源（main.rs:212-216）。第 9 步的 Summary 字段取自 `analyser.summary` 返回的 `EntrapmentSummary`：total / target / trap / ambiguous 与 L0-L4 五级计数（main.rs:281-294）。

## 4. 简化源码片段

clap 定义（main.rs:31-94，删注释/属性）:

```rust
#[derive(Parser)]
#[command(name = "entrapment", about = "...classify trap PSMs by homology")]
struct Cli { #[command(subcommand)] command: Commands }

#[derive(Subcommand)]
enum Commands {
    Analyze {
        #[arg(short, long)] results: String,
        #[arg(short, long)] config: String,
        #[arg(short, long)] target_fasta: String,
        #[arg(short, long)] format: Option<FormatArg>,
        #[arg(short, long, default_value = "output/entrapment")] out: String,
        #[arg(long)] mzml_dir: Option<PathBuf>,
    },
    Report  { #[arg(short = 'l', long)] classified: String,
              #[arg(short, long)] out: Option<String> },
    Inspect { #[arg(short, long)] peptide: String,
              #[arg(short, long)] target_fasta: String,
              #[arg(short, long)] config: Option<String> },
}
```

FormatArg -> ResultFormat（main.rs:88-104）:

```rust
#[derive(Clone, ValueEnum)]
enum FormatArg { DiannParquet, GenericTsv }   // CLI: diann-parquet / generic-tsv

fn to_result_format(&self) -> ResultFormat {
    match self {
        Self::DiannParquet => ResultFormat::DiannParquet,
        Self::GenericTsv   => ResultFormat::GenericTsv,
    }
}
```

run_analyze 骨架（main.rs:156-296）:

```rust
fn run_analyze(...) -> Result<(), Box<dyn std::error::Error>> {
    let config = EntrapmentConfig::from_yaml(config_path)?;          // 1
    let format = match format_arg {                                 // 2
        Some(fa) => fa.to_result_format(),
        None     => ResultFormat::from_path(results_path)?,
    };
    let psms = loader::load_psms(results_path, &format, None)?;      // 3
    let analyser = EntrapmentAnalyzer::new(config.clone(), fasta_path)?; // 4
    let mut classified = analyser.classify_all(&psms)?;             // 5
    if let Some(dir) = mzml_dir {                                   // 5b/5c 可选, 非致命
        match trace_provenance_batch(&mut classified, dir, &config) {
            Ok(n)  => println!("traced {n}"),
            Err(e) => eprintln!("Warning: {e}"),
        }
        // trace_multi_target_provenance(&classified, &psms, &groups, dir, &config, &out_dir)
    }
    std::fs::create_dir_all(&out_dir)?;                            // 6
    write_classified_tsv(&classified, &out_dir.join("classified.tsv"))?;   // 7
    write_razor_errors_tsv(&classified, &out_dir.join("razor_errors.tsv"))?;
    let summary = analyser.summary(&classified);
    write_run_metadata(&metadata, &out_dir.join("run_metadata.json"))?;
    report::render_report(&summary, &classified,                   // 8
        &out_dir.join("entrapment_report.html"))?;
    // 9: println! 打印 Summary ...
    Ok(())
}
```

## 5. 用法示例

跑完整分析（按扩展名探测格式，用默认输出目录 `output/entrapment`）:

```
cargo run -p protein-copilot-entrapment-cli --offline -- \
  analyze -r results.parquet -c entrapment.yaml -t target.fasta
```

显式指定格式 + 输出目录 + mzML 溯源（长 flag 写法）:

```
cargo run -p protein-copilot-entrapment-cli --offline -- \
  analyze --results report.tsv --config cfg.yaml --target-fasta db.fasta \
          --format generic-tsv --out output/run1 --mzml-dir data/mzml
```

单肽探针:

```
cargo run -p protein-copilot-entrapment-cli --offline -- \
  inspect -p SAMPLEPEPTIDEK -t target.fasta
```

调日志级别用 `RUST_LOG=debug`（main.rs:112-117）。装好的二进制名是 `entrapment`，即 `entrapment analyze -r ...`。

## 6. 测试入口

```
cargo test -p protein-copilot-entrapment-cli --offline
```

该 crate 只有 `src/main.rs`、无 `#[cfg(test)]` 模块，故实测输出为 `running 0 tests ... test result: ok. 0 passed`（已核验）。CLI 的端到端覆盖落在 `integration-tests` 与 `entrapment-analysis` 的库测试，本 bin crate 只做参数解析与分派、不持有需单测的逻辑。

回到 [README](README.md)。
