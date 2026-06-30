//! ProteinCopilot MCP Server — single binary entry point.
//!
//! Assembles all library crates and exposes them as MCP tools
//! over stdio transport for use with Copilot CLI / Claude Desktop.

use rmcp::transport::stdio;
use rmcp::ServiceExt;
use tracing_subscriber::{fmt, layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

mod catalog;
mod history;
mod tools;

use tools::ProteinCopilotServer;

#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage();
        return;
    }
    if args.iter().any(|a| a == "-V" || a == "--version") {
        println!("protein-copilot-mcp {}", env!("CARGO_PKG_VERSION"));
        return;
    }
    if args.iter().any(|a| a == "--list-tools") {
        let mut tools = ProteinCopilotServer::new().list_tools();
        tools.sort_by(|a, b| a.name.cmp(&b.name));
        if args.iter().any(|a| a == "--json") {
            match serde_json::to_string_pretty(&tools) {
                Ok(s) => println!("{s}"),
                Err(e) => {
                    eprintln!("failed to serialize tools: {e}");
                    std::process::exit(1);
                }
            }
        } else {
            print!("{}", catalog::format_catalog(&tools));
        }
        return;
    }

    // Initialize tracing (respects RUST_LOG env var)
    // PROTEIN_LOG_JSON=1 switches to structured JSON output
    let env_filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));

    let use_json = std::env::var("PROTEIN_LOG_JSON").is_ok_and(|v| v == "1");

    if use_json {
        tracing_subscriber::registry()
            .with(env_filter)
            .with(
                fmt::layer()
                    .json()
                    .with_writer(std::io::stderr)
                    .with_span_events(fmt::format::FmtSpan::CLOSE)
                    .with_target(true)
                    .with_timer(fmt::time::uptime()),
            )
            .init();
    } else {
        tracing_subscriber::registry()
            .with(env_filter)
            .with(
                fmt::layer()
                    .with_writer(std::io::stderr)
                    .with_span_events(fmt::format::FmtSpan::CLOSE)
                    .with_target(true)
                    .with_timer(fmt::time::uptime()),
            )
            .init();
    }

    tracing::info!("Starting ProteinCopilot MCP Server");

    let server = ProteinCopilotServer::new();
    let service = match server.serve(stdio()).await {
        Ok(s) => s,
        Err(e) => {
            tracing::error!("Failed to start MCP server: {e}");
            std::process::exit(1);
        }
    };
    if let Err(e) = service.waiting().await {
        tracing::error!("MCP server exited with error: {e}");
        std::process::exit(1);
    }
}

/// Prints CLI usage. The server normally runs as an MCP stdio service with no
/// arguments; these flags are convenience inspectors for the published binary.
fn print_usage() {
    println!(
        "protein-copilot-mcp {ver} — ProteinCopilot MCP Server

USAGE:
    protein-copilot-mcp [FLAGS]

With no flags, runs as an MCP server over stdio (for Copilot CLI / Claude Desktop).

FLAGS:
    --list-tools         Print the tool catalog (name, params, types, ranges, output) and exit
    --list-tools --json  Print the full tool JSON Schema (machine-readable) and exit
    -h, --help           Print this help and exit
    -V, --version        Print version and exit

ENV:
    RUST_LOG             Log level (default: info), e.g. RUST_LOG=debug
    PROTEIN_LOG_JSON=1   Emit logs as JSON",
        ver = env!("CARGO_PKG_VERSION")
    );
}
