//! ProteinCopilot MCP Server — single binary entry point.
//!
//! Assembles all library crates and exposes them as MCP tools
//! over stdio transport for use with Copilot CLI / Claude Desktop.

use rmcp::transport::stdio;
use rmcp::ServiceExt;
use tracing_subscriber::{fmt, layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

mod history;
mod tools;

use tools::ProteinCopilotServer;

#[tokio::main]
async fn main() {
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
