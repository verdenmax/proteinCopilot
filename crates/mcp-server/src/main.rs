//! ProteinCopilot MCP Server — single binary entry point.
//!
//! Assembles all library crates and exposes them as MCP tools
//! over stdio transport for use with Copilot CLI / Claude Desktop.

use rmcp::transport::stdio;
use rmcp::ServiceExt;
use tracing_subscriber::EnvFilter;

mod tools;

use tools::ProteinCopilotServer;

#[tokio::main]
async fn main() {
    // Initialize tracing (respects RUST_LOG env var)
    tracing_subscriber::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .with_writer(std::io::stderr)
        .init();

    tracing::info!("Starting ProteinCopilot MCP Server");

    let server = ProteinCopilotServer::new();
    let service = server.serve(stdio()).await.unwrap();
    service.waiting().await.unwrap();
}
