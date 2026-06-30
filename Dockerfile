# Multi-stage build for the ProteinCopilot MCP server (single binary).
#
#   docker build -t protein-copilot-mcp .
#   docker run -i --rm protein-copilot-mcp           # speaks MCP over stdio
#
# Note: the MCP server communicates over stdio; MCP clients normally spawn the
# binary directly. Use `-i` so the container keeps stdin open.

FROM rust:1.85-bookworm AS builder
WORKDIR /src
COPY . .
RUN cargo build --release -p protein-copilot-mcp-server

FROM debian:bookworm-slim
RUN useradd --create-home --user-group app
COPY --from=builder /src/target/release/protein-copilot-mcp /usr/local/bin/protein-copilot-mcp
USER app
ENV RUST_LOG=info
ENTRYPOINT ["protein-copilot-mcp"]
