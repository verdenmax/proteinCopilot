#!/usr/bin/env node
// ProteinCopilot MCP Server — npm launcher shim
//
// Detects the current OS/arch/libc and spawns the matching platform-specific
// binary from an optionalDependency installed alongside this package.
//
// Reference: Biome's bin/biome shim pattern
//   https://github.com/biomejs/biome

'use strict';

const { platform, arch, env, version, release } = require('process');
const { execSync, spawnSync } = require('child_process');

// ── musl detection (Linux only) ────────────────────────────────────────────
function isMusl() {
  if (platform !== 'linux') return false;
  try {
    const stderr = execSync('ldd --version', { stdio: ['pipe', 'pipe', 'pipe'] });
    return String(stderr).includes('musl');
  } catch (e) {
    return String(e.stderr || '').includes('musl');
  }
}

// ── platform → package-relative binary path ────────────────────────────────
const PLATFORM_MAP = {
  'linux-x64': isMusl()
    ? 'protein-copilot-mcp-linux-x64-musl/protein-copilot-mcp'
    : 'protein-copilot-mcp-linux-x64-gnu/protein-copilot-mcp',
  'darwin-x64':   'protein-copilot-mcp-darwin-x64/protein-copilot-mcp',
  'darwin-arm64': 'protein-copilot-mcp-darwin-arm64/protein-copilot-mcp',
  'win32-x64':    'protein-copilot-mcp-win32-x64/protein-copilot-mcp.exe',
};

const platformKey = `${platform}-${arch}`;
const binRelPath = PLATFORM_MAP[platformKey];

if (!binRelPath) {
  console.error(
    `protein-copilot-mcp: unsupported platform "${platformKey}".\n` +
    'Pre-built binaries are available for linux-x64, darwin-x64, darwin-arm64, and win32-x64.\n' +
    'You can also install from source via Cargo:\n' +
    '  cargo install --git https://github.com/verdenmax/proteinCopilot.git -p protein-copilot-mcp-server\n' +
    'Or download a binary from:\n' +
    '  https://github.com/verdenmax/proteinCopilot/releases'
  );
  process.exit(1);
}

let binPath;
try {
  binPath = require.resolve(binRelPath);
} catch (e) {
  console.error(
    `protein-copilot-mcp: failed to locate platform binary "${binRelPath}".\n` +
    'This usually means the platform-specific optional dependency was not installed.\n' +
    'Try re-running: npm install protein-copilot-mcp\n' +
    'If the problem persists, please report it at:\n' +
    '  https://github.com/verdenmax/proteinCopilot/issues'
  );
  process.exit(1);
}

// Ensure executable on Unix (npm may strip the execute bit)
if (platform !== 'win32') {
  try {
    const fs = require('fs');
    fs.chmodSync(binPath, 0o755);
  } catch {
    // best-effort; the file may already be executable
  }
}

// Forward all arguments + stdio to the Rust binary
const result = spawnSync(binPath, process.argv.slice(2), {
  stdio: 'inherit',
  env: {
    ...env,
    PROTEIN_DISTRIBUTION: 'npm',
    JS_RUNTIME_VERSION: version,
    JS_RUNTIME_NAME: release.name,
  },
});

if (result.error) {
  throw result.error;
}

process.exit(result.status ?? 1);
