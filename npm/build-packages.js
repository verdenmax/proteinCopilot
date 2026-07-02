#!/usr/bin/env node
// build-packages.js — assemble npm packages from CI-built platform binaries.
//
// Usage:
//   node npm/build-packages.js --version 0.1.0 --binaries-dir ./artifacts --output-dir npm/dist
//
// Expects artifacts directory laid out by actions/download-artifact@v4:
//   artifacts/npm-binary-x86_64-unknown-linux-gnu/protein-copilot-mcp
//   artifacts/npm-binary-x86_64-unknown-linux-musl/protein-copilot-mcp
//   artifacts/npm-binary-x86_64-apple-darwin/protein-copilot-mcp
//   artifacts/npm-binary-aarch64-apple-darwin/protein-copilot-mcp
//   artifacts/npm-binary-x86_64-pc-windows-msvc/protein-copilot-mcp.exe
//
// Produces:
//   <output-dir>/
//     protein-copilot-mcp/             # main package (shim + optionalDeps)
//     protein-copilot-mcp-linux-x64-gnu/
//     protein-copilot-mcp-linux-x64-musl/
//     protein-copilot-mcp-darwin-x64/
//     protein-copilot-mcp-darwin-arm64/
//     protein-copilot-mcp-win32-x64/

'use strict';

const fs = require('fs');
const path = require('path');
const { spawnSync } = require('child_process');

// ── CLI arg parsing ─────────────────────────────────────────────────────────
function parseArgs() {
  const args = process.argv.slice(2);
  const opts = { version: null, binariesDir: null, outputDir: null };
  for (let i = 0; i < args.length; i++) {
    if (args[i] === '--version' && i + 1 < args.length) {
      opts.version = args[++i];
    } else if (args[i] === '--binaries-dir' && i + 1 < args.length) {
      opts.binariesDir = args[++i];
    } else if (args[i] === '--output-dir' && i + 1 < args.length) {
      opts.outputDir = args[++i];
    }
  }

  if (!opts.version || !opts.binariesDir || !opts.outputDir) {
    console.error('Usage: node build-packages.js --version <semver> --binaries-dir <dir> --output-dir <dir>');
    process.exit(1);
  }

  if (!fs.existsSync(opts.binariesDir)) {
    console.error(`Binaries directory not found: ${opts.binariesDir}`);
    process.exit(1);
  }

  return opts;
}

// ── Rust target → npm platform package metadata ─────────────────────────────
// Each entry: { pkgName, os, cpu, libc?, binaryName }
const TARGET_MAP = {
  'x86_64-unknown-linux-gnu':  { pkgName: 'protein-copilot-mcp-linux-x64-gnu',  os: 'linux',  cpu: 'x64',   libc: 'glibc', binaryName: 'protein-copilot-mcp' },
  'x86_64-unknown-linux-musl': { pkgName: 'protein-copilot-mcp-linux-x64-musl', os: 'linux',  cpu: 'x64',   libc: 'musl',  binaryName: 'protein-copilot-mcp' },
  'x86_64-apple-darwin':      { pkgName: 'protein-copilot-mcp-darwin-x64',     os: 'darwin', cpu: 'x64',                 binaryName: 'protein-copilot-mcp' },
  'aarch64-apple-darwin':     { pkgName: 'protein-copilot-mcp-darwin-arm64',   os: 'darwin', cpu: 'arm64',               binaryName: 'protein-copilot-mcp' },
  'x86_64-pc-windows-msvc':   { pkgName: 'protein-copilot-mcp-win32-x64',      os: 'win32',  cpu: 'x64',                 binaryName: 'protein-copilot-mcp.exe' },
};

// ── Helpers ──────────────────────────────────────────────────────────────────
function copyFile(src, dest) {
  fs.copyFileSync(src, dest);
  console.log(`  cp ${src} -> ${dest}`);
}

function writeJson(filePath, obj) {
  fs.writeFileSync(filePath, JSON.stringify(obj, null, 2) + '\n');
  console.log(`  write ${filePath}`);
}

// ── Main ─────────────────────────────────────────────────────────────────────
function main() {
  const { version, binariesDir, outputDir } = parseArgs();
  const repoRoot = path.resolve(__dirname, '..');

  console.log(`Building npm packages v${version}`);
  console.log(`  binaries: ${binariesDir}`);
  console.log(`  output:   ${outputDir}`);

  // Clean output
  if (fs.existsSync(outputDir)) {
    fs.rmSync(outputDir, { recursive: true });
  }
  fs.mkdirSync(outputDir, { recursive: true });

  // Read platform template
  const templatePath = path.join(__dirname, 'platform-template', 'package.json');
  const template = JSON.parse(fs.readFileSync(templatePath, 'utf-8'));

  const mainPkgName = 'protein-copilot-mcp';
  const mainPkgDir = path.join(outputDir, mainPkgName);
  const optionalDeps = {};

  // Build each platform package
  for (const [target, meta] of Object.entries(TARGET_MAP)) {
    const artifactName = `npm-binary-${target}`;
    const binaryInDir = path.join(binariesDir, artifactName);
    const binarySrc = path.join(binaryInDir, meta.binaryName);

    if (!fs.existsSync(binarySrc)) {
      console.error(`  WARNING: binary not found at ${binarySrc} — skipping ${meta.pkgName}`);
      continue;
    }

    const pkgDir = path.join(outputDir, meta.pkgName);
    fs.mkdirSync(pkgDir, { recursive: true });

    // Copy binary
    const binaryDest = path.join(pkgDir, meta.binaryName);
    copyFile(binarySrc, binaryDest);

    // Ensure executable on Unix
    if (meta.os !== 'win32') {
      try { fs.chmodSync(binaryDest, 0o755); } catch { /* best-effort */ }
    }

    // Generate platform package.json
    const pkgJson = { ...template };
    pkgJson.name = meta.pkgName;
    pkgJson.version = version;
    pkgJson.os = [meta.os];
    pkgJson.cpu = [meta.cpu];
    if (meta.libc) {
      pkgJson.libc = [meta.libc];
    }

    writeJson(path.join(pkgDir, 'package.json'), pkgJson);
    optionalDeps[meta.pkgName] = version;

    console.log(`  built ${meta.pkgName}`);
  }

  // Build main package
  console.log(`\nBuilding main package: ${mainPkgName}`);
  fs.mkdirSync(path.join(mainPkgDir, 'bin'), { recursive: true });

  // Copy shim
  const shimSrc = path.join(__dirname, mainPkgName, 'bin', 'protein-copilot-mcp.js');
  copyFile(shimSrc, path.join(mainPkgDir, 'bin', 'protein-copilot-mcp.js'));

  // Read main template and fill in version + optionalDeps
  const mainTemplatePath = path.join(__dirname, mainPkgName, 'package.json');
  const mainPkgJson = JSON.parse(fs.readFileSync(mainTemplatePath, 'utf-8'));
  mainPkgJson.version = version;
  mainPkgJson.optionalDependencies = optionalDeps;

  writeJson(path.join(mainPkgDir, 'package.json'), mainPkgJson);

  // Copy README and LICENSE
  copyFile(path.join(repoRoot, 'README.md'), path.join(mainPkgDir, 'README.md'));
  copyFile(path.join(repoRoot, 'LICENSE'), path.join(mainPkgDir, 'LICENSE'));

  console.log(`\nDone! ${Object.keys(optionalDeps).length + 1} packages written to ${outputDir}`);
  console.log('Platform packages:');
  for (const pkg of Object.keys(optionalDeps)) {
    console.log(`  ${pkg}@${version}`);
  }
  console.log(`  ${mainPkgName}@${version}`);
}

main();
