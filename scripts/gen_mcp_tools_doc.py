#!/usr/bin/env python3
"""Generate docs/mcp-tools.md from the live MCP binary's tools/list (JSON Schema).

This keeps the public tool contract authoritative: it is produced from exactly the
inputSchema/outputSchema the binary sends to MCP clients, so it never drifts from code.

Usage:
    cargo build --release -p protein-copilot-mcp-server --offline
    python3 scripts/gen_mcp_tools_doc.py
"""
import json
import subprocess
import sys
import threading
import time
import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "target" / "release" / "protein-copilot-mcp"
OUT = ROOT / "docs" / "mcp-tools.md"

CATEGORIES = [
    ("读取谱图", ["read_spectra", "get_spectrum"]),
    ("参数推荐", ["recommend_params", "list_presets", "prepare_search"]),
    ("搜索生命周期", ["run_search", "get_search_status", "cancel_search", "check_engine", "diagnose_search"]),
    ("结果摘要与导出", ["generate_summary", "export_results", "list_searches"]),
    ("蛋白推断", ["infer_proteins"]),
    ("谱图注释与可视化", ["annotate_spectrum", "extract_xic"]),
    ("DIA 数据提取", ["extract_dia_precursors", "extract_spectrum_precursors", "get_dia_cache_status"]),
    ("外部结果导入", ["import_search_results"]),
    ("FASTA 数据库", ["list_databases", "download_database", "get_database_info"]),
    ("entrapment 分析", ["classify_entrapment_hits", "analyze_entrapment_stats", "find_similar_targets", "annotate_provenance"]),
]


def fetch_tools():
    if not BIN.exists():
        sys.exit(f"binary not found: {BIN}\nrun: cargo build --release -p protein-copilot-mcp-server --offline")
    p = subprocess.Popen([str(BIN)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True, bufsize=1)

    def send(o):
        p.stdin.write(json.dumps(o) + "\n")
        p.stdin.flush()

    def readline(timeout=20):
        box = {}
        t = threading.Thread(target=lambda: box.setdefault("l", p.stdout.readline()))
        t.daemon = True
        t.start()
        t.join(timeout)
        return box.get("l")

    send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
          "params": {"protocolVersion": "2024-11-05", "capabilities": {},
                     "clientInfo": {"name": "doc-gen", "version": "0"}}})
    readline()
    send({"jsonrpc": "2.0", "method": "notifications/initialized"})
    time.sleep(0.3)
    send({"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}})
    resp = readline()
    p.terminate()
    if not resp or '"tools"' not in resp:
        sys.exit(f"tools/list failed: {resp!r}")
    return json.loads(resp)["result"]["tools"]


def ref_name(ref):
    return ref.split("/")[-1]


def type_str(v):
    if not isinstance(v, dict):
        return "any"
    if "$ref" in v:
        n = ref_name(v["$ref"])
        return f"[`{n}`](#类型-{n.lower()})"
    for key in ("anyOf", "allOf", "oneOf"):
        if key in v:
            parts = [s for s in v[key] if not (isinstance(s, dict) and s.get("type") == "null")]
            nullable = len(parts) != len(v[key])
            rendered = " \\| ".join(type_str(s) for s in parts) if parts else "null"
            return rendered + (" *(可空)*" if nullable else "")
    t = v.get("type")
    if isinstance(t, list):
        base = [x for x in t if x != "null"]
        s = " \\| ".join(base) if base else "null"
        return s + (" *(可空)*" if "null" in t else "")
    if t == "array":
        return f"array&lt;{type_str(v.get('items', {}))}&gt;"
    if "enum" in v:
        return "enum(" + ", ".join(f"`{e}`" for e in v["enum"]) + ")"
    return t or "object"


def fmt_default(v):
    if isinstance(v, dict) and v.get("default") is not None:
        return f"`{json.dumps(v['default'], ensure_ascii=False)}`"
    return "—"


def render(tools):
    byname = {t["name"]: t for t in tools}
    out = []
    out.append("# MCP 工具参考（ProteinCopilot）\n")
    out.append("> 本文件由 `scripts/gen_mcp_tools_doc.py` 从 release 二进制的 `tools/list`"
               "（JSON Schema）**自动生成**，是面向调用方的权威接口契约——无需阅读源码。"
               "工具签名变更后请重跑该脚本。\n")
    out.append(f"> 工具总数：**{len(tools)}**。传输：JSON-RPC 2.0 over stdio。"
               f"生成时间：{datetime.date.today()}。\n")

    out.append("\n## 启动与接入\n")
    out.append("```bash\n# 直接运行已编译的二进制（推荐发布形态）\n./protein-copilot-mcp\n\n"
               "# 或从源码运行\ncargo run --release -p protein-copilot-mcp-server\n```\n")
    out.append("命令行自检（无需客户端，直接在终端查看工具契约）：\n")
    out.append("```bash\n./protein-copilot-mcp --list-tools          # 文本目录：参数/类型/范围/默认/输出\n"
               "./protein-copilot-mcp --list-tools --json   # 完整 JSON Schema（机器可读）\n"
               "./protein-copilot-mcp --help                # 用法\n```\n")
    out.append("在 MCP 客户端（Copilot CLI / Claude Desktop 等）中登记：\n")
    out.append("```json\n{\n  \"mcpServers\": {\n    \"protein-copilot\": {\n"
               "      \"command\": \"/path/to/protein-copilot-mcp\",\n"
               "      \"env\": { \"RUST_LOG\": \"info\" }\n    }\n  }\n}\n```\n")
    out.append("- 所有工具的输入/输出均为结构化 JSON，类型见下方「参数」表与「共享数据类型」。\n"
               "- 描述文本为二进制 `#[schemars]` 原文（即客户端实际收到的内容），故为英文。\n"
               "- 搜索为异步：`run_search` 立即返回 `run_id`，用 `get_search_status` 轮询，"
               "完成后 `generate_summary` / `export_results` / `infer_proteins`。\n")

    out.append("\n## 工具索引\n")
    for cat, names in CATEGORIES:
        line = ", ".join(f"[`{n}`](#{n})" for n in names if n in byname)
        out.append(f"- **{cat}**：{line}")
    listed = {n for _, ns in CATEGORIES for n in ns}
    extra = [t["name"] for t in tools if t["name"] not in listed]
    if extra:
        out.append("- **其它**：" + ", ".join(f"[`{n}`](#{n})" for n in extra))

    defs = {}
    for t in tools:
        for dn, dv in (t["inputSchema"].get("$defs") or {}).items():
            defs.setdefault(dn, dv)

    out.append("\n---\n\n## 工具详情\n")
    rendered_tools = set()

    def emit_tool(n):
        t = byname[n]
        s = t["inputSchema"]
        out.append(f"\n#### `{n}`\n")
        out.append((t["description"] or "").strip() + "\n")
        props = s.get("properties", {})
        required = s.get("required", [])
        if props:
            out.append("\n| 参数 | 必填 | 类型 | 默认 | 说明 |")
            out.append("|------|------|------|------|------|")
            for k, v in sorted(props.items(), key=lambda kv: (kv[0] not in required, kv[0])):
                req = "是" if k in required else "否"
                desc = (v.get("description") or "").replace("\n", " ").replace("|", "\\|").strip()
                out.append(f"| `{k}` | {req} | {type_str(v)} | {fmt_default(v)} | {desc} |")
        else:
            out.append("\n*无参数。*")
        osch = t.get("outputSchema") or {}
        otitle = osch.get("title") or osch.get("type") or "object"
        tail = f" — {osch.get('description', '').strip()}" if osch.get("description") else ""
        out.append(f"\n**输出**：`{otitle}`{tail}\n")
        rendered_tools.add(n)

    for cat, names in CATEGORIES:
        out.append(f"\n### {cat}\n")
        for n in names:
            if n in byname:
                emit_tool(n)
    leftover = [t["name"] for t in tools if t["name"] not in rendered_tools]
    if leftover:
        out.append("\n### 其它\n")
        for n in leftover:
            emit_tool(n)

    out.append("\n---\n\n## 共享数据类型\n")
    out.append("\n参数与输出中引用的复合类型定义如下（枚举列出全部取值，结构体列出字段）。\n")
    for dn in sorted(defs):
        dv = defs[dn]
        out.append(f"\n### 类型 {dn}\n")
        if dv.get("description"):
            out.append(dv["description"].strip() + "\n")
        if "enum" in dv:
            out.append("\n枚举取值：" + ", ".join(f"`{e}`" for e in dv["enum"]) + "\n")
        elif "oneOf" in dv or "anyOf" in dv:
            out.append("\n变体：\n")
            for var in (dv.get("oneOf") or dv.get("anyOf")):
                vdesc = (" — " + var["description"].strip()) if var.get("description") else ""
                if "const" in var:
                    out.append(f"- `{var['const']}`{vdesc}")
                elif var.get("enum"):
                    for e in var["enum"]:
                        out.append(f"- `{e}`{vdesc}")
                elif var.get("properties"):
                    for tag, tv in var["properties"].items():
                        inner = tv.get("properties") if isinstance(tv, dict) else None
                        if inner:
                            flds = "; ".join(
                                f"`{fk}`: {type_str(fv)}" + ("" if fk in tv.get("required", []) else " *(可选)*")
                                for fk, fv in inner.items())
                            out.append(f"- `{tag}`{vdesc}：{{ {flds} }}")
                        else:
                            out.append(f"- `{tag}`{vdesc}")
                elif var.get("type"):
                    out.append(f"- {type_str(var)}{vdesc}")
        elif dv.get("type") == "object" and dv.get("properties"):
            out.append("\n| 字段 | 必填 | 类型 | 默认 | 说明 |")
            out.append("|------|------|------|------|------|")
            req = dv.get("required", [])
            for k, v in dv["properties"].items():
                r = "是" if k in req else "否"
                desc = (v.get("description") or "").replace("\n", " ").replace("|", "\\|").strip()
                out.append(f"| `{k}` | {r} | {type_str(v)} | {fmt_default(v)} | {desc} |")
        else:
            out.append(f"\n类型：{type_str(dv)}\n")

    return "\n".join(out) + "\n"


def main():
    tools = fetch_tools()
    OUT.write_text(render(tools), encoding="utf-8")
    print(f"wrote {OUT} ({len(tools)} tools)")


if __name__ == "__main__":
    main()
