# Architecture HTML Showcase Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a presentation-friendly `docs/architecture.html` that explains ProteinCopilot's architecture with a scroll narrative, lightweight hotspot interactions, and a concise MCP explainer.

**Architecture:** The implementation stays in a single self-contained HTML file with inline CSS and vanilla JavaScript so it can be opened directly in a browser or served by any static file server. The page is organized as a presentation flow — hero, MCP explainer, layered architecture, workflow pipeline, and extensibility — with click/hover hotspots that update a shared detail panel instead of hiding core information behind interaction.

**Tech Stack:** HTML5, inline CSS, vanilla JavaScript, browser-native SVG/semantic layout, Markdown docs

**Spec:** `docs/superpowers/specs/2026-04-08-architecture-html-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `docs/architecture.html` | Self-contained architecture presentation page with all markup, styling, copy, and interaction logic |
| Modify | `README.md:96-97` | Add a discoverable link to the new architecture presentation page next to the detailed Markdown architecture doc |
| Modify | `docs/architecture.md:1-6` | Add a short note pointing readers to the presentation HTML while keeping the Markdown file as the detailed reference |

---

### Task 1: Scaffold the single-file architecture page

**Files:**
- Create: `docs/architecture.html`
- Test: manual smoke check of `docs/architecture.html`

- [ ] **Step 1: Run the missing-file smoke check**

Run:

```bash
test ! -f docs/architecture.html && echo "docs/architecture.html missing"
```

Expected: prints `docs/architecture.html missing`.

- [ ] **Step 2: Create `docs/architecture.html` with the shell, section anchors, and base theme**

Write this file:

```html
<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ProteinCopilot Architecture</title>
  <style>
    :root {
      --bg: #08111f;
      --bg-soft: #0f1b33;
      --surface: #ffffff;
      --surface-soft: #f8fafc;
      --ink: #0f172a;
      --ink-soft: #475569;
      --border: #dbe3ef;
      --llm: #7c3aed;
      --mcp: #0891b2;
      --rust: #2563eb;
      --ext: #ea580c;
      --shadow: 0 18px 40px rgba(15, 23, 42, 0.12);
      --radius: 24px;
      --max-width: 1180px;
    }

    * { box-sizing: border-box; }

    html { scroll-behavior: smooth; }

    body {
      margin: 0;
      font-family: Inter, "Noto Sans SC", "PingFang SC", sans-serif;
      color: var(--ink);
      background: linear-gradient(180deg, #020617 0%, #08111f 26%, #eef4ff 26%, #eef4ff 100%);
    }

    a { color: inherit; }

    .skip-link {
      position: absolute;
      left: 16px;
      top: -48px;
      padding: 10px 14px;
      border-radius: 999px;
      background: #fff;
    }

    .skip-link:focus { top: 16px; }

    .hero,
    .page-section,
    .page-footer {
      width: min(var(--max-width), calc(100% - 32px));
      margin: 0 auto;
    }

    .hero {
      padding: 56px 0 48px;
      color: #fff;
    }

    .hero-card {
      padding: 40px;
      border-radius: 32px;
      background:
        radial-gradient(circle at top right, rgba(124, 58, 237, 0.24), transparent 30%),
        radial-gradient(circle at left center, rgba(8, 145, 178, 0.28), transparent 24%),
        linear-gradient(135deg, #0f172a, #111c34 62%, #16233f);
      box-shadow: 0 24px 60px rgba(2, 6, 23, 0.45);
    }

    .eyebrow {
      margin: 0 0 12px;
      font-size: 13px;
      font-weight: 700;
      letter-spacing: 0.12em;
      text-transform: uppercase;
      opacity: 0.78;
    }

    .hero h1,
    .page-section h2 {
      margin: 0;
    }

    .hero p {
      max-width: 720px;
      line-height: 1.7;
    }

    .hero-nav {
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      margin-top: 24px;
    }

    .hero-nav a,
    .chip {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      padding: 10px 14px;
      border: 1px solid rgba(255, 255, 255, 0.16);
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.08);
      text-decoration: none;
    }

    main {
      padding: 8px 0 64px;
    }

    .page-section {
      margin-top: 24px;
      padding: 32px;
      border: 1px solid var(--border);
      border-radius: var(--radius);
      background: rgba(255, 255, 255, 0.9);
      box-shadow: var(--shadow);
    }

    .page-section p {
      color: var(--ink-soft);
      line-height: 1.7;
    }

    .section-grid {
      display: grid;
      gap: 24px;
    }

    .detail-panel {
      margin-top: 24px;
      padding: 20px;
      border-radius: 20px;
      background: var(--surface-soft);
      border: 1px solid var(--border);
    }

    .reveal {
      opacity: 1;
      transform: none;
    }

    @media (max-width: 800px) {
      .hero-card,
      .page-section { padding: 24px; }
    }
  </style>
</head>
<body>
  <a class="skip-link" href="#main">Skip to content</a>

  <header id="hero" class="hero">
    <div class="hero-card">
      <p class="eyebrow">ProteinCopilot</p>
      <h1>架构总览</h1>
      <p>Rust 负责确定性科学计算，LLM 负责编排、理解用户意图与解释结果。MCP Server 把两者连接成一条可扩展的蛋白质组学工作流。</p>
      <nav class="hero-nav" aria-label="Architecture sections">
        <a href="#mcp">MCP 是什么</a>
        <a href="#layers">分层架构</a>
        <a href="#pipeline">工作流</a>
        <a href="#extensibility">扩展路径</a>
      </nav>
    </div>
  </header>

  <main id="main">
    <section id="mcp" class="page-section reveal">
      <p class="eyebrow">Section 01</p>
      <h2>MCP 是什么？</h2>
      <p>MCP 是 LLM 客户端与 Rust 工具之间的结构化协议层：LLM 不直接做数值计算，而是通过标准化 tool 调用把任务交给确定性代码执行。</p>
    </section>

    <section id="layers" class="page-section reveal">
      <p class="eyebrow">Section 02</p>
      <h2>分层架构</h2>
      <p>这一层会展示用户 / MCP Client / LLM、ProteinCopilot MCP Server、Rust crates 与搜索引擎之间的职责边界。</p>
    </section>

    <section id="pipeline" class="page-section reveal">
      <p class="eyebrow">Section 03</p>
      <h2>端到端工作流</h2>
      <p>这一层会把 read spectra、recommend params、run search、generate/export report 串成一条从输入到结果的故事线。</p>
    </section>

    <section id="extensibility" class="page-section reveal">
      <p class="eyebrow">Section 04</p>
      <h2>可扩展性</h2>
      <p>这一层会说明如何新增搜索引擎 adapter、增加 MCP 模块，以及在未来把领域拆成更多 MCP Server。</p>
    </section>

    <aside class="page-section detail-panel" id="detail-panel" aria-live="polite">
      <p class="eyebrow">Detail</p>
      <h2 id="detail-title">讲解提示</h2>
      <p id="detail-body">点击架构图中的热点后，这里会显示对应模块的职责说明。即使不点击，页面本身也必须完整可读。</p>
    </aside>
  </main>

  <script>
    const detailData = {};
  </script>
</body>
</html>
```

- [ ] **Step 3: Run the section-anchor check**

Run:

```bash
rg -n '<header id="hero"|<section id="mcp"|<section id="layers"|<section id="pipeline"|<section id="extensibility"|id="detail-panel"' docs/architecture.html
```

Expected: one match for the hero plus five section/detail IDs.

- [ ] **Step 4: Commit**

```bash
git add docs/architecture.html
git commit -m "docs: scaffold architecture presentation page

Create the single-file architecture HTML shell with the approved section
order and base visual theme.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 2: Add the hero, MCP explainer, and layered architecture content

**Files:**
- Modify: `docs/architecture.html`
- Test: manual browser smoke check of `docs/architecture.html`

- [ ] **Step 1: Run the content-presence check before editing**

Run:

```bash
rg -n 'MCP 不是 LLM 本身|data-panel="llm-layer"|data-panel="mcp-server"|data-panel="core-crate"' docs/architecture.html
```

Expected: no matches.

- [ ] **Step 2: Replace the MCP and layered-architecture placeholders with real content and hotspots**

In `docs/architecture.html`, replace the current `#mcp` and `#layers` sections with:

```html
<section id="mcp" class="page-section reveal">
  <p class="eyebrow">Section 01</p>
  <h2>MCP 是什么？</h2>
  <p>MCP 不是 LLM 本身，而是 LLM Client 与 ProteinCopilot Rust 能力之间的结构化桥梁。它把“理解意图”和“执行确定性工具”分开，让调用过程可解释、可组合、可扩展。</p>
  <div class="section-grid mcp-grid">
    <article class="info-card">
      <h3>Client</h3>
      <p>Copilot CLI / Claude Desktop 读取 Agent 与 Prompt，并代表用户发起 tool 调用。</p>
    </article>
    <article class="info-card">
      <h3>Protocol</h3>
      <p>JSON-RPC 风格的结构化消息，把输入、输出和错误都稳定地交给 Rust 处理。</p>
    </article>
    <article class="info-card">
      <h3>Tool Bridge</h3>
      <p>ProteinCopilot MCP Server 把 read_spectra、run_search、generate_summary 等能力暴露给 LLM。</p>
    </article>
  </div>
</section>

<section id="layers" class="page-section reveal">
  <p class="eyebrow">Section 02</p>
  <h2>分层架构：从 Agent 到 crate 的责任边界</h2>
  <p>页面要清楚表达两条主线：LLM 负责编排和解释，Rust 负责所有确定性计算；同时 MCP Server 作为统一入口连接上层 Agent 和下层 crate。</p>
  <div class="section-grid architecture-layout">
    <div class="stack">
      <button class="hotspot layer llm-layer" type="button" data-panel="llm-layer" aria-controls="detail-panel" aria-pressed="false">
        <strong>用户 / MCP Client / Agent & Prompt</strong>
        <span>理解意图、规划步骤、解释结果</span>
      </button>
      <button class="hotspot layer mcp-layer" type="button" data-panel="mcp-server" aria-controls="detail-panel" aria-pressed="false">
        <strong>ProteinCopilot MCP Server</strong>
        <span>统一 tool 入口，负责协议层与能力注册</span>
      </button>
      <div class="crate-grid">
        <button class="hotspot crate-card" type="button" data-panel="core-crate" aria-controls="detail-panel" aria-pressed="false">core</button>
        <button class="hotspot crate-card" type="button" data-panel="spectrum-io-crate" aria-controls="detail-panel" aria-pressed="false">spectrum-io</button>
        <button class="hotspot crate-card" type="button" data-panel="param-recommend-crate" aria-controls="detail-panel" aria-pressed="false">param-recommend</button>
        <button class="hotspot crate-card" type="button" data-panel="search-engine-crate" aria-controls="detail-panel" aria-pressed="false">search-engine</button>
        <button class="hotspot crate-card" type="button" data-panel="report-crate" aria-controls="detail-panel" aria-pressed="false">report</button>
        <button class="hotspot crate-card" type="button" data-panel="dia-extraction-crate" aria-controls="detail-panel" aria-pressed="false">dia-extraction</button>
      </div>
      <button class="hotspot layer ext-layer" type="button" data-panel="search-adapter" aria-controls="detail-panel" aria-pressed="false">
        <strong>搜索引擎 / Adapter</strong>
        <span>SimpleSearch 当前可用，pFind adapter 预留扩展位</span>
      </button>
    </div>
    <div class="takeaway-card">
      <h3>核心结论</h3>
      <ul>
        <li>Rust 负责解析、校验、搜索、聚合、导出</li>
        <li>LLM 负责意图理解、工作流编排、自然语言解释</li>
        <li>MCP Server 是统一能力入口，不让 LLM 直接触碰数值计算</li>
      </ul>
    </div>
  </div>
</section>
```

- [ ] **Step 3: Run the content-presence check again**

Run:

```bash
rg -n 'MCP 不是 LLM 本身|data-panel="llm-layer"|data-panel="mcp-server"|data-panel="core-crate"|ProteinCopilot MCP Server' docs/architecture.html
```

Expected: matches for the MCP explainer copy and the hotspot attributes.

- [ ] **Step 4: Start a local static server and verify the first half of the page visually**

Run:

```bash
python -m http.server 8000 >/tmp/protein-copilot-architecture-html.log 2>&1 & echo $!
```

Expected: prints a numeric PID. Then open `http://localhost:8000/docs/architecture.html` and verify:

- the hero reads cleanly on a projector-sized viewport
- the MCP section is understandable in under one screen
- the layered architecture section clearly separates the LLM layer, the MCP layer, and the Rust crates

After verification, stop the server with:

```bash
kill <PID>
```

- [ ] **Step 5: Commit**

```bash
git add docs/architecture.html
git commit -m "docs: add architecture hero and layered stack

Fill in the hero, MCP quick explainer, and layered architecture content
with hotspot targets for the shared detail panel.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 3: Add the workflow and extensibility sections

**Files:**
- Modify: `docs/architecture.html`
- Test: manual browser smoke check of `docs/architecture.html`

- [ ] **Step 1: Run the workflow/extensibility presence check before editing**

Run:

```bash
rg -n 'read_spectra|recommend_params|run_search|generate_summary|data-panel="new-adapter"|data-panel="split-server"' docs/architecture.html
```

Expected: no matches.

- [ ] **Step 2: Replace the workflow and extensibility placeholders with the approved narrative**

In `docs/architecture.html`, replace the current `#pipeline` and `#extensibility` sections with:

```html
<section id="pipeline" class="page-section reveal">
  <p class="eyebrow">Section 03</p>
  <h2>端到端工作流：从谱图到可解释结果</h2>
  <p>这一段不要堆所有 API 细节，只强调平台如何把研究员的问题转成一条稳定的处理链路。</p>
  <div class="pipeline-steps">
    <button class="hotspot pipeline-step" type="button" data-panel="read-spectra" aria-controls="detail-panel" aria-pressed="false">
      <strong>01 · read_spectra</strong>
      <span>读取 mzML / mgf，先拿到数据摘要</span>
    </button>
    <button class="hotspot pipeline-step" type="button" data-panel="recommend-params" aria-controls="detail-panel" aria-pressed="false">
      <strong>02 · recommend_params</strong>
      <span>规则引擎输出参数建议，LLM 负责解释</span>
    </button>
    <button class="hotspot pipeline-step" type="button" data-panel="run-search" aria-controls="detail-panel" aria-pressed="false">
      <strong>03 · run_search</strong>
      <span>Rust 调用搜索引擎并产出结构化结果</span>
    </button>
    <button class="hotspot pipeline-step" type="button" data-panel="generate-report" aria-controls="detail-panel" aria-pressed="false">
      <strong>04 · generate_summary / export_results</strong>
      <span>输出统计摘要与文件结果，供 LLM 进一步解释</span>
    </button>
  </div>
</section>

<section id="extensibility" class="page-section reveal">
  <p class="eyebrow">Section 04</p>
  <h2>可扩展性：架构不是终点，而是平台起点</h2>
  <p>这一段要把“可以继续长大”讲清楚，让观众看到 search adapter、MCP 模块和 server 拆分都是被架构显式预留的。</p>
  <div class="section-grid extension-grid">
    <button class="hotspot extension-card" type="button" data-panel="new-adapter" aria-controls="detail-panel" aria-pressed="false">
      <strong>新增搜索引擎 adapter</strong>
      <span>MSFragger / Comet 只需要实现统一的 adapter 契约</span>
    </button>
    <button class="hotspot extension-card" type="button" data-panel="new-module" aria-controls="detail-panel" aria-pressed="false">
      <strong>新增 MCP 模块</strong>
      <span>qc、fdr、protein inference 都可以按 crate + tool 方式插入</span>
    </button>
    <button class="hotspot extension-card" type="button" data-panel="split-server" aria-controls="detail-panel" aria-pressed="false">
      <strong>未来拆成更多 MCP Server</strong>
      <span>某个领域需要独立部署时，可以从 library 提升为独立 server</span>
    </button>
  </div>
</section>
```

- [ ] **Step 3: Run the workflow/extensibility presence check again**

Run:

```bash
rg -n 'read_spectra|recommend_params|run_search|generate_summary|data-panel="new-adapter"|data-panel="new-module"|data-panel="split-server"' docs/architecture.html
```

Expected: matches for all four workflow steps and three extensibility cards.

- [ ] **Step 4: Run a browser smoke check for the second half of the page**

Run:

```bash
python -m http.server 8000 >/tmp/protein-copilot-architecture-html.log 2>&1 & echo $!
```

Expected: prints a numeric PID. Then open `http://localhost:8000/docs/architecture.html` and verify:

- the pipeline reads as a single story from input to report
- the extensibility cards clearly communicate adapter growth, module growth, and server splitting
- the page still feels presentation-first instead of reference-doc dense

After verification, stop the server with:

```bash
kill <PID>
```

- [ ] **Step 5: Commit**

```bash
git add docs/architecture.html
git commit -m "docs: add architecture workflow and extensibility sections

Add the pipeline story and the future-growth section so the page ends on
the platform's adapter and server-splitting expansion path.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 4: Add restrained interactions, detail-panel behavior, and responsive polish

**Files:**
- Modify: `docs/architecture.html`
- Test: manual browser interaction check of `docs/architecture.html`

- [ ] **Step 1: Run the interaction-hook check before editing**

Run:

```bash
rg -n 'IntersectionObserver|setActivePanel|aria-pressed|const detailData = \\{|classList\\.add\\("is-visible"' docs/architecture.html
```

Expected: no matches for the JS helpers and reveal state classes.

- [ ] **Step 2: Extend the CSS so the page looks like a restrained presentation page instead of a plain document**

Add or replace these style blocks inside the existing `<style>` tag:

```css
.mcp-grid,
.architecture-layout,
.extension-grid {
  grid-template-columns: repeat(3, minmax(0, 1fr));
}

.info-card,
.takeaway-card,
.extension-card,
.pipeline-step,
.crate-card,
.layer {
  padding: 18px;
  border: 1px solid var(--border);
  border-radius: 20px;
  background: #fff;
  color: var(--ink);
  box-shadow: 0 10px 24px rgba(15, 23, 42, 0.06);
}

.stack,
.pipeline-steps {
  display: grid;
  gap: 14px;
}

.crate-grid {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 12px;
}

.hotspot {
  width: 100%;
  border: 1px solid var(--border);
  cursor: pointer;
  text-align: left;
  transition: transform 180ms ease, box-shadow 180ms ease, border-color 180ms ease;
}

.hotspot:hover,
.hotspot:focus-visible {
  transform: translateY(-2px);
  box-shadow: 0 16px 32px rgba(15, 23, 42, 0.12);
}

.hotspot.is-active {
  border-color: var(--mcp);
  box-shadow: 0 0 0 3px rgba(8, 145, 178, 0.16);
}

.llm-layer { background: linear-gradient(135deg, rgba(124, 58, 237, 0.08), #fff); }
.mcp-layer { background: linear-gradient(135deg, rgba(8, 145, 178, 0.1), #fff); }
.ext-layer { background: linear-gradient(135deg, rgba(234, 88, 12, 0.08), #fff); }

.reveal {
  opacity: 0;
  transform: translateY(24px);
  transition: opacity 320ms ease, transform 320ms ease;
}

.reveal.is-visible {
  opacity: 1;
  transform: translateY(0);
}

@media (max-width: 960px) {
  .mcp-grid,
  .architecture-layout,
  .extension-grid,
  .crate-grid {
    grid-template-columns: 1fr;
  }
}
```

- [ ] **Step 3: Replace the empty `detailData` stub with real panel data and interaction logic**

Replace:

```html
<script>
  const detailData = {};
</script>
```

with:

```html
<script>
  const detailData = {
    "llm-layer": {
      title: "LLM 编排层",
      body: "负责理解用户意图、组织调用顺序、把 Rust 返回的结构化结果解释成人能读懂的结论。"
    },
    "mcp-server": {
      title: "ProteinCopilot MCP Server",
      body: "统一暴露 tool 接口，连接 Agent/Prompt 与下层 Rust crate，并保持交互契约结构化。"
    },
    "core-crate": {
      title: "core",
      body: "定义共享领域模型、错误类型和 adapter trait，是整个 workspace 的共享数据基础。"
    },
    "spectrum-io-crate": {
      title: "spectrum-io",
      body: "负责 mzML / mgf 的读取和摘要提取，是所有后续流程的输入入口。"
    },
    "param-recommend-crate": {
      title: "param-recommend",
      body: "规则引擎根据谱图特征推荐搜索参数，保持推荐逻辑确定性。"
    },
    "search-engine-crate": {
      title: "search-engine",
      body: "封装搜索引擎适配层与统一结果结构，是执行搜索的核心模块。"
    },
    "report-crate": {
      title: "report",
      body: "把搜索结果变成结构化摘要和导出文件，供 LLM 和用户继续消费。"
    },
    "dia-extraction-crate": {
      title: "dia-extraction",
      body: "为 DIA 工作流提取候选前体离子，把宽窗口数据接到后续搜索流程。"
    },
    "search-adapter": {
      title: "搜索引擎 / Adapter",
      body: "当前有 SimpleSearch，未来可以继续接入 pFind、MSFragger、Comet。"
    },
    "read-spectra": {
      title: "read_spectra",
      body: "先读取文件并总结数据特征，为后续推荐和搜索准备结构化输入。"
    },
    "recommend-params": {
      title: "recommend_params",
      body: "Rust 给出可重复的参数建议，LLM 再把建议解释成用户能理解的话。"
    },
    "run-search": {
      title: "run_search",
      body: "Rust 调用搜索引擎并返回统一的结果结构，避免 LLM 直接碰评分与数值计算。"
    },
    "generate-report": {
      title: "generate_summary / export_results",
      body: "输出摘要和文件结果，让后续解读建立在结构化数据之上。"
    },
    "new-adapter": {
      title: "新增 adapter",
      body: "通过统一的 SearchEngineAdapter 契约接入更多搜索引擎，而不改上层交互。"
    },
    "new-module": {
      title: "新增 MCP 模块",
      body: "qc、fdr 或 protein inference 都可以作为新 crate 和新 tool 独立加入。"
    },
    "split-server": {
      title: "拆分更多 MCP Server",
      body: "当某个领域需要独立部署时，可以把 library 抽成独立 server，而不推翻现有分层。"
    }
  };

  const titleNode = document.getElementById("detail-title");
  const bodyNode = document.getElementById("detail-body");
  const hotspots = document.querySelectorAll(".hotspot");

  function setActivePanel(panelId) {
    const detail = detailData[panelId];
    if (!detail) return;

    titleNode.textContent = detail.title;
    bodyNode.textContent = detail.body;

    hotspots.forEach((node) => {
      const active = node.dataset.panel === panelId;
      node.classList.toggle("is-active", active);
      node.setAttribute("aria-pressed", active ? "true" : "false");
    });
  }

  hotspots.forEach((node) => {
    node.addEventListener("click", () => setActivePanel(node.dataset.panel));
    node.addEventListener("mouseenter", () => {
      if (window.matchMedia("(pointer:fine)").matches) {
        setActivePanel(node.dataset.panel);
      }
    });
  });

  const observer = new IntersectionObserver((entries) => {
    entries.forEach((entry) => {
      if (entry.isIntersecting) {
        entry.target.classList.add("is-visible");
      }
    });
  }, { threshold: 0.2 });

  document.querySelectorAll(".reveal").forEach((node) => observer.observe(node));
  setActivePanel("mcp-server");
</script>
```

- [ ] **Step 4: Run the interaction-hook check again**

Run:

```bash
rg -n 'IntersectionObserver|setActivePanel|aria-pressed|const detailData = \\{|classList\\.add\\("is-visible"' docs/architecture.html
```

Expected: matches for the interaction helpers and reveal logic.

- [ ] **Step 5: Run a manual interaction smoke test**

Run:

```bash
python -m http.server 8000 >/tmp/protein-copilot-architecture-html.log 2>&1 & echo $!
```

Expected: prints a numeric PID. Then open `http://localhost:8000/docs/architecture.html` and verify:

- sections fade in as you scroll, but motion stays subtle
- clicking `ProteinCopilot MCP Server`, `search-engine`, and `新增 adapter` updates the detail panel
- the page still reads sensibly if you never click anything
- mobile-width layout stacks cleanly with no overlapping cards

After verification, stop the server with:

```bash
kill <PID>
```

- [ ] **Step 6: Commit**

```bash
git add docs/architecture.html
git commit -m "docs: add architecture interactions and presentation polish

Add restrained hotspot behavior, scroll reveals, detail-panel content, and
responsive styling to the architecture presentation page.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```

---

### Task 5: Link the new page into the existing docs and run the final smoke check

**Files:**
- Modify: `README.md:96-97`
- Modify: `docs/architecture.md:1-6`
- Test: manual browser smoke check of `docs/architecture.html`

- [ ] **Step 1: Run the link-presence check before editing**

Run:

```bash
rg -n 'architecture.html' README.md docs/architecture.md
```

Expected: no matches.

- [ ] **Step 2: Add the new presentation-page link to `README.md`**

Change the architecture reference block near the end of `README.md` to:

```md
详细计划：`tasks/001-mvp-proteomics-search-platform.md`
架构设计：`docs/architecture.md`
架构演示：`docs/architecture.html`
```

- [ ] **Step 3: Add a short note near the top of `docs/architecture.md`**

Insert this note below the introductory block quote:

```md
> 演示版可视化页面见：`docs/architecture.html`
>  
> `docs/architecture.md` 保留为完整的架构说明和决策记录。
```

- [ ] **Step 4: Run the link-presence check again**

Run:

```bash
rg -n 'architecture.html' README.md docs/architecture.md
```

Expected: one match in `README.md` and one match in `docs/architecture.md`.

- [ ] **Step 5: Run the final browser smoke check**

Run:

```bash
python -m http.server 8000 >/tmp/protein-copilot-architecture-html.log 2>&1 & echo $!
```

Expected: prints a numeric PID. Then open `http://localhost:8000/docs/architecture.html` and verify:

- the page loads with no missing styles or scripts
- the five-section narrative order matches the approved spec
- the MCP explanation is brief and understandable
- the layered architecture and extensibility sections land the two core messages

After verification, stop the server with:

```bash
kill <PID>
```

- [ ] **Step 6: Commit**

```bash
git add docs/architecture.html README.md docs/architecture.md
git commit -m "docs: link architecture presentation page

Expose the new architecture showcase from the existing README and detailed
architecture documentation.

Co-authored-by: Copilot <223556219+Copilot@users.noreply.github.com>"
```
