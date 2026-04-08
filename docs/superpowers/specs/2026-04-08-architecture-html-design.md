# Architecture HTML Showcase Design

## Problem

`docs/architecture.md` already captures the project architecture in detail, but it is optimized for
reading, not for live demos or presentation walkthroughs. The project needs a presentation-friendly
`docs/architecture.html` that quickly communicates the system shape, the Rust/LLM boundary, and the
platform's extensibility.

## Goals

- Create a presentation-first architecture page at `docs/architecture.html`
- Make the two primary messages obvious:
  - Rust performs deterministic scientific computation
  - LLM handles orchestration, intent understanding, and explanation
- Show that the architecture is designed to grow:
  - new search engine adapters can be added
  - new MCP modules can be introduced
  - individual domains can later be split into separate MCP servers
- Include a brief "what is MCP?" explanation near the top of the page
- Support live walkthroughs with light interaction, while staying readable as a static page
- Keep the implementation as a single self-contained HTML file with no build step

## Non-Goals

- Replacing `docs/architecture.md` as the canonical detailed architecture document
- Building a general-purpose docs site or SPA
- Mirroring every ADR, code sample, or tool definition from the Markdown docs
- Adding external assets, frameworks, or a front-end build pipeline

## Audience and Use Case

Primary audience:

- demo and presentation viewers
- collaborators who need a quick architectural overview

Primary usage mode:

- the page is shown during live walkthroughs and screen sharing

Secondary usage mode:

- the page can be opened directly by a developer for a fast self-guided overview

## Content Sources

The HTML page should synthesize current architecture information from:

- `docs/architecture.md`
- `README.md`
- `docs/mcp-tools.md`

The page content should stay aligned with the current crate layout and architecture principles:

- `core`
- `spectrum-io`
- `param-recommend`
- `search-engine`
- `report`
- `dia-extraction`
- `mcp-server`

## Proposed Page Structure

The approved narrative flow is:

1. **Hero**
   - Title: ProteinCopilot architecture
   - Short positioning statement for the platform
   - Immediate emphasis on the Rust/LLM responsibility split
2. **MCP quick explainer**
   - A compact explanation of MCP as the protocol bridge between the LLM client and Rust tools
   - Keep this brief and presentation-friendly, not protocol-spec heavy
3. **Layered architecture**
   - Show the stack from user/client → LLM orchestration → MCP server → Rust crates → search engine
   - Highlight the deterministic boundary between orchestration and computation
4. **End-to-end pipeline**
   - Show the main product flow: read spectra → recommend params → run search → generate/export report
   - Present this as the "from data to result" story
5. **Extensibility paths**
   - Show how the architecture evolves through new search adapters, new MCP modules, and future server splitting

This sequence is intentionally linear and presentation-friendly: explain the platform, define MCP,
show the architecture, show the workflow, then end on growth potential.

## Interaction Model

The page should use **selective interaction**, not a fully exploratory application.

Approved interaction principles:

- Scroll-based section reveal is allowed, but motion stays subtle
- Key diagrams support hover/click highlighting
- Clicking a crate or architecture node can reveal a short responsibility summary
- Extensibility nodes can reveal future evolution examples
- Interaction is additive only: the page must remain understandable without any clicking

Explicitly avoid:

- drag/pan canvases
- dense graph-navigation UIs
- heavy animation that competes with live narration
- framework-driven complexity for simple reveal behavior

## Visual Style

The approved direction is a **restrained presentation style**:

- dark hero section for impact
- lighter content sections for readability
- strong but controlled contrast for projector use
- consistent layer colors to distinguish:
  - user / client / LLM
  - MCP server
  - Rust crates
  - external search engine / future extension nodes

The page should feel more like a polished architecture walkthrough than a marketing landing page or
a plain technical document.

## Technical Approach

Implementation boundary approved by the user:

- output is a single file: `docs/architecture.html`
- use inline CSS and a small amount of vanilla JavaScript
- no framework, bundler, or build step
- no required companion assets

Recommended implementation details:

- use semantic HTML sections for the scroll narrative
- use inline SVG or simple HTML/CSS blocks for diagrams
- use `data-*` attributes plus small JS helpers for hotspot/highlight behavior
- keep all core information in the DOM by default so the page still works as static content
- ensure the file opens correctly both from a browser directly and from a simple static file server

## Content Mapping

The HTML page should visually emphasize these messages from the existing architecture docs:

1. **Rust vs LLM boundary**
   - Rust owns deterministic parsing, validation, search, aggregation, and export
   - LLM owns intent understanding, workflow orchestration, and natural-language explanation
2. **Unified MCP entry point**
   - the MCP server is the bridge between LLM clients and deterministic Rust capabilities
3. **Workspace modularity**
   - crate boundaries are explicit and intentionally layered
4. **End-to-end workflow**
   - the platform covers the path from raw spectra to interpretable results
5. **Future growth**
   - adapters and MCP module boundaries make evolution straightforward

## Acceptance Criteria

- A presenter can use the page without also opening `docs/architecture.md`
- A first-time viewer can understand what MCP does in this project within one short section
- The layered diagram clearly shows where LLM responsibilities stop and Rust responsibilities begin
- The pipeline section clearly communicates the product flow from input data to output report
- The extensibility section clearly shows at least:
  - adding more search engine adapters
  - adding more MCP modules
  - splitting a domain into a dedicated MCP server later
- The page remains readable even if the viewer does not interact with hotspots
- The page requires no build step and can live as a durable docs artifact

## Risks and Mitigations

- **Risk:** Too much detail turns the page back into a dense doc
  - **Mitigation:** Keep each section anchored around one presentation takeaway
- **Risk:** Interaction adds fragility without adding value
  - **Mitigation:** Limit interaction to hotspot highlighting and short detail reveals
- **Risk:** The MCP explanation becomes too abstract
  - **Mitigation:** Explain MCP through its role in this project, not as a generic protocol lecture

## Implementation Notes

- `docs/architecture.md` remains the detailed reference and source of truth
- `docs/architecture.html` is the visual presentation layer derived from that reference
- The page should prefer durable, hand-authored content over generated markup so it stays easy to
  revise during future architecture updates
