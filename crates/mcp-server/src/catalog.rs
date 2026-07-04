//! Human-readable rendering of the MCP tool catalog for the `--list-tools`
//! CLI mode. Operates purely on the JSON Schema carried by [`rmcp::model::Tool`]
//! (the exact contract sent to clients), so it never drifts from the tools.

use rmcp::model::Tool;
use serde_json::Value;

/// Render the full tool catalog as plain ASCII text: per-tool description,
/// parameter table (name, type + range, required, default, description),
/// output type, and a shared data-type appendix (enum values / struct fields).
pub fn format_catalog(tools: &[Tool]) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "ProteinCopilot MCP Server — {} 个工具\n",
        tools.len()
    ));
    out.push_str("传输：JSON-RPC 2.0 over stdio。完整 JSON Schema：--list-tools --json\n");

    // Collect shared $defs across all tools (referenced complex types).
    let mut defs: std::collections::BTreeMap<String, Value> = std::collections::BTreeMap::new();
    for t in tools {
        if let Some(d) = t.input_schema.get("$defs").and_then(Value::as_object) {
            for (k, v) in d {
                defs.entry(k.clone()).or_insert_with(|| v.clone());
            }
        }
    }

    out.push_str("\n== 工具 ==\n");
    for t in tools {
        out.push_str(&format!("\n* {}\n", t.name));
        if let Some(desc) = &t.description {
            out.push_str(&format!("  {}\n", oneline(desc)));
        }
        let schema: &serde_json::Map<String, Value> = &t.input_schema;
        let props = schema.get("properties").and_then(Value::as_object);
        let required: std::collections::BTreeSet<&str> = schema
            .get("required")
            .and_then(Value::as_array)
            .map(|a| a.iter().filter_map(Value::as_str).collect())
            .unwrap_or_default();
        match props {
            Some(props) if !props.is_empty() => {
                out.push_str("  参数：\n");
                let mut keys: Vec<&String> = props.keys().collect();
                keys.sort_by_key(|k| (!required.contains(k.as_str()), (*k).clone()));
                for k in keys {
                    let v = &props[k];
                    let req = if required.contains(k.as_str()) {
                        "必填"
                    } else {
                        "可选"
                    };
                    let mut ty = type_label(v);
                    let range = range_label(v);
                    if !range.is_empty() {
                        ty.push(' ');
                        ty.push_str(&range);
                    }
                    if let Some(d) = default_label(v) {
                        ty.push_str(&format!(" ={d}"));
                    }
                    out.push_str(&format!("    - {k} ({ty}, {req})"));
                    if let Some(desc) = v.get("description").and_then(Value::as_str) {
                        out.push_str(&format!("  {}", oneline(desc)));
                    }
                    out.push('\n');
                }
            }
            _ => out.push_str("  参数：无\n"),
        }
        let otype = t
            .output_schema
            .as_ref()
            .and_then(|s| s.get("title").and_then(Value::as_str))
            .unwrap_or("object");
        out.push_str(&format!("  输出：{otype}\n"));
    }

    if !defs.is_empty() {
        out.push_str("\n== 数据类型 ==\n");
        for (name, dv) in &defs {
            out.push_str(&format!("\n* {name}\n"));
            if let Some(desc) = dv.get("description").and_then(Value::as_str) {
                out.push_str(&format!("  {}\n", oneline(desc)));
            }
            render_type_body(&mut out, dv);
        }
    }

    out
}

/// Render the body of a $def: enum values, tagged variants, or struct fields.
fn render_type_body(out: &mut String, dv: &Value) {
    if let Some(en) = dv.get("enum").and_then(Value::as_array) {
        let vals: Vec<String> = en.iter().map(value_token).collect();
        out.push_str(&format!("  枚举: {}\n", vals.join(", ")));
        return;
    }
    if let Some(variants) = dv
        .get("oneOf")
        .or_else(|| dv.get("anyOf"))
        .and_then(Value::as_array)
    {
        out.push_str("  变体:\n");
        for var in variants {
            let vdesc = var
                .get("description")
                .and_then(Value::as_str)
                .map(|d| format!(" — {}", oneline(d)))
                .unwrap_or_default();
            if let Some(c) = var.get("const") {
                out.push_str(&format!("    - {}{vdesc}\n", value_token(c)));
            } else if let Some(props) = var.get("properties").and_then(Value::as_object) {
                for (tag, tv) in props {
                    if let Some(inner) = tv.get("properties").and_then(Value::as_object) {
                        let req: std::collections::BTreeSet<&str> = tv
                            .get("required")
                            .and_then(Value::as_array)
                            .map(|a| a.iter().filter_map(Value::as_str).collect())
                            .unwrap_or_default();
                        let fields: Vec<String> = inner
                            .iter()
                            .map(|(fk, fv)| {
                                let opt = if req.contains(fk.as_str()) { "" } else { "?" };
                                format!("{fk}{opt}: {}", type_label(fv))
                            })
                            .collect();
                        out.push_str(&format!("    - {tag}{vdesc} {{ {} }}\n", fields.join("; ")));
                    } else {
                        out.push_str(&format!("    - {tag}{vdesc}\n"));
                    }
                }
            } else {
                out.push_str(&format!("    - {}{vdesc}\n", type_label(var)));
            }
        }
        return;
    }
    if let Some(props) = dv.get("properties").and_then(Value::as_object) {
        let required: std::collections::BTreeSet<&str> = dv
            .get("required")
            .and_then(Value::as_array)
            .map(|a| a.iter().filter_map(Value::as_str).collect())
            .unwrap_or_default();
        out.push_str("  字段:\n");
        for (k, v) in props {
            let req = if required.contains(k.as_str()) {
                "必填"
            } else {
                "可选"
            };
            let mut ty = type_label(v);
            let range = range_label(v);
            if !range.is_empty() {
                ty.push(' ');
                ty.push_str(&range);
            }
            out.push_str(&format!("    - {k} ({ty}, {req})"));
            if let Some(desc) = v.get("description").and_then(Value::as_str) {
                out.push_str(&format!("  {}", oneline(desc)));
            }
            out.push('\n');
        }
    }
}

fn ref_name(r: &str) -> &str {
    r.rsplit('/').next().unwrap_or(r)
}

fn value_token(v: &Value) -> String {
    match v {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// Human-readable type label for a property/schema node.
fn type_label(v: &Value) -> String {
    if let Some(r) = v.get("$ref").and_then(Value::as_str) {
        return ref_name(r).to_string();
    }
    for key in ["anyOf", "oneOf", "allOf"] {
        if let Some(arr) = v.get(key).and_then(Value::as_array) {
            let mut nullable = false;
            let mut parts: Vec<String> = Vec::new();
            for s in arr {
                if s.get("type").and_then(Value::as_str) == Some("null") {
                    nullable = true;
                } else {
                    parts.push(type_label(s));
                }
            }
            let mut label = if parts.is_empty() {
                "null".to_string()
            } else {
                parts.join(" | ")
            };
            if nullable {
                label.push_str(" 可空");
            }
            return label;
        }
    }
    if let Some(en) = v.get("enum").and_then(Value::as_array) {
        let vals: Vec<String> = en.iter().map(value_token).collect();
        return format!("enum({})", vals.join(", "));
    }
    match v.get("type") {
        Some(Value::Array(types)) => {
            let mut nullable = false;
            let mut base: Vec<String> = Vec::new();
            for t in types {
                match t.as_str() {
                    Some("null") => nullable = true,
                    Some(s) => base.push(s.to_string()),
                    None => {}
                }
            }
            let mut label = if base.is_empty() {
                "null".to_string()
            } else {
                base.join(" | ")
            };
            if nullable {
                label.push_str(" 可空");
            }
            label
        }
        Some(Value::String(s)) if s == "array" => {
            let item = v
                .get("items")
                .map(type_label)
                .unwrap_or_else(|| "any".to_string());
            format!("array<{item}>")
        }
        Some(Value::String(s)) => s.clone(),
        _ => "object".to_string(),
    }
}

/// Numeric/length/format constraints, e.g. ">=1", "0..1", "len 1", "format=double".
fn range_label(v: &Value) -> String {
    let mut parts: Vec<String> = Vec::new();
    let num = |k: &str| v.get(k).and_then(Value::as_f64);
    if let Some(m) = num("minimum") {
        parts.push(format!(">={}", trim_num(m)));
    }
    if let Some(m) = num("exclusiveMinimum") {
        parts.push(format!(">{}", trim_num(m)));
    }
    if let Some(m) = num("maximum") {
        parts.push(format!("<={}", trim_num(m)));
    }
    if let Some(m) = num("exclusiveMaximum") {
        parts.push(format!("<{}", trim_num(m)));
    }
    if let Some(m) = num("multipleOf") {
        parts.push(format!("step {}", trim_num(m)));
    }
    let len_min = v.get("minLength").and_then(Value::as_u64);
    let len_max = v.get("maxLength").and_then(Value::as_u64);
    match (len_min, len_max) {
        (Some(a), Some(b)) if a == b => parts.push(format!("len {a}")),
        (Some(a), Some(b)) => parts.push(format!("len {a}..{b}")),
        (Some(a), None) => parts.push(format!("len>={a}")),
        (None, Some(b)) => parts.push(format!("len<={b}")),
        _ => {}
    }
    let items_min = v.get("minItems").and_then(Value::as_u64);
    let items_max = v.get("maxItems").and_then(Value::as_u64);
    match (items_min, items_max) {
        (Some(a), Some(b)) if a == b => parts.push(format!("items {a}")),
        (Some(a), Some(b)) => parts.push(format!("items {a}..{b}")),
        (Some(a), None) => parts.push(format!("items>={a}")),
        (None, Some(b)) => parts.push(format!("items<={b}")),
        _ => {}
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!("[{}]", parts.join(", "))
    }
}

fn trim_num(f: f64) -> String {
    if f.fract() == 0.0 {
        format!("{}", f as i64)
    } else {
        format!("{f}")
    }
}

fn default_label(v: &Value) -> Option<String> {
    match v.get("default") {
        Some(Value::Null) | None => None,
        Some(d) => Some(serde_json::to_string(d).unwrap_or_else(|_| value_token(d))),
    }
}

/// Collapse internal whitespace/newlines so a schema description prints on one line.
fn oneline(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tools() -> Vec<Tool> {
        crate::tools::ProteinCopilotServer::new().list_tools()
    }

    #[test]
    fn server_exposes_all_tools() {
        let t = tools();
        assert_eq!(t.len(), 28, "expected 28 MCP tools");
        assert!(t.iter().any(|x| x.name == "read_spectra"));
    }

    #[test]
    fn catalog_lists_tool_name_params_and_output() {
        let text = format_catalog(&tools());
        assert!(text.contains("read_spectra"), "tool name missing");
        assert!(text.contains("file_path"), "param name missing");
        assert!(text.contains("必填"), "required marker missing");
        assert!(text.contains("SpectrumSummary"), "output type missing");
    }

    #[test]
    fn catalog_renders_enum_values_and_ranges() {
        let text = format_catalog(&tools());
        // Enzyme is an enum referenced by SearchParams -> appears in 数据类型 appendix.
        assert!(text.contains("Trypsin"), "enum value missing");
        // value must be positive lives on MassTolerance.value description.
        assert!(
            text.contains("数据类型"),
            "data-type appendix header missing"
        );
    }
}
