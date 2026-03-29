//! Search engine registry — manages available search engine adapters.

use std::collections::HashMap;

use protein_copilot_core::engine::{EngineInfo, HealthStatus, SearchEngineAdapter};

/// Registry of available search engine adapters.
///
/// The MCP tool layer uses this to discover and select engines.
pub struct EngineRegistry {
    engines: HashMap<String, Box<dyn SearchEngineAdapter>>,
}

impl EngineRegistry {
    /// Creates a new empty registry.
    pub fn new() -> Self {
        Self {
            engines: HashMap::new(),
        }
    }

    /// Registers a search engine adapter.
    ///
    /// The engine is identified by its `engine_info().name`.
    pub fn register(&mut self, adapter: Box<dyn SearchEngineAdapter>) {
        let name = adapter.engine_info().name.clone();
        self.engines.insert(name, adapter);
    }

    /// Retrieves an adapter by engine name.
    pub fn get(&self, name: &str) -> Option<&dyn SearchEngineAdapter> {
        self.engines.get(name).map(|a| a.as_ref())
    }

    /// Lists all registered engines and their info.
    pub fn list_available(&self) -> Vec<EngineInfo> {
        self.engines.values().map(|a| a.engine_info()).collect()
    }

    /// Checks health of all registered engines.
    pub async fn health_check_all(
        &self,
    ) -> Vec<(
        EngineInfo,
        Result<HealthStatus, protein_copilot_core::error::CoreError>,
    )> {
        let mut results = Vec::new();
        for adapter in self.engines.values() {
            let info = adapter.engine_info();
            let status = adapter.health_check().await;
            results.push((info, status));
        }
        results
    }

    /// Returns the number of registered engines.
    pub fn len(&self) -> usize {
        self.engines.len()
    }

    /// Returns true if no engines are registered.
    pub fn is_empty(&self) -> bool {
        self.engines.is_empty()
    }
}

impl Default for EngineRegistry {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_registry() {
        let registry = EngineRegistry::new();
        assert!(registry.is_empty());
        assert_eq!(registry.len(), 0);
        assert!(registry.list_available().is_empty());
        assert!(registry.get("pFind").is_none());
    }
}
