///|
/// A registry entry for a fitted survival model.
pub(all) struct SurvivalModelRegistryEntry {
  name : String
  family : String
  sample_size : Int
  events : Int
  score : Double
  status : String
  registry_entry : Bool
} derive(Eq, Debug)

///|
/// A model registry used to select reproducible production candidates.
pub(all) struct SurvivalModelRegistry {
  entries : Array[SurvivalModelRegistryEntry]
  registry_version : String
  registry_row : Bool
} derive(Eq, Debug)

///|
/// Create an empty registry.
pub fn survival_model_registry() -> SurvivalModelRegistry {
  { entries: [], registry_version: "1", registry_row: true }
}

///|
/// Register a model entry.
pub fn SurvivalModelRegistry::register(
  self : SurvivalModelRegistry,
  entry : SurvivalModelRegistryEntry,
) -> SurvivalModelRegistry {
  let entries : Array[SurvivalModelRegistryEntry] = []
  for value in self.entries {
    entries.push(value)
  }
  entries.push(entry)
  { entries, registry_version: self.registry_version, registry_row: true }
}

///|
/// Register a parametric fit by AIC.
pub fn SurvivalModelRegistry::register_parametric(
  self : SurvivalModelRegistry,
  name : String,
  fit : ParametricFit,
) -> SurvivalModelRegistry {
  self.register({
    name,
    family: parametric_family_name(fit.family),
    sample_size: fit.sample_size,
    events: fit.event_count,
    score: fit.aic,
    status: if fit.converged {
      "converged"
    } else {
      "not_converged"
    },
    registry_entry: true,
  })
}

///|
/// Return the best entry by score.
pub fn SurvivalModelRegistry::best(
  self : SurvivalModelRegistry,
) -> SurvivalModelRegistryEntry? {
  let mut result : SurvivalModelRegistryEntry? = None
  for entry in self.entries {
    match result {
      None => result = Some(entry)
      Some(current) => if entry.score < current.score { result = Some(entry) }
    }
  }
  result
}

///|
/// Find a named model.
pub fn SurvivalModelRegistry::find(
  self : SurvivalModelRegistry,
  name : String,
) -> SurvivalModelRegistryEntry? {
  for entry in self.entries {
    if entry.name == name {
      return Some(entry)
    }
  }
  None
}

///|
/// Return converged entries.
pub fn SurvivalModelRegistry::converged(
  self : SurvivalModelRegistry,
) -> Array[SurvivalModelRegistryEntry] {
  let result : Array[SurvivalModelRegistryEntry] = []
  for entry in self.entries {
    if entry.status == "converged" {
      result.push(entry)
    }
  }
  result
}

///|
/// Export registry entries.
pub fn survival_model_registry_to_csv(
  registry : SurvivalModelRegistry,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("name,family,sample_size,events,score,status\n")
  for entry in registry.entries {
    builder.write_string(entry.name)
    builder.write_char(',')
    builder.write_string(entry.family)
    builder.write_char(',')
    builder.write_string(entry.sample_size.to_string())
    builder.write_char(',')
    builder.write_string(entry.events.to_string())
    builder.write_char(',')
    builder.write_string(entry.score.to_string())
    builder.write_char(',')
    builder.write_string(entry.status)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate registry invariants.
pub fn survival_model_registry_is_valid(
  registry : SurvivalModelRegistry,
) -> Bool {
  if registry.registry_version == "" {
    return false
  }
  for entry in registry.entries {
    if entry.name == "" ||
      entry.family == "" ||
      entry.sample_size < 0 ||
      entry.events < 0 ||
      entry.score.is_nan() {
      return false
    }
  }
  true
}
