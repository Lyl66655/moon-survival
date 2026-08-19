///|
/// Reproducibility configuration for a survival analysis.
pub(all) struct SurvivalReproducibilityConfig {
  seed : UInt
  package_version : String
  toolchain : String
  locale : String
  reproducibility_row : Bool
} derive(Eq, Debug)

///|
/// A deterministic analysis identifier.
pub(all) struct SurvivalAnalysisIdentifier {
  project : String
  analysis : String
  data_checksum : String
  identifier : String
} derive(Eq, Debug)

///|
/// A reproducibility manifest.
pub(all) struct SurvivalReproducibilityManifest {
  config : SurvivalReproducibilityConfig
  identifier : SurvivalAnalysisIdentifier
  source_rows : Int
  source_events : Int
  manifest_checksum : String
} derive(Eq, Debug)

///|
/// Create a reproducibility configuration.
pub fn survival_reproducibility_config(
  seed : UInt,
  package_version : String,
) -> SurvivalReproducibilityConfig {
  {
    seed,
    package_version,
    toolchain: "moon-stable",
    locale: "C",
    reproducibility_row: true,
  }
}

///|
/// Derive an analysis identifier from data and labels.
pub fn survival_analysis_identifier(
  project : String,
  analysis : String,
  data : SurvivalData,
) -> SurvivalAnalysisIdentifier {
  let checksum = survival_export_checksum(survival_data_to_csv(data))
  {
    project,
    analysis,
    data_checksum: checksum,
    identifier: project + ":" + analysis + ":" + checksum,
  }
}

///|
/// Build a reproducibility manifest.
pub fn survival_reproducibility_manifest(
  config : SurvivalReproducibilityConfig,
  project : String,
  analysis : String,
  data : SurvivalData,
) -> SurvivalReproducibilityManifest {
  let identifier = survival_analysis_identifier(project, analysis, data)
  let content = identifier.identifier +
    config.package_version +
    config.toolchain +
    config.locale
  {
    config,
    identifier,
    source_rows: data.length(),
    source_events: data.event_count(),
    manifest_checksum: survival_export_checksum(content),
  }
}

///|
/// Derive a replicate seed without global mutable state.
pub fn survival_replicate_seed(seed : UInt, replicate : Int) -> UInt {
  seed + replicate.reinterpret_as_uint() * 2654435761U
}

///|
/// Derive a fold seed.
pub fn survival_fold_seed(seed : UInt, fold : Int) -> UInt {
  survival_replicate_seed(seed, fold + 1000)
}

///|
/// Derive a bootstrap seed.
pub fn survival_bootstrap_seed(seed : UInt, replicate : Int) -> UInt {
  survival_replicate_seed(seed, replicate + 100000)
}

///|
/// Create an integer fingerprint of a vector.
pub fn survival_integer_fingerprint(values : Array[Int]) -> String {
  let mut value = 7
  for item in values {
    value = (value * 31 + item) % 1000000007
  }
  value.to_string()
}

///|
/// Create a floating-point fingerprint by decimal serialization.
pub fn survival_numeric_fingerprint(values : Array[Double]) -> String {
  let builder = StringBuilder::new()
  for value in values {
    builder.write_string(value.to_string())
    builder.write_char('|')
  }
  survival_export_checksum(builder.to_string())
}

///|
/// Check whether a configuration is complete.
pub fn survival_reproducibility_is_valid(
  config : SurvivalReproducibilityConfig,
) -> Bool {
  config.package_version != "" && config.toolchain != "" && config.locale != ""
}

///|
/// Check a manifest's consistency.
pub fn survival_manifest_is_valid(
  manifest : SurvivalReproducibilityManifest,
) -> Bool {
  survival_reproducibility_is_valid(manifest.config) &&
  manifest.identifier.identifier != "" &&
  manifest.source_rows >= manifest.source_events &&
  manifest.source_events >= 0 &&
  manifest.manifest_checksum != ""
}

///|
/// Serialize a reproducibility manifest.
pub fn survival_reproducibility_to_csv(
  manifest : SurvivalReproducibilityManifest,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("field,value\n")
  builder.write_string("seed,")
  builder.write_string(manifest.config.seed.to_string())
  builder.write_char('\n')
  builder.write_string("package_version,")
  builder.write_string(manifest.config.package_version)
  builder.write_char('\n')
  builder.write_string("toolchain,")
  builder.write_string(manifest.config.toolchain)
  builder.write_char('\n')
  builder.write_string("locale,")
  builder.write_string(manifest.config.locale)
  builder.write_char('\n')
  builder.write_string("identifier,")
  builder.write_string(manifest.identifier.identifier)
  builder.write_char('\n')
  builder.write_string("source_rows,")
  builder.write_string(manifest.source_rows.to_string())
  builder.write_char('\n')
  builder.write_string("source_events,")
  builder.write_string(manifest.source_events.to_string())
  builder.write_char('\n')
  builder.write_string("manifest_checksum,")
  builder.write_string(manifest.manifest_checksum)
  builder.write_char('\n')
  builder.to_string()
}
