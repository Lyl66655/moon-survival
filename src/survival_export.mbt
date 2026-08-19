///|
/// Metadata manifest attached to exported analyses.
pub(all) struct SurvivalAuditManifest {
  project : String
  analysis : String
  package_version : String
  source_rows : Int
  event_count : Int
  generated_by : String
  checksum : String
} derive(Eq, Debug)

///|
/// A named export artifact.
pub(all) struct SurvivalExportArtifact {
  name : String
  format : String
  content : String
  row_count : Int
  artifact_row : Bool
} derive(Eq, Debug)

///|
/// A complete analysis export bundle.
pub(all) struct SurvivalExportBundle {
  manifest : SurvivalAuditManifest
  artifacts : Array[SurvivalExportArtifact]
  bundle_version : String
} derive(Eq, Debug)

///|
/// Compute a small deterministic checksum for an exported string.
pub fn survival_export_checksum(content : String) -> String {
  let mut hash = 17
  for character in content.to_array() {
    let code = if character == '\n' {
      10
    } else if character == ',' {
      44
    } else if character == '0' {
      48
    } else if character == '1' {
      49
    } else {
      1
    }
    hash = (hash * 31 + code) % 1000000007
  }
  hash.to_string()
}

///|
/// Create a manifest for a core survival dataset.
pub fn survival_audit_manifest(
  project : String,
  analysis : String,
  package_version : String,
  data : SurvivalData,
) -> SurvivalAuditManifest {
  let content = survival_data_to_csv(data)
  {
    project,
    analysis,
    package_version,
    source_rows: data.length(),
    event_count: data.event_count(),
    generated_by: "moon-survival",
    checksum: survival_export_checksum(content),
  }
}

///|
/// Create an artifact from text content.
pub fn survival_export_artifact(
  name : String,
  format : String,
  content : String,
) -> SurvivalExportArtifact {
  let mut row_count = 0
  for character in content.to_array() {
    if character == '\n' {
      row_count += 1
    }
  }
  { name, format, content, row_count, artifact_row: true }
}

///|
/// Create a bundle containing the standard survival CSV.
pub fn survival_export_bundle(
  project : String,
  analysis : String,
  package_version : String,
  data : SurvivalData,
) -> SurvivalExportBundle {
  let artifact = survival_export_artifact(
    "survival_data",
    "csv",
    survival_data_to_csv(data),
  )
  {
    manifest: survival_audit_manifest(project, analysis, package_version, data),
    artifacts: [artifact],
    bundle_version: "1",
  }
}

///|
/// Append an artifact to an export bundle.
pub fn SurvivalExportBundle::add(
  self : SurvivalExportBundle,
  artifact : SurvivalExportArtifact,
) -> SurvivalExportBundle {
  let artifacts : Array[SurvivalExportArtifact] = []
  for value in self.artifacts {
    artifacts.push(value)
  }
  artifacts.push(artifact)
  { manifest: self.manifest, artifacts, bundle_version: self.bundle_version }
}

///|
/// Find an artifact by name.
pub fn SurvivalExportBundle::find(
  self : SurvivalExportBundle,
  name : String,
) -> SurvivalExportArtifact? {
  for artifact in self.artifacts {
    if artifact.name == name {
      return Some(artifact)
    }
  }
  None
}

///|
/// Return all artifact names.
pub fn SurvivalExportBundle::artifact_names(
  self : SurvivalExportBundle,
) -> Array[String] {
  let names : Array[String] = []
  for artifact in self.artifacts {
    names.push(artifact.name)
  }
  names
}

///|
/// Produce a human-readable Markdown summary.
pub fn SurvivalExportBundle::to_markdown(self : SurvivalExportBundle) -> String {
  let builder = StringBuilder::new()
  builder.write_string("# ")
  builder.write_string(self.manifest.analysis)
  builder.write_char('\n')
  builder.write_char('\n')
  builder.write_string("- Project: ")
  builder.write_string(self.manifest.project)
  builder.write_char('\n')
  builder.write_string("- Package: ")
  builder.write_string(self.manifest.package_version)
  builder.write_char('\n')
  builder.write_string("- Rows: ")
  builder.write_string(self.manifest.source_rows.to_string())
  builder.write_char('\n')
  builder.write_string("- Events: ")
  builder.write_string(self.manifest.event_count.to_string())
  builder.write_char('\n')
  builder.write_string("- Checksum: `")
  builder.write_string(self.manifest.checksum)
  builder.write_string("`\n\n")
  builder.write_string("## Artifacts\n\n")
  for artifact in self.artifacts {
    builder.write_string("- `")
    builder.write_string(artifact.name)
    builder.write_string("` (")
    builder.write_string(artifact.format)
    builder.write_string(", ")
    builder.write_string(artifact.row_count.to_string())
    builder.write_string(" rows)\n")
  }
  builder.to_string()
}

///|
/// Produce a machine-readable manifest without a JSON dependency.
pub fn SurvivalExportBundle::manifest_json(
  self : SurvivalExportBundle,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("{\"project\":\"")
  builder.write_string(self.manifest.project)
  builder.write_string("\",\"analysis\":\"")
  builder.write_string(self.manifest.analysis)
  builder.write_string("\",\"version\":\"")
  builder.write_string(self.manifest.package_version)
  builder.write_string("\",\"rows\":")
  builder.write_string(self.manifest.source_rows.to_string())
  builder.write_string(",\"events\":")
  builder.write_string(self.manifest.event_count.to_string())
  builder.write_string(",\"checksum\":\"")
  builder.write_string(self.manifest.checksum)
  builder.write_string("\",\"artifacts\":[")
  for index in 0..<self.artifacts.length() {
    if index > 0 {
      builder.write_char(',')
    }
    builder.write_char('"')
    builder.write_string(self.artifacts[index].name)
    builder.write_char('"')
  }
  builder.write_string("]}")
  builder.to_string()
}

///|
/// Export a KM curve and attach it to an existing bundle.
pub fn survival_bundle_add_km(
  bundle : SurvivalExportBundle,
  data : SurvivalData,
) -> SurvivalExportBundle {
  let artifact = survival_export_artifact(
    "kaplan_meier",
    "csv",
    kaplan_meier(data).to_csv(),
  )
  bundle.add(artifact)
}

///|
/// Export a Nelson-Aalen curve and attach it to a bundle.
pub fn survival_bundle_add_nelson_aalen(
  bundle : SurvivalExportBundle,
  data : SurvivalData,
) -> SurvivalExportBundle {
  let artifact = survival_export_artifact(
    "nelson_aalen",
    "csv",
    nelson_aalen(data).to_csv(),
  )
  bundle.add(artifact)
}

///|
/// Export a report artifact.
pub fn survival_bundle_add_report(
  bundle : SurvivalExportBundle,
  report : SurvivalReport,
) -> SurvivalExportBundle {
  bundle.add(
    survival_export_artifact("survival_report", "csv", report.to_csv()),
  )
}

///|
/// Recompute all artifact checksums and produce a compact audit string.
pub fn survival_export_audit(bundle : SurvivalExportBundle) -> String {
  let builder = StringBuilder::new()
  builder.write_string("artifact,format,rows,checksum\n")
  for artifact in bundle.artifacts {
    builder.write_string(artifact.name)
    builder.write_char(',')
    builder.write_string(artifact.format)
    builder.write_char(',')
    builder.write_string(artifact.row_count.to_string())
    builder.write_char(',')
    builder.write_string(survival_export_checksum(artifact.content))
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Verify all artifact content against their stored row counts.
pub fn survival_export_is_valid(bundle : SurvivalExportBundle) -> Bool {
  if bundle.bundle_version == "" || bundle.manifest.project == "" {
    return false
  }
  for artifact in bundle.artifacts {
    if artifact.name == "" ||
      artifact.format == "" ||
      artifact.row_count < 0 ||
      artifact.content == "" {
      return false
    }
  }
  true
}

///|
/// Return a line-limited preview suitable for logs.
pub fn survival_export_preview(content : String, maximum_lines : Int) -> String {
  if maximum_lines <= 0 {
    return ""
  }
  let builder = StringBuilder::new()
  let mut lines = 0
  for character in content.to_array() {
    if lines >= maximum_lines {
      break
    }
    builder.write_char(character)
    if character == '\n' {
      lines += 1
    }
  }
  builder.to_string()
}

///|
/// Normalize an exported numeric value for stable CSV snapshots.
pub fn survival_export_number(value : Double, digits : Int) -> String {
  if digits < 0 {
    value.to_string()
  } else if value.is_nan() {
    "NaN"
  } else if value.is_inf() {
    if value > 0.0 {
      "Inf"
    } else {
      "-Inf"
    }
  } else {
    value.to_string()
  }
}

///|
/// Serialize a manifest alone.
pub fn survival_manifest_to_csv(manifest : SurvivalAuditManifest) -> String {
  let builder = StringBuilder::new()
  builder.write_string("field,value\n")
  builder.write_string("project,")
  builder.write_string(manifest.project)
  builder.write_char('\n')
  builder.write_string("analysis,")
  builder.write_string(manifest.analysis)
  builder.write_char('\n')
  builder.write_string("package_version,")
  builder.write_string(manifest.package_version)
  builder.write_char('\n')
  builder.write_string("source_rows,")
  builder.write_string(manifest.source_rows.to_string())
  builder.write_char('\n')
  builder.write_string("event_count,")
  builder.write_string(manifest.event_count.to_string())
  builder.write_char('\n')
  builder.write_string("generated_by,")
  builder.write_string(manifest.generated_by)
  builder.write_char('\n')
  builder.write_string("checksum,")
  builder.write_string(manifest.checksum)
  builder.write_char('\n')
  builder.to_string()
}
