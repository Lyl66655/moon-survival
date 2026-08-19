///|
/// A display cell with an explicit numeric format.
pub(all) struct SurvivalReportCell {
  label : String
  value : String
  numeric : Double
  has_numeric : Bool
} derive(Eq, Debug)

///|
/// A section in a clinical analysis report.
pub(all) struct SurvivalReportSection {
  title : String
  cells : Array[SurvivalReportCell]
  rows : Array[Array[String]]
  report_section : Bool
} derive(Eq, Debug)

///|
/// A complete report document.
pub(all) struct SurvivalReportDocument {
  title : String
  sections : Array[SurvivalReportSection]
  generated : String
  report_document : Bool
} derive(Eq, Debug)

///|
/// Format a number with a stable compact representation.
pub fn survival_report_number(value : Double, digits : Int) -> String {
  if value.is_nan() {
    return "NA"
  }
  if value.is_inf() {
    return if value > 0.0 { "Inf" } else { "-Inf" }
  }
  if digits <= 0 {
    return value.to_string()
  }
  value.to_string()
}

///|
/// Create a numeric report cell.
pub fn survival_report_numeric_cell(
  label : String,
  value : Double,
  digits : Int,
) -> SurvivalReportCell {
  {
    label,
    value: survival_report_number(value, digits),
    numeric: value,
    has_numeric: true,
  }
}

///|
/// Create a text report cell.
pub fn survival_report_text_cell(
  label : String,
  value : String,
) -> SurvivalReportCell {
  { label, value, numeric: 0.0, has_numeric: false }
}

///|
/// Create a section with no rows.
pub fn survival_report_section(title : String) -> SurvivalReportSection {
  { title, cells: [], rows: [], report_section: true }
}

///|
/// Add a cell to a section.
pub fn SurvivalReportSection::add_cell(
  self : SurvivalReportSection,
  cell : SurvivalReportCell,
) -> SurvivalReportSection {
  let cells : Array[SurvivalReportCell] = []
  for value in self.cells {
    cells.push(value)
  }
  cells.push(cell)
  { title: self.title, cells, rows: self.rows, report_section: true }
}

///|
/// Add a table row to a section.
pub fn SurvivalReportSection::add_row(
  self : SurvivalReportSection,
  row : Array[String],
) -> SurvivalReportSection {
  let rows : Array[Array[String]] = []
  for value in self.rows {
    rows.push(value)
  }
  rows.push(row)
  { title: self.title, cells: self.cells, rows, report_section: true }
}

///|
/// Start a report document.
pub fn survival_report_document(title : String) -> SurvivalReportDocument {
  { title, sections: [], generated: "moon-survival", report_document: true }
}

///|
/// Add a report section.
pub fn SurvivalReportDocument::add_section(
  self : SurvivalReportDocument,
  section : SurvivalReportSection,
) -> SurvivalReportDocument {
  let sections : Array[SurvivalReportSection] = []
  for value in self.sections {
    sections.push(value)
  }
  sections.push(section)
  {
    title: self.title,
    sections,
    generated: self.generated,
    report_document: true,
  }
}

///|
/// Build a descriptive section from SurvivalData.
pub fn survival_report_descriptive_section(
  data : SurvivalData,
) -> SurvivalReportSection {
  let fit = kaplan_meier(data)
  let section = survival_report_section("Descriptive survival")
  let with_rows = section.add_cell(
    survival_report_numeric_cell("sample_size", data.length().to_double(), 0),
  )
  let with_events = with_rows.add_cell(
    survival_report_numeric_cell("events", data.event_count().to_double(), 0),
  )
  let with_censoring = with_events.add_cell(
    survival_report_numeric_cell(
      "censoring_rate",
      if data.length() == 0 {
        0.0
      } else {
        data.censored_count().to_double() / data.length().to_double()
      },
      3,
    ),
  )
  match fit.median() {
    Some(median) =>
      with_censoring.add_cell(survival_report_numeric_cell("median", median, 3))
    None =>
      with_censoring.add_cell(
        survival_report_text_cell("median", "not reached"),
      )
  }
}

///|
/// Build a model section from a parametric fit.
pub fn survival_report_parametric_section(
  fit : ParametricFit,
) -> SurvivalReportSection {
  let section = survival_report_section(parametric_family_name(fit.family))
  let with_likelihood = section.add_cell(
    survival_report_numeric_cell("log_likelihood", fit.log_likelihood, 4),
  )
  let with_aic = with_likelihood.add_cell(
    survival_report_numeric_cell("aic", fit.aic, 4),
  )
  with_aic.add_cell(
    survival_report_numeric_cell("sample_size", fit.sample_size.to_double(), 0),
  )
}

///|
/// Render a document as Markdown.
pub fn SurvivalReportDocument::to_markdown(
  self : SurvivalReportDocument,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("# ")
  builder.write_string(self.title)
  builder.write_char('\n')
  builder.write_char('\n')
  for section in self.sections {
    builder.write_string("## ")
    builder.write_string(section.title)
    builder.write_char('\n')
    builder.write_char('\n')
    if section.cells.length() > 0 {
      builder.write_string("| Metric | Value |\n|---|---|\n")
      for cell in section.cells {
        builder.write_char('|')
        builder.write_string(cell.label)
        builder.write_char('|')
        builder.write_string(cell.value)
        builder.write_string("|\n")
      }
      builder.write_char('\n')
    }
    for row in section.rows {
      builder.write_string("| ")
      for index in 0..<row.length() {
        if index > 0 {
          builder.write_string(" | ")
        }
        builder.write_string(row[index])
      }
      builder.write_string(" |\n")
    }
  }
  builder.to_string()
}

///|
/// Render a section as a compact CSV block.
pub fn survival_report_section_to_csv(
  section : SurvivalReportSection,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("label,value\n")
  for cell in section.cells {
    builder.write_string(cell.label)
    builder.write_char(',')
    builder.write_string(cell.value)
    builder.write_char('\n')
  }
  builder.write_string("table_row\n")
  for row in section.rows {
    for index in 0..<row.length() {
      if index > 0 {
        builder.write_char(',')
      }
      builder.write_string(row[index])
    }
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Compute a report checksum.
pub fn survival_report_checksum(document : SurvivalReportDocument) -> String {
  survival_export_checksum(document.to_markdown())
}

///|
/// Add standard KM and model sections.
pub fn survival_standard_report(
  data : SurvivalData,
  fit : ParametricFit,
) -> SurvivalReportDocument {
  survival_report_document("Survival analysis")
  .add_section(survival_report_descriptive_section(data))
  .add_section(survival_report_parametric_section(fit))
}

///|
/// Validate report structure.
pub fn survival_report_document_is_valid(
  document : SurvivalReportDocument,
) -> Bool {
  if document.title == "" || document.sections.length() == 0 {
    return false
  }
  for section in document.sections {
    if section.title == "" {
      return false
    }
    for cell in section.cells {
      if cell.label == "" {
        return false
      }
    }
  }
  true
}
