///|
/// A delimiter configuration for text exports.
pub(all) struct SurvivalExportFormat {
  name : String
  delimiter : Char
  line_ending : String
  include_header : Bool
  format_row : Bool
} derive(Eq, Debug)

///|
/// A standard CSV format.
pub fn survival_csv_format() -> SurvivalExportFormat {
  {
    name: "csv",
    delimiter: ',',
    line_ending: "\n",
    include_header: true,
    format_row: true,
  }
}

///|
/// A tab-separated format.
pub fn survival_tsv_format() -> SurvivalExportFormat {
  {
    name: "tsv",
    delimiter: '\t',
    line_ending: "\n",
    include_header: true,
    format_row: true,
  }
}

///|
/// Escape a field for a delimiter format.
pub fn survival_export_escape(
  value : String,
  format : SurvivalExportFormat,
) -> String {
  if value.contains(format.delimiter.to_string()) ||
    value.contains("\n") ||
    value.contains("\"") {
    "\"" + value.replace(old="\"", new="\"\"") + "\""
  } else {
    value
  }
}

///|
/// Serialize a numeric matrix with a selected format.
pub fn survival_matrix_export(
  matrix : Array[Array[Double]],
  headers : Array[String],
  format : SurvivalExportFormat,
) -> String {
  let builder = StringBuilder::new()
  if format.include_header {
    for i in 0..<headers.length() {
      if i > 0 {
        builder.write_char(format.delimiter)
      }
      builder.write_string(survival_export_escape(headers[i], format))
    }
    builder.write_string(format.line_ending)
  }
  for row in matrix {
    for i in 0..<row.length() {
      if i > 0 {
        builder.write_char(format.delimiter)
      }
      builder.write_string(row[i].to_string())
    }
    builder.write_string(format.line_ending)
  }
  builder.to_string()
}

///|
/// Serialize one row as a JSON-like object.
pub fn survival_row_json(
  names : Array[String],
  values : Array[Double],
) -> String {
  let builder = StringBuilder::new()
  builder.write_char('{')
  let count = if names.length() < values.length() {
    names.length()
  } else {
    values.length()
  }
  for i in 0..<count {
    if i > 0 {
      builder.write_char(',')
    }
    builder.write_char('"')
    builder.write_string(names[i])
    builder.write_string("\":")
    builder.write_string(values[i].to_string())
  }
  builder.write_char('}')
  builder.to_string()
}

///|
/// Serialize a numeric matrix as line-delimited JSON-like objects.
pub fn survival_matrix_json(
  names : Array[String],
  matrix : Array[Array[Double]],
) -> String {
  let builder = StringBuilder::new()
  for row in matrix {
    builder.write_string(survival_row_json(names, row))
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Build an HTML table for a report preview.
pub fn survival_matrix_html(
  headers : Array[String],
  matrix : Array[Array[Double]],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("<table><thead><tr>")
  for header in headers {
    builder.write_string("<th>")
    builder.write_string(header)
    builder.write_string("</th>")
  }
  builder.write_string("</tr></thead><tbody>")
  for row in matrix {
    builder.write_string("<tr>")
    for value in row {
      builder.write_string("<td>")
      builder.write_string(value.to_string())
      builder.write_string("</td>")
    }
    builder.write_string("</tr>")
  }
  builder.write_string("</tbody></table>")
  builder.to_string()
}

///|
/// Convert a format to a stable label.
pub fn survival_export_format_name(format : SurvivalExportFormat) -> String {
  format.name
}

///|
/// Check an export format.
pub fn survival_export_format_is_valid(format : SurvivalExportFormat) -> Bool {
  format.name != "" && format.line_ending != ""
}
