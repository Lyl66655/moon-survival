///|
/// CSV parsing failures with row and column context retained by the caller.
pub(all) enum CsvError {
  EmptyInput
  UnterminatedQuote
  RowWidthMismatch
  MissingColumn(String)
  InvalidNumber(String)
  InvalidEventValue
  DataFailure(DataError)
} derive(Eq, Debug)

///|
/// A small RFC 4180-compatible table representation.
pub(all) struct CsvTable {
  headers : Array[String]
  rows : Array[Array[String]]
} derive(Eq, Debug)

///|
/// Parse CSV with quoted fields, escaped quotes, CRLF, and configurable
/// delimiter. Empty trailing fields are preserved.
pub fn csv_parse(
  text : String,
  delimiter? : Char = ',',
) -> Result[CsvTable, CsvError] {
  if text.is_empty() {
    return Err(EmptyInput)
  }
  let characters = text.to_array()
  let rows : Array[Array[String]] = []
  let row : Array[String] = []
  let field = StringBuilder::new()
  let mut in_quotes = false
  let mut index = 0
  while index < characters.length() {
    let character = characters[index]
    if in_quotes {
      if character == '"' {
        if index + 1 < characters.length() && characters[index + 1] == '"' {
          field.write_char('"')
          index += 1
        } else {
          in_quotes = false
        }
      } else {
        field.write_char(character)
      }
    } else if character == '"' && field.is_empty() {
      in_quotes = true
    } else if character == delimiter {
      row.push(field.to_string())
      field.reset()
    } else if character == '\n' || character == '\r' {
      if character == '\r' &&
        index + 1 < characters.length() &&
        characters[index + 1] == '\n' {
        index += 1
      }
      row.push(field.to_string())
      field.reset()
      if row.length() > 1 || !row[0].trim().is_empty() {
        let copied = []
        for value in row {
          copied.push(value)
        }
        rows.push(copied)
      }
      row.clear()
    } else {
      field.write_char(character)
    }
    index += 1
  }
  if in_quotes {
    return Err(UnterminatedQuote)
  }
  row.push(field.to_string())
  if row.length() > 1 || !row[0].trim().is_empty() {
    let copied = []
    for value in row {
      copied.push(value)
    }
    rows.push(copied)
  }
  if rows.length() == 0 {
    return Err(EmptyInput)
  }
  let headers = []
  for header in rows[0] {
    headers.push(header.trim().to_owned())
  }
  let data_rows : Array[Array[String]] = []
  for i in 1..<rows.length() {
    if rows[i].length() != headers.length() {
      return Err(RowWidthMismatch)
    }
    let values = []
    for value in rows[i] {
      values.push(value.trim().to_owned())
    }
    data_rows.push(values)
  }
  Ok({ headers, rows: data_rows })
}

///|
/// Find a column by name, returning its zero-based table index.
pub fn CsvTable::column(self : CsvTable, name : String) -> Int? {
  for i in 0..<self.headers.length() {
    if self.headers[i] == name {
      return Some(i)
    }
  }
  None
}

///|
/// Parse a numeric column into doubles.
pub fn CsvTable::double_column(
  self : CsvTable,
  name : String,
) -> Result[Array[Double], CsvError] {
  let column = match self.column(name) {
    Some(value) => value
    None => return Err(MissingColumn(name))
  }
  let result = []
  for row in self.rows {
    match parse_csv_double(row[column]) {
      Ok(value) => result.push(value)
      Err(error) => return Err(error)
    }
  }
  Ok(result)
}

///|
/// Parse an event indicator column.
pub fn CsvTable::event_column(
  self : CsvTable,
  name : String,
) -> Result[Array[Int], CsvError] {
  let column = match self.column(name) {
    Some(value) => value
    None => return Err(MissingColumn(name))
  }
  let result = []
  for row in self.rows {
    match parse_csv_int(row[column]) {
      Ok(value) => {
        if value != 0 && value != 1 {
          return Err(InvalidEventValue)
        }
        result.push(value)
      }
      Err(error) => return Err(error)
    }
  }
  Ok(result)
}

///|
/// Load ordinary survival data from named CSV columns.
pub fn survival_data_from_csv(
  text : String,
  time_column : String,
  event_column : String,
) -> Result[SurvivalData, CsvError] {
  let table = match csv_parse(text) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let times = match table.double_column(time_column) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let events = match table.event_column(event_column) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  match SurvivalData::try_new(times, events) {
    Ok(value) => Ok(value)
    Err(error) => Err(DataFailure(error))
  }
}

///|
/// Load a Cox design matrix from named columns.
pub fn cox_data_from_csv(
  text : String,
  time_column : String,
  event_column : String,
  covariate_columns : Array[String],
) -> Result[CoxData, CsvError] {
  let table = match csv_parse(text) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let times = match table.double_column(time_column) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let events = match table.event_column(event_column) {
    Ok(value) => value
    Err(error) => return Err(error)
  }
  let columns = []
  for name in covariate_columns {
    match table.double_column(name) {
      Ok(value) => columns.push(value)
      Err(error) => return Err(error)
    }
  }
  let covariates = []
  for i in 0..<times.length() {
    let row = []
    for column in columns {
      row.push(column[i])
    }
    covariates.push(row)
  }
  match CoxData::try_new(times, events, covariates) {
    Ok(value) => Ok(value)
    Err(error) => Err(DataFailure(error))
  }
}

///|
/// Escape one CSV cell.
pub fn csv_escape(value : String) -> String {
  if value.contains_char(',') ||
    value.contains_char('"') ||
    value.contains_char('\n') ||
    value.contains_char('\r') {
    let builder = StringBuilder::new()
    builder.write_char('"')
    for character in value.to_array() {
      if character == '"' {
        builder.write_char('"')
        builder.write_char('"')
      } else {
        builder.write_char(character)
      }
    }
    builder.write_char('"')
    builder.to_string()
  } else {
    value
  }
}

///|
/// Serialize ordinary survival observations to a stable CSV fixture.
pub fn survival_data_to_csv(data : SurvivalData) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,event\n")
  for observation in data.to_array() {
    builder.write_string("\{observation.time},\{observation.event}\n")
  }
  builder.to_string()
}

///|
/// Serialize a table with a header and quoted cells.
pub fn csv_stringify(table : CsvTable, delimiter? : Char = ',') -> String {
  let builder = StringBuilder::new()
  write_csv_row(builder, table.headers, delimiter)
  for row in table.rows {
    write_csv_row(builder, row, delimiter)
  }
  builder.to_string()
}

///|
/// Parse a floating-point CSV field.
fn parse_csv_double(value : String) -> Result[Double, CsvError] {
  let parsed = try @string.parse_double(value.to_string_view()) catch {
    _ => None
  } noraise {
    value => Some(value)
  }
  match parsed {
    Some(number) => Ok(number)
    None => Err(InvalidNumber(value))
  }
}

///|
/// Parse an integer CSV field.
fn parse_csv_int(value : String) -> Result[Int, CsvError] {
  let parsed = try @string.parse_int(value.to_string_view()) catch {
    _ => None
  } noraise {
    value => Some(value)
  }
  match parsed {
    Some(number) => Ok(number)
    None => Err(InvalidNumber(value))
  }
}

///|
/// Write a row with configurable delimiter and escaping.
fn write_csv_row(
  builder : StringBuilder,
  row : Array[String],
  delimiter : Char,
) -> Unit {
  for i in 0..<row.length() {
    if i > 0 {
      builder.write_char(delimiter)
    }
    builder.write_string(csv_escape(row[i]))
  }
  builder.write_string("\n")
}
