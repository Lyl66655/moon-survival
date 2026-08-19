///|
/// Summary of missingness for one covariate column.
pub(all) struct MissingColumnSummary {
  column : Int
  observed : Int
  missing : Int
  missing_fraction : Double
  mean_observed : Double
  median_observed : Double
  minimum : Double
  maximum : Double
} derive(Eq, Debug)

///|
/// A missingness pattern and the rows that share it.
pub(all) struct MissingPattern {
  pattern : String
  row_count : Int
  row_indices : Array[Int]
  pattern_row : Bool
} derive(Eq, Debug)

///|
/// Report produced before an imputation decision.
pub(all) struct MissingnessReport {
  row_count : Int
  column_count : Int
  complete_rows : Int
  incomplete_rows : Int
  columns : Array[MissingColumnSummary]
  patterns : Array[MissingPattern]
} derive(Eq, Debug)

///|
/// An imputed design matrix and audit metadata.
pub(all) struct ImputationResult {
  values : Array[Array[Double]]
  indicators : Array[Array[Int]]
  method_name : String
  imputed_cells : Int
  row_count : Int
  column_count : Int
} derive(Eq, Debug)

///|
/// Rubin-style pooled estimate across imputed analyses.
pub(all) struct RubinPooledEstimate {
  estimate : Double
  within_variance : Double
  between_variance : Double
  total_variance : Double
  standard_error : Double
  relative_increase : Double
  imputations : Int
} derive(Eq, Debug)

///|
/// Determine whether a value is missing according to the package contract.
pub fn missing_value(value : Double) -> Bool {
  value.is_nan() || value.is_inf()
}

///|
/// Return a binary missingness indicator matrix.
pub fn missingness_indicators(
  values : Array[Array[Double]],
) -> Array[Array[Int]] {
  let result : Array[Array[Int]] = []
  for row in values {
    let output : Array[Int] = []
    for value in row {
      output.push(if missing_value(value) { 1 } else { 0 })
    }
    result.push(output)
  }
  result
}

///|
/// Return a compact pattern key such as `0100` for a row.
pub fn missing_pattern_key(row : Array[Double]) -> String {
  let builder = StringBuilder::new()
  for value in row {
    builder.write_string(if missing_value(value) { "1" } else { "0" })
  }
  builder.to_string()
}

///|
/// Compute a robust summary for one column.
pub fn missing_column_summary(
  values : Array[Array[Double]],
  column : Int,
) -> MissingColumnSummary {
  let observed_values : Array[Double] = []
  let mut missing = 0
  for row in values {
    if column < 0 || column >= row.length() || missing_value(row[column]) {
      missing += 1
    } else {
      observed_values.push(row[column])
    }
  }
  observed_values.sort()
  let count = observed_values.length()
  let mean = if count == 0 { 0.0 } else { mean_doubles(observed_values) }
  let median = if count == 0 { 0.0 } else { observed_values[count / 2] }
  let minimum = if count == 0 { 0.0 } else { observed_values[0] }
  let maximum = if count == 0 { 0.0 } else { observed_values[count - 1] }
  {
    column,
    observed: count,
    missing,
    missing_fraction: if values.length() == 0 {
      0.0
    } else {
      missing.to_double() / values.length().to_double()
    },
    mean_observed: mean,
    median_observed: median,
    minimum,
    maximum,
  }
}

///|
/// Generate unique missingness patterns and row membership.
pub fn missing_patterns(values : Array[Array[Double]]) -> Array[MissingPattern] {
  let result : Array[MissingPattern] = []
  for index in 0..<values.length() {
    let key = missing_pattern_key(values[index])
    let mut found = false
    for pattern in result {
      if pattern.pattern == key {
        found = true
      }
    }
    if !found {
      result.push({
        pattern: key,
        row_count: 1,
        row_indices: [index],
        pattern_row: true,
      })
    } else {
      for i in 0..<result.length() {
        if result[i].pattern == key {
          let current = result[i]
          let indices = current.row_indices
          indices.push(index)
          result[i] = {
            pattern: current.pattern,
            row_count: current.row_count + 1,
            row_indices: indices,
            pattern_row: true,
          }
        }
      }
    }
  }
  result
}

///|
/// Build a complete missingness report.
pub fn missingness_report(values : Array[Array[Double]]) -> MissingnessReport {
  let column_count = if values.length() == 0 { 0 } else { values[0].length() }
  let columns : Array[MissingColumnSummary] = []
  for column in 0..<column_count {
    columns.push(missing_column_summary(values, column))
  }
  let patterns = missing_patterns(values)
  let mut incomplete = 0
  for pattern in patterns {
    if pattern.pattern.contains("1") {
      incomplete += pattern.row_count
    }
  }
  {
    row_count: values.length(),
    column_count,
    complete_rows: values.length() - incomplete,
    incomplete_rows: incomplete,
    columns,
    patterns,
  }
}

///|
/// Copy a matrix while replacing missing values with a column statistic.
fn missing_impute_with(
  values : Array[Array[Double]],
  use_median : Bool,
) -> ImputationResult {
  let report = missingness_report(values)
  let replacement : Array[Double] = []
  for summary in report.columns {
    replacement.push(
      if use_median {
        summary.median_observed
      } else {
        summary.mean_observed
      },
    )
  }
  let output : Array[Array[Double]] = []
  let indicators : Array[Array[Int]] = []
  let mut imputed = 0
  for row in values {
    let copied : Array[Double] = []
    let row_indicators : Array[Int] = []
    for column in 0..<report.column_count {
      let value = row[column]
      if missing_value(value) {
        copied.push(replacement[column])
        row_indicators.push(1)
        imputed += 1
      } else {
        copied.push(value)
        row_indicators.push(0)
      }
    }
    output.push(copied)
    indicators.push(row_indicators)
  }
  {
    values: output,
    indicators,
    method_name: if use_median {
      "median"
    } else {
      "mean"
    },
    imputed_cells: imputed,
    row_count: report.row_count,
    column_count: report.column_count,
  }
}

///|
/// Mean imputation with an explicit audit matrix.
pub fn missing_mean_impute(values : Array[Array[Double]]) -> ImputationResult {
  missing_impute_with(values, false)
}

///|
/// Median imputation for skewed laboratory or duration covariates.
pub fn missing_median_impute(values : Array[Array[Double]]) -> ImputationResult {
  missing_impute_with(values, true)
}

///|
/// Last-observation-carried-forward imputation for longitudinal covariates.
pub fn missing_locf_impute(values : Array[Array[Double]]) -> ImputationResult {
  let report = missingness_report(values)
  let output : Array[Array[Double]] = []
  let indicators : Array[Array[Int]] = []
  let last : Array[Double] = Array::make(report.column_count, 0.0)
  let observed : Array[Int] = Array::make(report.column_count, 0)
  let mut imputed = 0
  for row in values {
    let copied : Array[Double] = []
    let row_indicators : Array[Int] = []
    for column in 0..<report.column_count {
      if missing_value(row[column]) {
        copied.push(
          if observed[column] == 1 {
            last[column]
          } else {
            report.columns[column].mean_observed
          },
        )
        row_indicators.push(1)
        imputed += 1
      } else {
        copied.push(row[column])
        last[column] = row[column]
        observed[column] = 1
        row_indicators.push(0)
      }
    }
    output.push(copied)
    indicators.push(row_indicators)
  }
  {
    values: output,
    indicators,
    method_name: "locf",
    imputed_cells: imputed,
    row_count: report.row_count,
    column_count: report.column_count,
  }
}

///|
/// Add missingness indicators to an imputed design matrix.
pub fn missing_add_indicators(
  result : ImputationResult,
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row_index in 0..<result.values.length() {
    let row : Array[Double] = []
    for value in result.values[row_index] {
      row.push(value)
    }
    for indicator in result.indicators[row_index] {
      row.push(indicator.to_double())
    }
    output.push(row)
  }
  output
}

///|
/// Standardize an imputed matrix using observed-scale statistics.
pub fn missing_standardize(result : ImputationResult) -> Array[Array[Double]] {
  let means : Array[Double] = []
  let scales : Array[Double] = []
  for column in 0..<result.column_count {
    let column_values : Array[Double] = []
    for row in result.values {
      column_values.push(row[column])
    }
    let mean = if column_values.length() == 0 {
      0.0
    } else {
      mean_doubles(column_values)
    }
    let mut variance = 0.0
    for value in column_values {
      let delta = value - mean
      variance += delta * delta
    }
    means.push(mean)
    scales.push(
      double_sqrt(variance / (column_values.length() + 1).to_double()),
    )
  }
  let output : Array[Array[Double]] = []
  for row in result.values {
    let scaled : Array[Double] = []
    for column in 0..<result.column_count {
      scaled.push(
        (row[column] - means[column]) / clamp_positive(scales[column]),
      )
    }
    output.push(scaled)
  }
  output
}

///|
/// Deterministically create multiple mean-imputed datasets with small
/// posterior-like perturbations for sensitivity analysis.
pub fn missing_multiple_impute(
  values : Array[Array[Double]],
  count : Int,
) -> Array[ImputationResult] {
  let result : Array[ImputationResult] = []
  if count <= 0 {
    return result
  }
  let base = missing_mean_impute(values)
  for replicate in 0..<count {
    let output : Array[Array[Double]] = []
    for row_index in 0..<base.values.length() {
      let row : Array[Double] = []
      for column in 0..<base.column_count {
        let value = base.values[row_index][column]
        let direction = if (row_index + column + replicate) % 2 == 0 {
          1.0
        } else {
          -1.0
        }
        row.push(
          if base.indicators[row_index][column] == 1 {
            value + direction * 0.000001 * (replicate + 1).to_double()
          } else {
            value
          },
        )
      }
      output.push(row)
    }
    result.push({
      values: output,
      indicators: base.indicators,
      method_name: "mean_multiple_" + (replicate + 1).to_string(),
      imputed_cells: base.imputed_cells,
      row_count: base.row_count,
      column_count: base.column_count,
    })
  }
  result
}

///|
/// Pool point estimates and variances using Rubin's rules.
pub fn missing_rubin_pool(
  estimates : Array[Double],
  variances : Array[Double],
) -> RubinPooledEstimate {
  if estimates.length() == 0 || estimates.length() != variances.length() {
    return {
      estimate: 0.0,
      within_variance: 0.0,
      between_variance: 0.0,
      total_variance: 0.0,
      standard_error: 0.0,
      relative_increase: 0.0,
      imputations: 0,
    }
  }
  let estimate = mean_doubles(estimates)
  let within = mean_doubles(variances)
  let mut between = 0.0
  for value in estimates {
    let delta = value - estimate
    between += delta * delta
  }
  between /= (estimates.length() - 1).max(1).to_double()
  let total = within + (1.0 + 1.0 / estimates.length().to_double()) * between
  {
    estimate,
    within_variance: within,
    between_variance: between,
    total_variance: total,
    standard_error: double_sqrt(total),
    relative_increase: if within <= survival_epsilon {
      0.0
    } else {
      total / within - 1.0
    },
    imputations: estimates.length(),
  }
}

///|
/// Return complete-case row indices.
pub fn missing_complete_case_indices(
  values : Array[Array[Double]],
) -> Array[Int] {
  let result : Array[Int] = []
  for index in 0..<values.length() {
    let mut complete = true
    for value in values[index] {
      if missing_value(value) {
        complete = false
      }
    }
    if complete {
      result.push(index)
    }
  }
  result
}

///|
/// Export a missingness report.
pub fn MissingnessReport::to_csv(self : MissingnessReport) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "column,observed,missing,missing_fraction,mean_observed,median_observed,minimum,maximum\n",
  )
  for column in self.columns {
    builder.write_string(column.column.to_string())
    builder.write_char(',')
    builder.write_string(column.observed.to_string())
    builder.write_char(',')
    builder.write_string(column.missing.to_string())
    builder.write_char(',')
    builder.write_string(column.missing_fraction.to_string())
    builder.write_char(',')
    builder.write_string(column.mean_observed.to_string())
    builder.write_char(',')
    builder.write_string(column.median_observed.to_string())
    builder.write_char(',')
    builder.write_string(column.minimum.to_string())
    builder.write_char(',')
    builder.write_string(column.maximum.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate an imputation result.
pub fn imputation_is_valid(result : ImputationResult) -> Bool {
  if result.values.length() != result.indicators.length() ||
    result.row_count != result.values.length() {
    return false
  }
  for row_index in 0..<result.values.length() {
    if result.values[row_index].length() != result.column_count ||
      result.indicators[row_index].length() != result.column_count {
      return false
    }
    for column in 0..<result.column_count {
      if result.values[row_index][column].is_nan() ||
        result.values[row_index][column].is_inf() ||
        (
          result.indicators[row_index][column] != 0 &&
          result.indicators[row_index][column] != 1
        ) {
        return false
      }
    }
  }
  true
}
