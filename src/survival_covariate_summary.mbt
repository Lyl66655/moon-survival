///|
/// A numeric covariate summary.
pub(all) struct SurvivalCovariateSummary {
  column : Int
  count : Int
  missing : Int
  mean : Double
  standard_deviation : Double
  minimum : Double
  lower_quartile : Double
  median : Double
  upper_quartile : Double
  maximum : Double
} derive(Eq, Debug)

///|
/// Summarize each numeric covariate column.
pub fn survival_covariate_summaries(
  values : Array[Array[Double]],
) -> Array[SurvivalCovariateSummary] {
  let result : Array[SurvivalCovariateSummary] = []
  let width = if values.length() == 0 { 0 } else { values[0].length() }
  for column in 0..<width {
    let observed : Array[Double] = []
    let mut missing = 0
    for row in values {
      if column >= row.length() || row[column].is_nan() || row[column].is_inf() {
        missing += 1
      } else {
        observed.push(row[column])
      }
    }
    let sorted = survival_ext_sorted(observed)
    result.push({
      column,
      count: observed.length(),
      missing,
      mean: survival_ext_mean(observed),
      standard_deviation: survival_ext_standard_deviation(observed),
      minimum: survival_ext_min(observed),
      lower_quartile: survival_ext_quantile(sorted, 0.25),
      median: survival_ext_quantile(sorted, 0.5),
      upper_quartile: survival_ext_quantile(sorted, 0.75),
      maximum: survival_ext_max(observed),
    })
  }
  result
}

///|
/// Compute a column correlation matrix.
pub fn survival_covariate_correlation(
  values : Array[Array[Double]],
) -> Array[Array[Double]] {
  let summaries = survival_covariate_summaries(values)
  let width = summaries.length()
  let result : Array[Array[Double]] = []
  for first in 0..<width {
    let row : Array[Double] = []
    for second in 0..<width {
      let first_values : Array[Double] = []
      let second_values : Array[Double] = []
      for source in values {
        if first < source.length() &&
          second < source.length() &&
          !source[first].is_nan() &&
          !source[second].is_nan() {
          first_values.push(source[first])
          second_values.push(source[second])
        }
      }
      let first_mean = survival_ext_mean(first_values)
      let second_mean = survival_ext_mean(second_values)
      let mut covariance = 0.0
      let mut first_variance = 0.0
      let mut second_variance = 0.0
      for i in 0..<first_values.length() {
        let left = first_values[i] - first_mean
        let right = second_values[i] - second_mean
        covariance += left * right
        first_variance += left * left
        second_variance += right * right
      }
      row.push(
        survival_ext_safe_ratio(
          covariance,
          double_sqrt(first_variance * second_variance),
        ),
      )
    }
    result.push(row)
  }
  result
}

///|
/// Return unique observed levels for a categorical column.
pub fn survival_covariate_levels(
  values : Array[Array[Double]],
  column : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  for row in values {
    if column < row.length() &&
      !row[column].is_nan() &&
      !result.contains(row[column]) {
      result.push(row[column])
    }
  }
  result.sort()
  result
}

///|
/// Count missing cells by row.
pub fn survival_covariate_missing_by_row(
  values : Array[Array[Double]],
) -> Array[Int] {
  let result : Array[Int] = []
  for row in values {
    let mut missing = 0
    for value in row {
      if value.is_nan() || value.is_inf() {
        missing += 1
      }
    }
    result.push(missing)
  }
  result
}

///|
/// Standardize one covariate column.
pub fn survival_covariate_standardize(
  values : Array[Array[Double]],
  column : Int,
) -> Array[Double] {
  let source : Array[Double] = []
  for row in values {
    if column < row.length() {
      source.push(row[column])
    }
  }
  let mean = survival_ext_mean(source)
  let scale = survival_ext_standard_deviation(source)
  let result : Array[Double] = []
  for value in source {
    result.push(survival_ext_safe_ratio(value - mean, scale))
  }
  result
}

///|
/// Export covariate summaries.
pub fn survival_covariate_summaries_to_csv(
  summaries : Array[SurvivalCovariateSummary],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "column,count,missing,mean,standard_deviation,minimum,lower_quartile,median,upper_quartile,maximum\n",
  )
  for summary in summaries {
    builder.write_string(summary.column.to_string())
    builder.write_char(',')
    builder.write_string(summary.count.to_string())
    builder.write_char(',')
    builder.write_string(summary.missing.to_string())
    builder.write_char(',')
    builder.write_string(summary.mean.to_string())
    builder.write_char(',')
    builder.write_string(summary.standard_deviation.to_string())
    builder.write_char(',')
    builder.write_string(summary.minimum.to_string())
    builder.write_char(',')
    builder.write_string(summary.lower_quartile.to_string())
    builder.write_char(',')
    builder.write_string(summary.median.to_string())
    builder.write_char(',')
    builder.write_string(summary.upper_quartile.to_string())
    builder.write_char(',')
    builder.write_string(summary.maximum.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate covariate summaries.
pub fn survival_covariate_summary_is_valid(
  summary : SurvivalCovariateSummary,
) -> Bool {
  summary.count >= 0 &&
  summary.missing >= 0 &&
  summary.minimum <= summary.lower_quartile &&
  summary.lower_quartile <= summary.median &&
  summary.median <= summary.upper_quartile &&
  summary.upper_quartile <= summary.maximum
}
