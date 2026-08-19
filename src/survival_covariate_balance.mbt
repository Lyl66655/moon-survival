///|
/// A covariate balance result.
pub(all) struct SurvivalBalanceRow {
  covariate : Int
  mean_first : Double
  mean_second : Double
  standardized_difference : Double
  variance_ratio : Double
  balanced : Bool
  balance_row : Bool
} derive(Eq, Debug)

///|
/// Calculate balance for a numeric design matrix.
pub fn survival_covariate_balance(
  values : Array[Array[Double]],
  groups : Array[Int],
  threshold : Double,
) -> Array[SurvivalBalanceRow] {
  let result : Array[SurvivalBalanceRow] = []
  if values.length() == 0 {
    return result
  }
  let width = values[0].length()
  for column in 0..<width {
    let first : Array[Double] = []
    let second : Array[Double] = []
    for i in 0..<values.length() {
      if i < groups.length() && column < values[i].length() {
        if groups[i] == 0 {
          first.push(values[i][column])
        } else {
          second.push(values[i][column])
        }
      }
    }
    let first_variance = survival_ext_variance(first)
    let second_variance = survival_ext_variance(second)
    let standard_difference = survival_ext_safe_ratio(
      survival_ext_mean(first) - survival_ext_mean(second),
      double_sqrt((first_variance + second_variance) / 2.0),
    )
    let ratio = survival_ext_safe_ratio(first_variance, second_variance)
    result.push({
      covariate: column,
      mean_first: survival_ext_mean(first),
      mean_second: survival_ext_mean(second),
      standardized_difference: standard_difference,
      variance_ratio: ratio,
      balanced: standard_difference.abs() <= threshold,
      balance_row: true,
    })
  }
  result
}

///|
/// Compute the number of balanced covariates.
pub fn survival_balanced_count(rows : Array[SurvivalBalanceRow]) -> Int {
  let mut count = 0
  for row in rows {
    if row.balanced {
      count += 1
    }
  }
  count
}

///|
/// Compute the maximum absolute standardized difference.
pub fn survival_max_balance_difference(
  rows : Array[SurvivalBalanceRow],
) -> Double {
  let mut maximum = 0.0
  for row in rows {
    if row.standardized_difference.abs() > maximum {
      maximum = row.standardized_difference.abs()
    }
  }
  maximum
}

///|
/// Compare balance tables before and after weighting.
pub fn survival_balance_improvement(
  before : Array[SurvivalBalanceRow],
  after : Array[SurvivalBalanceRow],
) -> Double {
  survival_max_balance_difference(before) -
  survival_max_balance_difference(after)
}

///|
/// Export balance rows.
pub fn survival_balance_to_csv(rows : Array[SurvivalBalanceRow]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "covariate,mean_first,mean_second,standardized_difference,variance_ratio,balanced\n",
  )
  for row in rows {
    builder.write_string(row.covariate.to_string())
    builder.write_char(',')
    builder.write_string(row.mean_first.to_string())
    builder.write_char(',')
    builder.write_string(row.mean_second.to_string())
    builder.write_char(',')
    builder.write_string(row.standardized_difference.to_string())
    builder.write_char(',')
    builder.write_string(row.variance_ratio.to_string())
    builder.write_char(',')
    builder.write_string(row.balanced.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate balance rows.
pub fn survival_balance_is_valid(rows : Array[SurvivalBalanceRow]) -> Bool {
  for row in rows {
    if row.covariate < 0 ||
      row.standardized_difference.is_nan() ||
      row.variance_ratio.is_nan() {
      return false
    }
  }
  true
}
