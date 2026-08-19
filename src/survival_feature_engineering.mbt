///|
/// Deterministic feature engineering helpers for survival model inputs.
pub(all) struct SurvivalFeatureColumn {
  name : String
  values : Array[Double]
  missing_count : Int
  mean : Double
  standard_deviation : Double
  feature_row : Int
} derive(Eq, Debug)

///|
pub(all) struct SurvivalFeatureBin {
  bin_row : Int
  lower : Double
  upper : Double
  count : Int
  mean : Double
  event_rate : Double
  feature_row : Int
} derive(Eq, Debug)

///|
/// Winsorize a numeric feature at lower and upper probabilities.
pub fn survival_feature_winsorize(
  values : Array[Double],
  lower_probability : Double,
  upper_probability : Double,
) -> Array[Double] {
  let lower = survival_ext_quantile(values, lower_probability)
  let upper = survival_ext_quantile(
    values,
    upper_probability.max(lower_probability),
  )
  let output : Array[Double] = []
  for value in values {
    if value.is_nan() {
      output.push(value)
    } else {
      output.push(value.max(lower).min(upper))
    }
  }
  output
}

///|
/// Standardize a feature while preserving missing values as NaN.
pub fn survival_feature_zscore(values : Array[Double]) -> Array[Double] {
  let mean = survival_ext_mean(values)
  let standard_deviation = survival_ext_standard_deviation(values)
  let output : Array[Double] = []
  for value in values {
    if value.is_nan() {
      output.push(value)
    } else if standard_deviation <= survival_epsilon {
      output.push(0.0)
    } else {
      output.push((value - mean) / standard_deviation)
    }
  }
  output
}

///|
/// Create a feature column with stable summary statistics.
pub fn survival_feature_column(
  name : String,
  values : Array[Double],
  feature_row : Int,
) -> SurvivalFeatureColumn {
  let mut missing = 0
  for value in values {
    if value.is_nan() {
      missing += 1
    }
  }
  {
    name,
    values,
    missing_count: missing,
    mean: survival_ext_mean(values),
    standard_deviation: survival_ext_standard_deviation(values),
    feature_row,
  }
}

///|
/// Assign values to equal-width bins and summarize events in each bin.
pub fn survival_feature_bins(
  values : Array[Double],
  events : Array[Int],
  bin_count : Int,
) -> Array[SurvivalFeatureBin] {
  let output : Array[SurvivalFeatureBin] = []
  if bin_count <= 0 || values.length() == 0 {
    return output
  }
  let minimum = survival_ext_min(values)
  let maximum = survival_ext_max(values)
  let width = if maximum <= minimum {
    1.0
  } else {
    (maximum - minimum) / bin_count.to_double()
  }
  let counts = Array::make(bin_count, 0)
  let event_counts = Array::make(bin_count, 0)
  let sums = Array::make(bin_count, 0.0)
  let count = values.length().min(events.length())
  for i = 0; i < count; i = i + 1 {
    let raw = ((values[i] - minimum) / width).floor().to_int()
    let bin = raw.max(0).min(bin_count - 1)
    counts[bin] += 1
    sums[bin] += values[i]
    if events[i] == 1 {
      event_counts[bin] += 1
    }
  }
  for bin = 0; bin < bin_count; bin = bin + 1 {
    let lower = minimum + bin.to_double() * width
    let upper = if bin + 1 == bin_count { maximum } else { lower + width }
    let mean = if counts[bin] == 0 {
      0.0
    } else {
      sums[bin] / counts[bin].to_double()
    }
    let rate = if counts[bin] == 0 {
      0.0
    } else {
      event_counts[bin].to_double() / counts[bin].to_double()
    }
    output.push({
      bin_row: bin,
      lower,
      upper,
      count: counts[bin],
      mean,
      event_rate: rate,
      feature_row: bin,
    })
  }
  output
}

///|
/// Add an intercept column to a matrix of feature rows.
pub fn survival_feature_add_intercept(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    let expanded : Array[Double] = [1.0]
    for value in row {
      expanded.push(value)
    }
    output.push(expanded)
  }
  output
}

///|
/// Apply a log1p transform to nonnegative time or exposure values.
pub fn survival_feature_log1p(values : Array[Double]) -> Array[Double] {
  let output : Array[Double] = []
  for value in values {
    output.push(if value <= -1.0 { 0.0 } else { double_log(1.0 + value) })
  }
  output
}

///|
/// Return one-hot indicator columns for integer categories.
pub fn survival_feature_one_hot(
  values : Array[Int],
  levels : Array[Int],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for value in values {
    let row : Array[Double] = []
    for level in levels {
      row.push(if value == level { 1.0 } else { 0.0 })
    }
    output.push(row)
  }
  output
}

///|
pub fn survival_feature_bins_to_csv(bins : Array[SurvivalFeatureBin]) -> String {
  let mut out = "bin_row,lower,upper,count,mean,event_rate,feature_row\n"
  for row in bins {
    out = out +
      row.bin_row.to_string() +
      "," +
      row.lower.to_string() +
      "," +
      row.upper.to_string() +
      "," +
      row.count.to_string() +
      "," +
      row.mean.to_string() +
      "," +
      row.event_rate.to_string() +
      "," +
      row.feature_row.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_feature_bins_are_valid(
  bins : Array[SurvivalFeatureBin],
) -> Bool {
  for row in bins {
    if row.lower > row.upper ||
      row.count < 0 ||
      row.mean.is_nan() ||
      row.event_rate < 0.0 ||
      row.event_rate > 1.0 {
      return false
    }
  }
  true
}
