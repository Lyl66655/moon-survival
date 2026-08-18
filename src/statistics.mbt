///|
/// Descriptive statistics for numeric vectors.
pub(all) struct NumericSummary {
  count : Int
  mean : Double
  variance : Double
  standard_deviation : Double
  minimum : Double
  first_quartile : Double
  median : Double
  third_quartile : Double
  maximum : Double
} derive(Eq, Debug)

///|
/// Return a sorted copy of a numeric vector.
pub fn sorted_doubles(values : Array[Double]) -> Array[Double] {
  let result = []
  for value in values {
    result.push(value)
  }
  result.sort()
  result
}

///|
/// Linear-interpolated empirical quantile.
pub fn quantile_doubles(
  values : Array[Double],
  probability : Double,
) -> Double? {
  if values.length() == 0 {
    return None
  }
  let sorted = sorted_doubles(values)
  let p = clamp_probability(probability)
  if sorted.length() == 1 {
    return Some(sorted[0])
  }
  let position = p * (sorted.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 >= sorted.length() { lower } else { lower + 1 }
  let fraction = position - lower.to_double()
  Some(sorted[lower] + fraction * (sorted[upper] - sorted[lower]))
}

///|
/// Sample variance with Bessel correction.
pub fn sample_variance(values : Array[Double]) -> Double {
  if values.length() < 2 {
    0.0
  } else {
    let mean = mean_doubles(values)
    let mut sum = 0.0
    for value in values {
      let difference = value - mean
      sum += difference * difference
    }
    sum / (values.length() - 1).to_double()
  }
}

///|
/// Population variance.
pub fn population_variance(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mean = mean_doubles(values)
    let mut sum = 0.0
    for value in values {
      let difference = value - mean
      sum += difference * difference
    }
    sum / values.length().to_double()
  }
}

///|
/// Pearson correlation of two vectors.
pub fn correlation(a : Array[Double], b : Array[Double]) -> Double {
  if a.length() != b.length() || a.length() < 2 {
    panic()
  }
  let mean_a = mean_doubles(a)
  let mean_b = mean_doubles(b)
  let mut numerator = 0.0
  let mut sum_a = 0.0
  let mut sum_b = 0.0
  for i in 0..<a.length() {
    let da = a[i] - mean_a
    let db = b[i] - mean_b
    numerator += da * db
    sum_a += da * da
    sum_b += db * db
  }
  if sum_a <= survival_epsilon || sum_b <= survival_epsilon {
    0.0
  } else {
    numerator / double_sqrt(sum_a * sum_b)
  }
}

///|
/// Covariance of two vectors.
pub fn covariance(a : Array[Double], b : Array[Double]) -> Double {
  if a.length() != b.length() || a.length() < 2 {
    panic()
  }
  let mean_a = mean_doubles(a)
  let mean_b = mean_doubles(b)
  let mut total = 0.0
  for i in 0..<a.length() {
    total += (a[i] - mean_a) * (b[i] - mean_b)
  }
  total / (a.length() - 1).to_double()
}

///|
/// Compute a complete numeric summary.
pub fn numeric_summary(values : Array[Double]) -> NumericSummary {
  if values.length() == 0 {
    return {
      count: 0,
      mean: 0.0,
      variance: 0.0,
      standard_deviation: 0.0,
      minimum: 0.0,
      first_quartile: 0.0,
      median: 0.0,
      third_quartile: 0.0,
      maximum: 0.0,
    }
  }
  let sorted = sorted_doubles(values)
  {
    count: values.length(),
    mean: mean_doubles(values),
    variance: sample_variance(values),
    standard_deviation: double_sqrt(sample_variance(values)),
    minimum: sorted[0],
    first_quartile: quantile_doubles(sorted, 0.25).unwrap_or(0.0),
    median: quantile_doubles(sorted, 0.5).unwrap_or(0.0),
    third_quartile: quantile_doubles(sorted, 0.75).unwrap_or(0.0),
    maximum: sorted[sorted.length() - 1],
  }
}

///|
/// Column means of a rectangular matrix.
pub fn column_means(matrix : Array[Array[Double]]) -> Array[Double] {
  if matrix.length() == 0 {
    return []
  }
  let columns = matrix[0].length()
  let means = zeros(columns)
  for row in matrix {
    if row.length() != columns {
      panic()
    }
    for j in 0..<columns {
      means[j] += row[j]
    }
  }
  for j in 0..<columns {
    means[j] /= matrix.length().to_double()
  }
  means
}

///|
/// Column standard deviations of a rectangular matrix.
pub fn column_standard_deviations(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let means = column_means(matrix)
  let deviations = zeros(means.length())
  if matrix.length() < 2 {
    return deviations
  }
  for row in matrix {
    for j in 0..<means.length() {
      let difference = row[j] - means[j]
      deviations[j] += difference * difference
    }
  }
  for j in 0..<means.length() {
    deviations[j] = double_sqrt(
      deviations[j] / (matrix.length() - 1).to_double(),
    )
  }
  deviations
}

///|
/// Center and scale columns, returning means and scales for later prediction.
pub(all) struct Standardization {
  means : Array[Double]
  scales : Array[Double]
  values : Array[Array[Double]]
} derive(Eq, Debug)

///|
/// Standardize each covariate column using sample standard deviations.
pub fn standardize_matrix(matrix : Array[Array[Double]]) -> Standardization {
  let means = column_means(matrix)
  let scales = column_standard_deviations(matrix)
  let values = []
  for row in matrix {
    let output = []
    for j in 0..<means.length() {
      let scale = if scales[j] <= survival_epsilon { 1.0 } else { scales[j] }
      output.push((row[j] - means[j]) / scale)
    }
    values.push(output)
  }
  { means, scales, values }
}

///|
/// Apply an existing standardization to one row.
pub fn Standardization::transform(
  self : Standardization,
  row : Array[Double],
) -> Array[Double] {
  if row.length() != self.means.length() {
    panic()
  }
  let result = []
  for j in 0..<row.length() {
    let scale = if self.scales[j] <= survival_epsilon {
      1.0
    } else {
      self.scales[j]
    }
    result.push((row[j] - self.means[j]) / scale)
  }
  result
}

///|
/// Compute the Brier score for event occurrence by a horizon.
pub fn brier_score(
  data : SurvivalData,
  predictions : Array[Double],
  horizon : Double,
) -> Double {
  if data.length() != predictions.length() {
    panic()
  }
  let mut total = 0.0
  for i in 0..<data.length() {
    let observation = data.get(i)
    let outcome = if observation.time <= horizon && observation.event == 1 {
      1.0
    } else {
      0.0
    }
    let difference = clamp_probability(predictions[i]) - outcome
    total += difference * difference
  }
  if data.length() == 0 {
    0.0
  } else {
    total / data.length().to_double()
  }
}

///|
/// Integrated Brier score across a vector of horizons.
pub fn integrated_brier_score(
  data : SurvivalData,
  predictions : Array[Array[Double]],
  horizons : Array[Double],
) -> Double {
  if predictions.length() != horizons.length() {
    panic()
  }
  if horizons.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for i in 0..<horizons.length() {
    total += brier_score(data, predictions[i], horizons[i])
  }
  total / horizons.length().to_double()
}

///|
/// Schoenfeld-style residual summary for a single Cox covariate.
pub(all) struct ResidualSummary {
  mean : Double
  variance : Double
  standard_deviation : Double
  minimum : Double
  maximum : Double
} derive(Eq, Debug)

///|
/// Summarize a residual vector.
pub fn residual_summary(values : Array[Double]) -> ResidualSummary {
  if values.length() == 0 {
    {
      mean: 0.0,
      variance: 0.0,
      standard_deviation: 0.0,
      minimum: 0.0,
      maximum: 0.0,
    }
  } else {
    let sorted = sorted_doubles(values)
    let variance = sample_variance(values)
    {
      mean: mean_doubles(values),
      variance,
      standard_deviation: double_sqrt(variance),
      minimum: sorted[0],
      maximum: sorted[sorted.length() - 1],
    }
  }
}

///|
/// Estimate an empirical event rate by fixed-width time bins.
pub(all) struct EventRateBin {
  start : Double
  end : Double
  n_at_risk : Int
  events : Int
  event_rate : Double
} derive(Eq, Debug)

///|
/// Build a descriptive event-rate histogram for reliability dashboards.
pub fn event_rate_bins(
  data : SurvivalData,
  start : Double,
  end : Double,
  bins : Int,
) -> Array[EventRateBin] {
  if bins <= 0 || end <= start {
    return []
  }
  let result = []
  let width = (end - start) / bins.to_double()
  for bin in 0..<bins {
    let bin_start = start + bin.to_double() * width
    let bin_end = if bin == bins - 1 { end } else { bin_start + width }
    let mut risk = 0
    let mut events = 0
    for observation in data.to_array() {
      if observation.time >= bin_start {
        risk += 1
      }
      if observation.event == 1 &&
        observation.time >= bin_start &&
        observation.time < bin_end {
        events += 1
      }
    }
    result.push({
      start: bin_start,
      end: bin_end,
      n_at_risk: risk,
      events,
      event_rate: if risk == 0 {
        0.0
      } else {
        events.to_double() / risk.to_double()
      },
    })
  }
  result
}
