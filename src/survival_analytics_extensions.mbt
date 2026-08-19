///|
/// Return the smallest finite value in a vector.
pub fn survival_ext_min(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut result = values[0]
  for value in values {
    if !value.is_nan() && value < result {
      result = value
    }
  }
  result
}

///|
/// Return the largest finite value in a vector.
pub fn survival_ext_max(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut result = values[0]
  for value in values {
    if !value.is_nan() && value > result {
      result = value
    }
  }
  result
}

///|
/// Sum finite values while skipping missing entries.
pub fn survival_ext_sum(values : Array[Double]) -> Double {
  let mut total = 0.0
  for value in values {
    if !value.is_nan() && !value.is_inf() {
      total += value
    }
  }
  total
}

///|
/// Mean finite values.
pub fn survival_ext_mean(values : Array[Double]) -> Double {
  let mut total = 0.0
  let mut count = 0
  for value in values {
    if !value.is_nan() && !value.is_inf() {
      total += value
      count += 1
    }
  }
  if count == 0 {
    0.0
  } else {
    total / count.to_double()
  }
}

///|
/// Sample variance with a stable empty convention.
pub fn survival_ext_variance(values : Array[Double]) -> Double {
  let mean = survival_ext_mean(values)
  let mut total = 0.0
  let mut count = 0
  for value in values {
    if !value.is_nan() && !value.is_inf() {
      let delta = value - mean
      total += delta * delta
      count += 1
    }
  }
  if count < 2 {
    0.0
  } else {
    total / (count - 1).to_double()
  }
}

///|
/// Standard deviation of finite values.
pub fn survival_ext_standard_deviation(values : Array[Double]) -> Double {
  double_sqrt(survival_ext_variance(values))
}

///|
/// Return a sorted copy of a vector.
pub fn survival_ext_sorted(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  for value in values {
    result.push(value)
  }
  result.sort()
  result
}

///|
/// Return a linear-interpolated quantile.
pub fn survival_ext_quantile(
  values : Array[Double],
  probability : Double,
) -> Double {
  let sorted = survival_ext_sorted(values)
  if sorted.length() == 0 {
    return 0.0
  }
  let position = clamp_probability(probability) *
    (sorted.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 >= sorted.length() { lower } else { lower + 1 }
  sorted[lower] +
  (sorted[upper] - sorted[lower]) * (position - lower.to_double())
}

///|
/// Normalize a vector to a unit sum.
pub fn survival_ext_normalize(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let total = survival_ext_sum(values)
  for value in values {
    result.push(
      if total.abs() <= survival_epsilon {
        0.0
      } else {
        value / total
      },
    )
  }
  result
}

///|
/// Cumulative sum.
pub fn survival_ext_cumulative_sum(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let mut total = 0.0
  for value in values {
    total += value
    result.push(total)
  }
  result
}

///|
/// Cumulative product.
pub fn survival_ext_cumulative_product(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let mut product = 1.0
  for value in values {
    product *= value
    result.push(product)
  }
  result
}

///|
/// First differences.
pub fn survival_ext_differences(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  if values.length() < 2 {
    return result
  }
  for index in 1..<values.length() {
    result.push(values[index] - values[index - 1])
  }
  result
}

///|
/// Trapezoid integration on a grid.
pub fn survival_ext_trapezoid_area(
  times : Array[Double],
  values : Array[Double],
) -> Double {
  let count = if times.length() < values.length() {
    times.length()
  } else {
    values.length()
  }
  if count < 2 {
    return 0.0
  }
  let mut area = 0.0
  for index in 1..<count {
    area += (times[index] - times[index - 1]) *
      (values[index] + values[index - 1]) /
      2.0
  }
  area
}

///|
/// Right-continuous step integration.
pub fn survival_ext_step_area(
  times : Array[Double],
  values : Array[Double],
  horizon : Double,
) -> Double {
  if horizon <= 0.0 || times.length() == 0 || values.length() == 0 {
    return 0.0
  }
  let mut area = 0.0
  let mut previous = 0.0
  let mut value = values[0]
  for index in 0..<times.length() {
    if index >= values.length() || times[index] >= horizon {
      break
    }
    area += (times[index] - previous) * value
    previous = times[index]
    value = values[index]
  }
  area + (horizon - previous) * value
}

///|
/// Linear interpolation on an ascending grid.
pub fn survival_ext_interpolate(
  times : Array[Double],
  values : Array[Double],
  query : Double,
) -> Double {
  let count = if times.length() < values.length() {
    times.length()
  } else {
    values.length()
  }
  if count == 0 {
    return 0.0
  }
  if query <= times[0] {
    return values[0]
  }
  for index in 1..<count {
    if query <= times[index] {
      let width = times[index] - times[index - 1]
      if width.abs() <= survival_epsilon {
        return values[index]
      }
      let fraction = (query - times[index - 1]) / width
      return values[index - 1] + fraction * (values[index] - values[index - 1])
    }
  }
  values[count - 1]
}

///|
/// Interpolate a complete query grid.
pub fn survival_ext_interpolate_many(
  times : Array[Double],
  values : Array[Double],
  queries : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  for query in queries {
    result.push(survival_ext_interpolate(times, values, query))
  }
  result
}

///|
/// Clamp a vector to an interval.
pub fn survival_ext_clamp_vector(
  values : Array[Double],
  lower : Double,
  upper : Double,
) -> Array[Double] {
  let result : Array[Double] = []
  for value in values {
    result.push(
      if value < lower {
        lower
      } else if value > upper {
        upper
      } else {
        value
      },
    )
  }
  result
}

///|
/// Apply a log transform safely.
pub fn survival_ext_log_vector(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  for value in values {
    result.push(double_log(value))
  }
  result
}

///|
/// Apply a stable exponential transform.
pub fn survival_ext_exp_vector(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  for value in values {
    result.push(double_exp(value))
  }
  result
}

///|
/// Compute a logit transform.
pub fn survival_ext_logit(probability : Double) -> Double {
  let p = clamp_probability(probability)
  double_log(p / clamp_positive(1.0 - p))
}

///|
/// Compute an inverse logit.
pub fn survival_ext_inv_logit(value : Double) -> Double {
  clamp_probability(1.0 / (1.0 + double_exp(-value)))
}

///|
/// Compute entropy of a probability vector.
pub fn survival_ext_entropy(probabilities : Array[Double]) -> Double {
  let mut entropy = 0.0
  for probability in probabilities {
    let p = clamp_positive(probability)
    entropy -= p * double_log(p)
  }
  entropy
}

///|
/// Compute the Gini impurity of a probability vector.
pub fn survival_ext_gini(probabilities : Array[Double]) -> Double {
  let mut squared = 0.0
  for probability in probabilities {
    squared += probability * probability
  }
  1.0 - squared
}

///|
/// Weighted mean.
pub fn survival_ext_weighted_mean(
  values : Array[Double],
  weights : Array[Double],
) -> Double {
  if values.length() != weights.length() {
    return 0.0
  }
  let mut total = 0.0
  let mut denominator = 0.0
  for index in 0..<values.length() {
    total += values[index] * weights[index]
    denominator += weights[index]
  }
  if denominator.abs() <= survival_epsilon {
    0.0
  } else {
    total / denominator
  }
}

///|
/// Weighted variance.
pub fn survival_ext_weighted_variance(
  values : Array[Double],
  weights : Array[Double],
) -> Double {
  let mean = survival_ext_weighted_mean(values, weights)
  let mut total = 0.0
  let mut denominator = 0.0
  let count = if values.length() < weights.length() {
    values.length()
  } else {
    weights.length()
  }
  for index in 0..<count {
    let delta = values[index] - mean
    total += weights[index] * delta * delta
    denominator += weights[index]
  }
  if denominator <= survival_epsilon {
    0.0
  } else {
    total / denominator
  }
}

///|
/// Compute ranks with average handling for ties.
pub fn survival_ext_ranks(values : Array[Double]) -> Array[Double] {
  let result = Array::make(values.length(), 0.0)
  for index in 0..<values.length() {
    let less = values.filter(fn(value) { value < values[index] }).length()
    let equal = values.filter(fn(value) { value == values[index] }).length()
    result[index] = less.to_double() + (equal + 1).to_double() / 2.0
  }
  result
}

///|
/// Generate integer order indices.
pub fn survival_ext_order(values : Array[Double]) -> Array[Int] {
  let result : Array[Int] = []
  for index in 0..<values.length() {
    result.push(index)
  }
  result.sort_by(fn(a, b) {
    if values[a] < values[b] {
      -1
    } else if values[a] > values[b] {
      1
    } else {
      a.compare(b)
    }
  })
  result
}

///|
/// Select values by integer indices.
pub fn survival_ext_select(
  values : Array[Double],
  indices : Array[Int],
) -> Array[Double] {
  let result : Array[Double] = []
  for index in indices {
    if index >= 0 && index < values.length() {
      result.push(values[index])
    }
  }
  result
}

///|
/// Return unique sorted values.
pub fn survival_ext_unique(values : Array[Double]) -> Array[Double] {
  let sorted = survival_ext_sorted(values)
  let result : Array[Double] = []
  for value in sorted {
    if result.length() == 0 || result[result.length() - 1] != value {
      result.push(value)
    }
  }
  result
}

///|
/// Generate an inclusive regular grid.
pub fn survival_ext_grid(
  start : Double,
  end : Double,
  points : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  if points <= 0 {
    return result
  }
  if points == 1 {
    result.push(start)
    return result
  }
  for index in 0..<points {
    result.push(
      start + (end - start) * index.to_double() / (points - 1).to_double(),
    )
  }
  result
}

///|
/// Return the nearest grid value.
pub fn survival_ext_nearest(grid : Array[Double], query : Double) -> Double {
  if grid.length() == 0 {
    return 0.0
  }
  let mut best = grid[0]
  let mut distance = (grid[0] - query).abs()
  for value in grid {
    if (value - query).abs() < distance {
      best = value
      distance = (value - query).abs()
    }
  }
  best
}

///|
/// Return the index of the nearest grid value.
pub fn survival_ext_nearest_index(grid : Array[Double], query : Double) -> Int {
  if grid.length() == 0 {
    return -1
  }
  let mut best = 0
  let mut distance = (grid[0] - query).abs()
  for index in 1..<grid.length() {
    if (grid[index] - query).abs() < distance {
      best = index
      distance = (grid[index] - query).abs()
    }
  }
  best
}

///|
/// Safe ratio with an explicit zero denominator convention.
pub fn survival_ext_safe_ratio(
  numerator : Double,
  denominator : Double,
) -> Double {
  if denominator.abs() <= survival_epsilon {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Compute a vector dot product with length checking.
pub fn survival_ext_dot(left : Array[Double], right : Array[Double]) -> Double {
  if left.length() != right.length() {
    return 0.0
  }
  dot_product(left, right)
}

///|
/// Add vectors with length checking.
pub fn survival_ext_add(
  left : Array[Double],
  right : Array[Double],
) -> Array[Double] {
  if left.length() != right.length() {
    return []
  }
  let result : Array[Double] = []
  for index in 0..<left.length() {
    result.push(left[index] + right[index])
  }
  result
}

///|
/// Scale a vector.
pub fn survival_ext_scale(
  values : Array[Double],
  factor : Double,
) -> Array[Double] {
  let result : Array[Double] = []
  for value in values {
    result.push(value * factor)
  }
  result
}

///|
/// Euclidean norm.
pub fn survival_ext_norm(values : Array[Double]) -> Double {
  double_sqrt(survival_ext_dot(values, values))
}

///|
/// Count finite observations.
pub fn survival_ext_finite_count(values : Array[Double]) -> Int {
  let mut count = 0
  for value in values {
    if !value.is_nan() && !value.is_inf() {
      count += 1
    }
  }
  count
}

///|
/// Return a boolean finite mask as integers.
pub fn survival_ext_finite_mask(values : Array[Double]) -> Array[Int] {
  let result : Array[Int] = []
  for value in values {
    result.push(if value.is_nan() || value.is_inf() { 0 } else { 1 })
  }
  result
}

///|
/// Return a vector of all zeros.
pub fn survival_ext_zeros(length : Int) -> Array[Double] {
  Array::make(if length < 0 { 0 } else { length }, 0.0)
}

///|
/// Return a vector of all ones.
pub fn survival_ext_ones(length : Int) -> Array[Double] {
  Array::make(if length < 0 { 0 } else { length }, 1.0)
}

///|
/// Compute the maximum absolute difference between two vectors.
pub fn survival_ext_max_distance(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  if left.length() != right.length() {
    return 1.0e300
  }
  let mut result = 0.0
  for index in 0..<left.length() {
    if (left[index] - right[index]).abs() > result {
      result = (left[index] - right[index]).abs()
    }
  }
  result
}

///|
/// Return whether a sequence is non-decreasing.
pub fn survival_ext_is_monotone(values : Array[Double]) -> Bool {
  for index in 1..<values.length() {
    if values[index] < values[index - 1] {
      return false
    }
  }
  true
}

///|
/// Return a monotone cumulative maximum.
pub fn survival_ext_cumulative_max(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let mut current = -1.0e300
  for value in values {
    if value > current {
      current = value
    }
    result.push(current)
  }
  result
}

///|
/// Return a monotone cumulative minimum.
pub fn survival_ext_cumulative_min(values : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let mut current = 1.0e300
  for value in values {
    if value < current {
      current = value
    }
    result.push(current)
  }
  result
}

///|
/// Compute a rolling mean with a fixed window.
fn survival_ext_window(
  values : Array[Double],
  start : Int,
  end : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  for index in start..<end {
    if index >= 0 && index < values.length() {
      result.push(values[index])
    }
  }
  result
}

///|
/// Compute a rolling mean with a fixed window.
pub fn survival_ext_rolling_mean(
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  if window <= 0 {
    return result
  }
  for index in 0..<values.length() {
    let start = if index + 1 < window { 0 } else { index + 1 - window }
    result.push(
      survival_ext_mean(survival_ext_window(values, start, index + 1)),
    )
  }
  result
}

///|
/// Compute a rolling maximum.
pub fn survival_ext_rolling_max(
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  if window <= 0 {
    return result
  }
  for index in 0..<values.length() {
    let start = if index + 1 < window { 0 } else { index + 1 - window }
    result.push(survival_ext_max(survival_ext_window(values, start, index + 1)))
  }
  result
}

///|
/// Compute a rolling minimum.
pub fn survival_ext_rolling_min(
  values : Array[Double],
  window : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  if window <= 0 {
    return result
  }
  for index in 0..<values.length() {
    let start = if index + 1 < window { 0 } else { index + 1 - window }
    result.push(survival_ext_min(survival_ext_window(values, start, index + 1)))
  }
  result
}

///|
/// Compute a finite difference rate.
pub fn survival_ext_rate(
  times : Array[Double],
  values : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  let count = if times.length() < values.length() {
    times.length()
  } else {
    values.length()
  }
  if count < 2 {
    return result
  }
  for index in 1..<count {
    result.push(
      survival_ext_safe_ratio(
        values[index] - values[index - 1],
        times[index] - times[index - 1],
      ),
    )
  }
  result
}

///|
/// Return an odds ratio from two probabilities.
pub fn survival_ext_odds_ratio(first : Double, second : Double) -> Double {
  survival_ext_safe_ratio(first * (1.0 - second), second * (1.0 - first))
}

///|
/// Convert a hazard ratio to an absolute risk using a baseline risk.
pub fn survival_ext_absolute_risk(
  baseline_risk : Double,
  hazard_ratio : Double,
) -> Double {
  let baseline = clamp_probability(baseline_risk)
  let hazard = clamp_positive(hazard_ratio)
  clamp_probability(1.0 - double_exp(hazard * double_log(1.0 - baseline)))
}

///|
/// Return a number of equal-width time bins.
pub fn survival_ext_time_bin(time : Double, width : Double, bins : Int) -> Int {
  if width <= 0.0 || bins <= 0 || time < 0.0 {
    return -1
  }
  let index = (time / width).to_int()
  if index >= bins {
    bins - 1
  } else {
    index
  }
}

///|
/// Return a binary event outcome through a horizon.
pub fn survival_ext_horizon_outcome(
  data : SurvivalData,
  horizon : Double,
) -> Array[Int] {
  let result : Array[Int] = []
  for observation in data.to_array() {
    result.push(
      if observation.event == 1 && observation.time <= horizon {
        1
      } else {
        0
      },
    )
  }
  result
}

///|
/// Return a binary at-risk indicator at a horizon.
pub fn survival_ext_at_risk_outcome(
  data : SurvivalData,
  horizon : Double,
) -> Array[Int] {
  let result : Array[Int] = []
  for observation in data.to_array() {
    result.push(if observation.time >= horizon { 1 } else { 0 })
  }
  result
}

///|
/// Serialize a numeric vector.
pub fn survival_ext_values_to_csv(values : Array[Double]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("index,value\n")
  for index in 0..<values.length() {
    builder.write_string(index.to_string())
    builder.write_char(',')
    builder.write_string(values[index].to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
