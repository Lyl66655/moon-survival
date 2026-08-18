///|
/// Numerical constants used by the estimators.
let survival_epsilon : Double = 0.000000000001

///|
/// A numerically stable exponential with a practical range guard.
pub fn double_exp(x : Double) -> Double {
  if x > 700.0 {
    @math.exp(700.0)
  } else if x < -745.0 {
    0.0
  } else {
    @math.exp(x)
  }
}

///|
/// Natural logarithm that never evaluates the logarithm of zero in a
/// likelihood calculation.
pub fn double_log(x : Double) -> Double {
  if x <= survival_epsilon {
    @math.ln(survival_epsilon)
  } else {
    @math.ln(x)
  }
}

///|
/// Square root helper retained as part of the original public API.
pub fn double_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    x.sqrt()
  }
}

///|
/// Absolute difference for tolerance based comparisons.
pub fn double_distance(a : Double, b : Double) -> Double {
  (a - b).abs()
}

///|
/// Return whether two floating point values are close under a mixed
/// relative/absolute tolerance.
pub fn double_is_close(
  a : Double,
  b : Double,
  relative_tolerance? : Double = 0.0000001,
  absolute_tolerance? : Double = 0.0000000001,
) -> Bool {
  let scale = if a.abs() > b.abs() { a.abs() } else { b.abs() }
  (a - b).abs() <= absolute_tolerance ||
  (a - b).abs() <= relative_tolerance * scale
}

///|
/// Clamp a probability to its mathematically valid interval.
pub fn clamp_probability(x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else if x > 1.0 {
    1.0
  } else {
    x
  }
}

///|
/// Clamp a positive quantity away from zero before it is used as a divisor.
pub fn clamp_positive(x : Double, floor? : Double = survival_epsilon) -> Double {
  if x < floor {
    floor
  } else {
    x
  }
}

///|
/// Sum a vector without allocating an intermediate collection.
pub fn sum_doubles(values : Array[Double]) -> Double {
  let mut total = 0.0
  for value in values {
    total += value
  }
  total
}

///|
/// Mean of a vector. Empty vectors have the neutral value zero.
pub fn mean_doubles(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    sum_doubles(values) / values.length().to_double()
  }
}

///|
/// Weighted sum for inverse probability and reliability calculations.
pub fn weighted_sum(values : Array[Double], weights : Array[Double]) -> Double {
  if values.length() != weights.length() {
    panic()
  }
  let mut total = 0.0
  for i in 0..<values.length() {
    total += values[i] * weights[i]
  }
  total
}

///|
/// Weighted mean with a safe zero-weight result.
pub fn weighted_mean(values : Array[Double], weights : Array[Double]) -> Double {
  let denominator = sum_doubles(weights)
  if denominator <= survival_epsilon {
    0.0
  } else {
    weighted_sum(values, weights) / denominator
  }
}

///|
/// Dot product of two vectors.
pub fn dot_product(a : Array[Double], b : Array[Double]) -> Double {
  if a.length() != b.length() {
    panic()
  }
  let mut total = 0.0
  for i in 0..<a.length() {
    total += a[i] * b[i]
  }
  total
}

///|
/// Euclidean norm of a vector.
pub fn vector_norm(values : Array[Double]) -> Double {
  double_sqrt(dot_product(values, values))
}

///|
/// Elementwise vector addition.
pub fn vector_add(a : Array[Double], b : Array[Double]) -> Array[Double] {
  if a.length() != b.length() {
    panic()
  }
  let result = []
  for i in 0..<a.length() {
    result.push(a[i] + b[i])
  }
  result
}

///|
/// Elementwise vector subtraction.
pub fn vector_subtract(a : Array[Double], b : Array[Double]) -> Array[Double] {
  if a.length() != b.length() {
    panic()
  }
  let result = []
  for i in 0..<a.length() {
    result.push(a[i] - b[i])
  }
  result
}

///|
/// Scale every component of a vector.
pub fn vector_scale(values : Array[Double], scale : Double) -> Array[Double] {
  let result = []
  for value in values {
    result.push(value * scale)
  }
  result
}

///|
/// Return a zero vector of the requested dimension.
pub fn zeros(length : Int) -> Array[Double] {
  let result = []
  for _ in 0..<length {
    result.push(0.0)
  }
  result
}

///|
/// Return an identity matrix.
pub fn identity_matrix(size : Int) -> Array[Array[Double]] {
  let matrix = []
  for i in 0..<size {
    let row = []
    for j in 0..<size {
      row.push(if i == j { 1.0 } else { 0.0 })
    }
    matrix.push(row)
  }
  matrix
}

///|
/// Copy a vector while applying a scalar transform.
pub fn map_doubles(
  values : Array[Double],
  f : (Double) -> Double,
) -> Array[Double] {
  let result = []
  for value in values {
    result.push(f(value))
  }
  result
}

///|
/// Stable log-sum-exp for a vector of log weights.
pub fn log_sum_exp(values : Array[Double]) -> Double {
  if values.length() == 0 {
    double_log(0.0)
  } else {
    let mut maximum = values[0]
    for i in 1..<values.length() {
      if values[i] > maximum {
        maximum = values[i]
      }
    }
    let mut total = 0.0
    for value in values {
      total += double_exp(value - maximum)
    }
    maximum + double_log(total)
  }
}

///|
/// A deterministic pseudo-random generator used by examples and benchmarks.
pub(all) struct DeterministicRng {
  mut state : UInt
} derive(Debug, Eq)

///|
/// Create a deterministic generator. A zero seed is mapped to a non-zero
/// state so that the recurrence never becomes stuck.
pub fn DeterministicRng::new(seed : UInt) -> DeterministicRng {
  { state: if seed == 0U { 2654435761U } else { seed } }
}

///|
/// Advance the generator and return a uniformly distributed 32-bit value.
pub fn DeterministicRng::next_u32(self : DeterministicRng) -> UInt {
  self.state = self.state * 1664525U + 1013904223U
  self.state
}

///|
/// Return a deterministic uniform value in [0, 1).
pub fn DeterministicRng::uniform(self : DeterministicRng) -> Double {
  self.next_u32().to_double() / 4294967296.0
}

///|
/// Return a deterministic standard normal value using Box-Muller.
pub fn DeterministicRng::normal(self : DeterministicRng) -> Double {
  let u1 = clamp_positive(self.uniform(), floor=0.0000000001)
  let u2 = self.uniform()
  double_sqrt(-2.0 * double_log(u1)) * @math.cos(2.0 * 3.141592653589793 * u2)
}

///|
/// Convert an integer array to a floating point array.
pub fn ints_to_doubles(values : Array[Int]) -> Array[Double] {
  let result = []
  for value in values {
    result.push(value.to_double())
  }
  result
}

///|
/// Produce a monotonically increasing sequence of evenly spaced values.
pub fn linspace(start : Double, end : Double, count : Int) -> Array[Double] {
  if count <= 0 {
    []
  } else if count == 1 {
    [start]
  } else {
    let result = []
    let step = (end - start) / (count - 1).to_double()
    for i in 0..<count {
      result.push(start + i.to_double() * step)
    }
    result
  }
}

///|
/// Numerically stable logistic function.
pub fn logistic(x : Double) -> Double {
  if x >= 0.0 {
    let z = double_exp(-x)
    1.0 / (1.0 + z)
  } else {
    let z = double_exp(x)
    z / (1.0 + z)
  }
}

///|
/// Logit transform with clipping at the probability boundaries.
pub fn logit(p : Double) -> Double {
  let q = clamp_probability(p)
  if q <= 0.0 {
    -700.0
  } else if q >= 1.0 {
    700.0
  } else {
    double_log(q / (1.0 - q))
  }
}

///|
/// Standard normal density.
pub fn normal_density(x : Double) -> Double {
  double_exp(-0.5 * x * x) / double_sqrt(2.0 * 3.141592653589793)
}

///|
/// Approximation to the standard normal CDF with maximum error suitable for
/// confidence intervals in ordinary survival-analysis reports.
pub fn normal_cdf(x : Double) -> Double {
  let sign = if x < 0.0 { -1.0 } else { 1.0 }
  let z = x.abs() / double_sqrt(2.0)
  let t = 1.0 / (1.0 + 0.3275911 * z)
  let polynomial = 1.0 -
    (
      (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) *
      t +
      0.254829592
    ) *
    t *
    double_exp(-z * z)
  0.5 * (1.0 + sign * polynomial)
}

///|
/// Approximate a standard normal quantile using a rational approximation.
pub fn normal_quantile(p : Double) -> Double {
  let q = if p < 0.0000001 {
    0.0000001
  } else if p > 0.9999999 {
    0.9999999
  } else {
    p
  }
  let a = [
    -39.6968302866538, 220.946098424521, -275.928510446969, 138.357751867269, -30.6647980661472,
    2.50662827745924,
  ]
  let b = [
    -54.4760987982241, 161.585836858041, -155.698979859887, 66.8013118877197, -13.2806815528857,
  ]
  let c = [
    -0.00778489400243029, -0.322396458041136, -2.40075827716184, -2.54973253934373,
    4.37466414146497, 2.93816398269878,
  ]
  let d = [
    0.00778469570904146, 0.32246712907004, 2.445134137143, 3.75440866190742,
  ]
  if q < 0.02425 {
    let r = double_sqrt(-2.0 * double_log(q))
    let numerator = ((((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r +
      c[5]
    let denominator = (((d[0] * r + d[1]) * r + d[2]) * r + d[3]) * r + 1.0
    numerator / denominator
  } else if q > 0.97575 {
    let r = double_sqrt(-2.0 * double_log(1.0 - q))
    let numerator = ((((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r +
      c[5]
    let denominator = (((d[0] * r + d[1]) * r + d[2]) * r + d[3]) * r + 1.0
    -numerator / denominator
  } else {
    let r = q - 0.5
    let s = r * r
    let numerator = (
        ((((a[0] * s + a[1]) * s + a[2]) * s + a[3]) * s + a[4]) * s + a[5]
      ) *
      r
    let denominator = ((((b[0] * s + b[1]) * s + b[2]) * s + b[3]) * s + b[4]) *
      s +
      1.0
    numerator / denominator
  }
}
