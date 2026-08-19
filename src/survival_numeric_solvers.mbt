///|
/// Stable numerical primitives shared by model and reporting code.
pub(all) struct SurvivalLinearFit {
  intercept : Double
  slope : Double
  residual_sum_squares : Double
  r_squared : Double
  sample_size : Int
  solver_row : Int
} derive(Eq, Debug)

///|
pub(all) struct SurvivalRootResult {
  root : Double
  value_at_root : Double
  iterations : Int
  converged : Bool
  solver_row : Int
} derive(Eq, Debug)

///|
/// Stable log(sum(exp(values))).
pub fn survival_numeric_logsumexp(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let maximum = survival_ext_max(values)
  if maximum.is_inf() || maximum.is_nan() {
    return maximum
  }
  let mut total = 0.0
  for value in values {
    total += double_exp(value - maximum)
  }
  maximum + double_log(total)
}

///|
/// Stable softmax probabilities.
pub fn survival_numeric_softmax(values : Array[Double]) -> Array[Double] {
  let output : Array[Double] = []
  if values.length() == 0 {
    return output
  }
  let normalizer = survival_numeric_logsumexp(values)
  for value in values {
    output.push(double_exp(value - normalizer))
  }
  output
}

///|
pub fn survival_numeric_logit(probability : Double) -> Double {
  let value = clamp_probability(probability)
  if value <= 0.0 {
    -37.0
  } else if value >= 1.0 {
    37.0
  } else {
    double_log(value / (1.0 - value))
  }
}

///|
pub fn survival_numeric_expit(value : Double) -> Double {
  if value >= 0.0 {
    1.0 / (1.0 + double_exp(-value))
  } else {
    let exponential = double_exp(value)
    exponential / (1.0 + exponential)
  }
}

///|
/// Ordinary least squares fit for one covariate.
pub fn survival_numeric_linear_fit(
  x : Array[Double],
  y : Array[Double],
) -> SurvivalLinearFit {
  let count = x.length().min(y.length())
  if count == 0 {
    return {
      intercept: 0.0,
      slope: 0.0,
      residual_sum_squares: 0.0,
      r_squared: 0.0,
      sample_size: 0,
      solver_row: 0,
    }
  }
  let x_values : Array[Double] = []
  let y_values : Array[Double] = []
  for i = 0; i < count; i = i + 1 {
    x_values.push(x[i])
    y_values.push(y[i])
  }
  let mean_x = survival_ext_mean(x_values)
  let mean_y = survival_ext_mean(y_values)
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i = 0; i < count; i = i + 1 {
    numerator += (x[i] - mean_x) * (y[i] - mean_y)
    denominator += (x[i] - mean_x) * (x[i] - mean_x)
  }
  let slope = if denominator <= survival_epsilon {
    0.0
  } else {
    numerator / denominator
  }
  let intercept = mean_y - slope * mean_x
  let mut residual = 0.0
  let mut total = 0.0
  for i = 0; i < count; i = i + 1 {
    let error = y[i] - (intercept + slope * x[i])
    residual += error * error
    let centered = y[i] - mean_y
    total += centered * centered
  }
  let r_squared = if total <= survival_epsilon {
    0.0
  } else {
    1.0 - residual / total
  }
  {
    intercept,
    slope,
    residual_sum_squares: residual,
    r_squared,
    sample_size: count,
    solver_row: count,
  }
}

///|
pub fn survival_numeric_linear_predict(
  fit : SurvivalLinearFit,
  x : Double,
) -> Double {
  fit.intercept + fit.slope * x
}

///|
/// Linear interpolation on a sorted grid, with endpoint extrapolation.
pub fn survival_numeric_interpolate(
  x : Array[Double],
  y : Array[Double],
  query : Double,
) -> Double {
  let count = x.length().min(y.length())
  if count == 0 {
    return 0.0
  }
  if count == 1 {
    return y[0]
  }
  if query <= x[0] {
    return y[0]
  }
  for i = 1; i < count; i = i + 1 {
    if query <= x[i] {
      let width = x[i] - x[i - 1]
      if width.abs() <= survival_epsilon {
        return y[i]
      }
      let fraction = (query - x[i - 1]) / width
      return y[i - 1] + fraction * (y[i] - y[i - 1])
    }
  }
  y[count - 1]
}

///|
/// Central finite difference of a scalar sequence on an equally spaced grid.
pub fn survival_numeric_derivative(
  values : Array[Double],
  spacing : Double,
) -> Array[Double] {
  let output : Array[Double] = []
  if values.length() == 0 || spacing == 0.0 {
    return output
  }
  for i = 0; i < values.length(); i = i + 1 {
    let derivative = if i == 0 {
      (values[if values.length() > 1 { 1 } else { 0 }] - values[0]) / spacing
    } else if i + 1 >= values.length() {
      (values[i] - values[i - 1]) / spacing
    } else {
      (values[i + 1] - values[i - 1]) / (2.0 * spacing)
    }
    output.push(derivative)
  }
  output
}

///|
/// Solve a two-by-two linear system, returning None for a singular matrix.
pub fn survival_numeric_solve_2x2(
  a11 : Double,
  a12 : Double,
  a21 : Double,
  a22 : Double,
  b1 : Double,
  b2 : Double,
) -> (Double, Double, Bool) {
  let determinant = a11 * a22 - a12 * a21
  if determinant.abs() <= survival_epsilon {
    (0.0, 0.0, false)
  } else {
    (
      (b1 * a22 - a12 * b2) / determinant,
      (a11 * b2 - b1 * a21) / determinant,
      true,
    )
  }
}

///|
pub fn survival_numeric_linear_fit_is_valid(fit : SurvivalLinearFit) -> Bool {
  fit.sample_size >= 0 &&
  fit.residual_sum_squares >= 0.0 &&
  fit.r_squared <= 1.0
}

///|
pub fn survival_numeric_root_is_valid(root : SurvivalRootResult) -> Bool {
  root.iterations >= 0 && root.value_at_root.is_nan() == false
}
