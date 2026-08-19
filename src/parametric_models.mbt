///|
/// Supported fully parametric survival distributions.
pub(all) enum ParametricFamily {
  ExponentialFamily
  WeibullFamily
  LogNormalFamily
  GompertzFamily
} derive(Eq, Debug)

///|
/// A fitted parametric survival model with a common prediction interface.
pub(all) struct ParametricFit {
  family : ParametricFamily
  parameters : Array[Double]
  log_likelihood : Double
  aic : Double
  bic : Double
  sample_size : Int
  event_count : Int
  iterations : Int
  converged : Bool
} derive(Eq, Debug)

///|
/// One point in a parametric survival curve.
pub(all) struct ParametricCurvePoint {
  time : Double
  survival : Double
  hazard : Double
  density : Double
  cumulative_hazard : Double
} derive(Eq, Debug)

///|
/// A confidence interval on a model parameter.
pub(all) struct ParametricInterval {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
} derive(Eq, Debug)

///|
/// Compare fitted models by an information criterion.
pub(all) struct ParametricModelComparison {
  family : ParametricFamily
  log_likelihood : Double
  aic : Double
  bic : Double
  delta_aic : Double
  weight : Double
} derive(Eq, Debug)

///|
/// Return a stable human-readable distribution name.
pub fn parametric_family_name(family : ParametricFamily) -> String {
  match family {
    ExponentialFamily => "exponential"
    WeibullFamily => "weibull"
    LogNormalFamily => "log-normal"
    GompertzFamily => "gompertz"
  }
}

///|
/// Return the number of free parameters for a family.
pub fn parametric_parameter_count(family : ParametricFamily) -> Int {
  match family {
    ExponentialFamily => 1
    WeibullFamily => 2
    LogNormalFamily => 2
    GompertzFamily => 2
  }
}

///|
/// Validate the positive-time portion of a survival data set.
pub fn parametric_data_is_valid(data : SurvivalData) -> Bool {
  for observation in data.to_array() {
    if observation.time < 0.0 ||
      observation.time.is_nan() ||
      observation.time.is_inf() {
      return false
    }
  }
  true
}

///|
/// Calculate a censored exponential log likelihood.
pub fn parametric_exponential_log_likelihood(
  data : SurvivalData,
  rate : Double,
) -> Double {
  if rate <= 0.0 || !parametric_data_is_valid(data) {
    return double_log(0.0)
  }
  let mut result = 0.0
  for observation in data.to_array() {
    result += observation.event.to_double() * double_log(rate) -
      rate * observation.time
  }
  result
}

///|
/// Calculate a censored Weibull log likelihood.
pub fn parametric_weibull_log_likelihood(
  data : SurvivalData,
  shape : Double,
  scale : Double,
) -> Double {
  if shape <= 0.0 || scale <= 0.0 || !parametric_data_is_valid(data) {
    return double_log(0.0)
  }
  let mut result = 0.0
  for observation in data.to_array() {
    let scaled = @math.pow(clamp_positive(observation.time) / scale, shape)
    if observation.event == 1 {
      result += double_log(shape) -
        double_log(scale) +
        (shape - 1.0) * double_log(clamp_positive(observation.time) / scale) -
        scaled
    } else {
      result -= scaled
    }
  }
  result
}

///|
/// Calculate a censored log-normal log likelihood.
pub fn parametric_lognormal_log_likelihood(
  data : SurvivalData,
  location : Double,
  scale : Double,
) -> Double {
  if scale <= 0.0 || !parametric_data_is_valid(data) {
    return double_log(0.0)
  }
  let mut result = 0.0
  for observation in data.to_array() {
    let log_time = double_log(clamp_positive(observation.time))
    let standardized = (log_time - location) / scale
    if observation.event == 1 {
      result += double_log(normal_density(standardized)) -
        double_log(clamp_positive(observation.time)) -
        double_log(scale)
    } else {
      result += double_log(normal_cdf(-standardized))
    }
  }
  result
}

///|
/// Calculate a censored Gompertz log likelihood.
pub fn parametric_gompertz_log_likelihood(
  data : SurvivalData,
  shape : Double,
  rate : Double,
) -> Double {
  if rate <= 0.0 || !parametric_data_is_valid(data) {
    return double_log(0.0)
  }
  let mut result = 0.0
  for observation in data.to_array() {
    let time = observation.time
    let cumulative = if shape.abs() <= survival_epsilon {
      rate * time
    } else {
      rate / shape * (double_exp(shape * time) - 1.0)
    }
    if observation.event == 1 {
      result += double_log(rate) + shape * time - cumulative
    } else {
      result -= cumulative
    }
  }
  result
}

///|
/// Fit an exponential model and expose it through the common fit record.
pub fn fit_parametric_exponential(data : SurvivalData) -> ParametricFit {
  let fit = exponential_fit(data)
  let parameters = [fit.rate]
  let log_likelihood = fit.log_likelihood
  parametric_make_fit(
    ExponentialFamily,
    parameters,
    log_likelihood,
    data,
    1,
    fit.rate > 0.0 || data.event_count() == 0,
  )
}

///|
/// Fit a Weibull model using the existing censored likelihood solver.
pub fn fit_parametric_weibull(data : SurvivalData) -> ParametricFit {
  let fit = weibull_fit(data)
  parametric_make_fit(
    WeibullFamily,
    [fit.shape, fit.scale],
    fit.log_likelihood,
    data,
    fit.iterations,
    fit.converged || data.event_count() == 0,
  )
}

///|
/// Fit a log-normal model with a censored EM-style location/scale update.
pub fn fit_parametric_lognormal(
  data : SurvivalData,
  max_iterations? : Int = 60,
) -> ParametricFit {
  if data.length() == 0 || !parametric_data_is_valid(data) {
    return parametric_make_fit(LogNormalFamily, [0.0, 1.0], 0.0, data, 0, false)
  }
  let mut location = parametric_log_time_mean(data)
  let mut scale = parametric_log_time_scale(data, location)
  let mut converged = false
  let mut iterations = 0
  let limit = if max_iterations <= 0 { 1 } else { max_iterations }
  for iteration in 0..<limit {
    let moments = parametric_lognormal_expected_moments(data, location, scale)
    let next_location = moments[0]
    let next_scale = double_sqrt(moments[1])
    let movement = (next_location - location).abs() + (next_scale - scale).abs()
    location = next_location
    scale = clamp_positive(next_scale, floor=0.000001)
    iterations = iteration + 1
    if movement <= 0.000001 {
      converged = true
      break
    }
  }
  let log_likelihood = parametric_lognormal_log_likelihood(
    data, location, scale,
  )
  parametric_make_fit(
    LogNormalFamily,
    [location, scale],
    log_likelihood,
    data,
    iterations,
    converged,
  )
}

///|
/// Fit a Gompertz model by a stable two-parameter grid/Newton procedure.
pub fn fit_parametric_gompertz(
  data : SurvivalData,
  max_iterations? : Int = 50,
) -> ParametricFit {
  let exponential = exponential_fit(data)
  if data.length() == 0 || exponential.rate <= 0.0 {
    return parametric_make_fit(
      GompertzFamily,
      [0.0, exponential.rate],
      0.0,
      data,
      0,
      false,
    )
  }
  let mut shape = 0.0
  let mut rate = exponential.rate
  let mut converged = false
  let mut iterations = 0
  let limit = if max_iterations <= 0 { 1 } else { max_iterations }
  for iteration in 0..<limit {
    let step = 0.0001
    let left = parametric_gompertz_log_likelihood(data, shape - step, rate)
    let center = parametric_gompertz_log_likelihood(data, shape, rate)
    let right = parametric_gompertz_log_likelihood(data, shape + step, rate)
    let score = (right - left) / (2.0 * step)
    let curvature = (right - 2.0 * center + left) / (step * step)
    if curvature.abs() <= survival_epsilon {
      break
    }
    let update = score / curvature
    let candidate = shape - update
    shape = if candidate < -2.0 {
      -2.0
    } else if candidate > 2.0 {
      2.0
    } else {
      candidate
    }
    let cumulative = parametric_gompertz_exposure(data, shape)
    rate = if cumulative <= survival_epsilon {
      exponential.rate
    } else {
      data.event_count().to_double() / cumulative
    }
    iterations = iteration + 1
    if update.abs() <= 0.000001 {
      converged = true
      break
    }
  }
  parametric_make_fit(
    GompertzFamily,
    [shape, rate],
    parametric_gompertz_log_likelihood(data, shape, rate),
    data,
    iterations,
    converged,
  )
}

///|
/// Fit all supported distributions for model screening.
pub fn fit_parametric_families(data : SurvivalData) -> Array[ParametricFit] {
  [
    fit_parametric_exponential(data),
    fit_parametric_weibull(data),
    fit_parametric_lognormal(data),
    fit_parametric_gompertz(data),
  ]
}

///|
/// Predict survival under a fitted parametric model.
pub fn ParametricFit::survival_at(
  self : ParametricFit,
  time : Double,
) -> Double {
  if time <= 0.0 {
    return 1.0
  }
  match self.family {
    ExponentialFamily => double_exp(-self.parameters[0] * time)
    WeibullFamily => {
      let shape = self.parameters[0]
      let scale = self.parameters[1]
      if scale <= 0.0 {
        1.0
      } else {
        double_exp(-@math.pow(time / scale, shape))
      }
    }
    LogNormalFamily => {
      let location = self.parameters[0]
      let scale = self.parameters[1]
      normal_cdf(-(double_log(time) - location) / scale)
    }
    GompertzFamily => {
      let shape = self.parameters[0]
      let rate = self.parameters[1]
      let cumulative = if shape.abs() <= survival_epsilon {
        rate * time
      } else {
        rate / shape * (double_exp(shape * time) - 1.0)
      }
      double_exp(-cumulative)
    }
  }
}

///|
/// Predict cumulative hazard under a fitted model.
pub fn ParametricFit::cumulative_hazard_at(
  self : ParametricFit,
  time : Double,
) -> Double {
  -double_log(self.survival_at(time))
}

///|
/// Predict the instantaneous hazard.
pub fn ParametricFit::hazard_at(self : ParametricFit, time : Double) -> Double {
  if time <= 0.0 {
    return 0.0
  }
  match self.family {
    ExponentialFamily => self.parameters[0]
    WeibullFamily => {
      let shape = self.parameters[0]
      let scale = self.parameters[1]
      if scale <= 0.0 {
        0.0
      } else {
        shape / scale * @math.pow(time / scale, shape - 1.0)
      }
    }
    LogNormalFamily => {
      let location = self.parameters[0]
      let scale = self.parameters[1]
      let z = (double_log(time) - location) / scale
      normal_density(z) / (time * scale * self.survival_at(time))
    }
    GompertzFamily => self.parameters[1] * double_exp(self.parameters[0] * time)
  }
}

///|
/// Predict the density of the event time.
pub fn ParametricFit::density_at(self : ParametricFit, time : Double) -> Double {
  self.survival_at(time) * self.hazard_at(time)
}

///|
/// Predict a failure-time quantile.
pub fn ParametricFit::quantile(
  self : ParametricFit,
  probability : Double,
) -> Double {
  let p = clamp_probability(probability)
  match self.family {
    ExponentialFamily =>
      if self.parameters[0] <= 0.0 {
        0.0
      } else {
        -double_log(1.0 - p) / self.parameters[0]
      }
    WeibullFamily =>
      if self.parameters[1] <= 0.0 {
        0.0
      } else {
        self.parameters[1] *
        @math.pow(-double_log(1.0 - p), 1.0 / self.parameters[0])
      }
    LogNormalFamily =>
      double_exp(self.parameters[0] + self.parameters[1] * normal_quantile(p))
    GompertzFamily => {
      let shape = self.parameters[0]
      let rate = self.parameters[1]
      if rate <= 0.0 {
        0.0
      } else if shape.abs() <= survival_epsilon {
        -double_log(1.0 - p) / rate
      } else {
        double_log(1.0 + shape * -double_log(1.0 - p) / rate) / shape
      }
    }
  }
}

///|
/// Numerically integrate survival to obtain a restricted mean.
pub fn ParametricFit::restricted_mean(
  self : ParametricFit,
  horizon : Double,
  points? : Int = 256,
) -> Double {
  if horizon <= 0.0 {
    return 0.0
  }
  let count = if points < 2 { 2 } else { points }
  let width = horizon / (count - 1).to_double()
  let mut total = 0.0
  for i in 0..<count {
    let value = self.survival_at(i.to_double() * width)
    let coefficient = if i == 0 || i == count - 1 { 0.5 } else { 1.0 }
    total += coefficient * value
  }
  total * width
}

///|
/// Generate a curve over a supplied time grid.
pub fn ParametricFit::curve(
  self : ParametricFit,
  times : Array[Double],
) -> Array[ParametricCurvePoint] {
  let result : Array[ParametricCurvePoint] = []
  for time in times {
    result.push({
      time,
      survival: self.survival_at(time),
      hazard: self.hazard_at(time),
      density: self.density_at(time),
      cumulative_hazard: self.cumulative_hazard_at(time),
    })
  }
  result
}

///|
/// Build a parameter interval using a supplied standard error and confidence.
pub fn parametric_interval(
  estimate : Double,
  standard_error : Double,
  confidence? : Double = 0.95,
) -> ParametricInterval {
  let level = if confidence <= 0.0 || confidence >= 1.0 {
    0.95
  } else {
    confidence
  }
  let z = normal_quantile(0.5 + level / 2.0)
  {
    estimate,
    standard_error,
    lower: estimate - z * standard_error,
    upper: estimate + z * standard_error,
  }
}

///|
/// Select the best model by AIC and return delta-AIC weights.
pub fn parametric_compare_by_aic(
  fits : Array[ParametricFit],
) -> Array[ParametricModelComparison] {
  if fits.length() == 0 {
    return []
  }
  let mut minimum = fits[0].aic
  for fit in fits {
    if fit.aic < minimum {
      minimum = fit.aic
    }
  }
  let mut denominator = 0.0
  for fit in fits {
    denominator += double_exp(-0.5 * (fit.aic - minimum))
  }
  let result : Array[ParametricModelComparison] = []
  for fit in fits {
    let weight = if denominator <= survival_epsilon {
      0.0
    } else {
      double_exp(-0.5 * (fit.aic - minimum)) / denominator
    }
    result.push({
      family: fit.family,
      log_likelihood: fit.log_likelihood,
      aic: fit.aic,
      bic: fit.bic,
      delta_aic: fit.aic - minimum,
      weight,
    })
  }
  result
}

///|
/// Return the fit with the smallest AIC.
pub fn select_parametric_model(fits : Array[ParametricFit]) -> ParametricFit? {
  if fits.length() == 0 {
    return None
  }
  let mut best = fits[0]
  for fit in fits {
    if fit.aic < best.aic {
      best = fit
    }
  }
  Some(best)
}

///|
/// Approximate the expected event time from a fitted model.
pub fn ParametricFit::mean_time(
  self : ParametricFit,
  horizon? : Double = 1000.0,
) -> Double {
  self.restricted_mean(horizon)
}

///|
/// Render model comparisons as CSV for review and audit.
pub fn parametric_comparisons_to_csv(
  comparisons : Array[ParametricModelComparison],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("family,log_likelihood,aic,bic,delta_aic,weight\n")
  for comparison in comparisons {
    builder.write_string(
      "\{parametric_family_name(comparison.family)},\{comparison.log_likelihood},\{comparison.aic},\{comparison.bic},\{comparison.delta_aic},\{comparison.weight}\n",
    )
  }
  builder.to_string()
}

///|
/// Compute a model's Brier score on right-censored observations.
pub fn ParametricFit::brier_score(
  self : ParametricFit,
  data : SurvivalData,
  horizon : Double,
) -> Double {
  let mut total = 0.0
  let mut count = 0
  for observation in data.to_array() {
    let prediction = 1.0 - self.survival_at(horizon)
    let outcome = if observation.time <= horizon && observation.event == 1 {
      1.0
    } else if observation.time > horizon {
      0.0
    } else {
      continue
    }
    total += (prediction - outcome) * (prediction - outcome)
    count += 1
  }
  if count == 0 {
    0.0
  } else {
    total / count.to_double()
  }
}

///|
fn parametric_make_fit(
  family : ParametricFamily,
  parameters : Array[Double],
  log_likelihood : Double,
  data : SurvivalData,
  iterations : Int,
  converged : Bool,
) -> ParametricFit {
  let parameter_count = parameters.length()
  let n = data.length()
  {
    family,
    parameters,
    log_likelihood,
    aic: 2.0 * parameter_count.to_double() - 2.0 * log_likelihood,
    bic: if n <= 0 {
      0.0
    } else {
      parameter_count.to_double() * double_log(n.to_double()) -
      2.0 * log_likelihood
    },
    sample_size: n,
    event_count: data.event_count(),
    iterations,
    converged,
  }
}

///|
fn parametric_log_time_mean(data : SurvivalData) -> Double {
  let values : Array[Double] = []
  for observation in data.to_array() {
    values.push(double_log(clamp_positive(observation.time)))
  }
  mean_doubles(values)
}

///|
fn parametric_log_time_scale(data : SurvivalData, location : Double) -> Double {
  let values : Array[Double] = []
  for observation in data.to_array() {
    values.push(double_log(clamp_positive(observation.time)))
  }
  let mut total = 0.0
  for value in values {
    let difference = value - location
    total += difference * difference
  }
  if values.length() < 2 {
    1.0
  } else {
    double_sqrt(total / (values.length() - 1).to_double())
  }
}

///|
fn parametric_lognormal_expected_moments(
  data : SurvivalData,
  location : Double,
  scale : Double,
) -> Array[Double] {
  let values : Array[Double] = []
  for observation in data.to_array() {
    let observed = double_log(clamp_positive(observation.time))
    if observation.event == 1 {
      values.push(observed)
    } else {
      let z = (observed - location) / scale
      let tail = clamp_positive(normal_cdf(-z), floor=0.000001)
      let inverse_mills = normal_density(z) / tail
      values.push(observed + scale * inverse_mills)
    }
  }
  let next_location = mean_doubles(values)
  let mut variance = 0.0
  for value in values {
    let difference = value - next_location
    variance += difference * difference
  }
  let next_variance = if values.length() < 2 {
    1.0
  } else {
    variance / (values.length() - 1).to_double()
  }
  [next_location, next_variance]
}

///|
fn parametric_gompertz_exposure(data : SurvivalData, shape : Double) -> Double {
  let mut total = 0.0
  for observation in data.to_array() {
    total += if shape.abs() <= survival_epsilon {
      observation.time
    } else {
      (double_exp(shape * observation.time) - 1.0) / shape
    }
  }
  total
}
