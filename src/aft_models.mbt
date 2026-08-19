///|
/// Accelerated-failure-time distributions supported by the regression layer.
/// The coefficient vector is interpreted on log time: a positive coefficient
/// delays the event when the corresponding covariate increases.
pub(all) enum AftDistribution {
  AftLogNormal
  AftLogLogistic
  AftWeibull
} derive(Eq, Debug)

///|
/// Validation failures for censored AFT input.
pub(all) enum AftError {
  EmptyData
  LengthMismatch
  NegativeTime
  InvalidTime
  InvalidEvent
  CovariateLengthMismatch
  InvalidCovariate
  InvalidIteration
} derive(Eq, Debug)

///|
/// One right-censored observation for an AFT model.
pub(all) struct AftObservation {
  time : Double
  event : Int
  covariates : Array[Double]
  aft_row : Bool
} derive(Eq, Debug)

///|
/// Sorted AFT input with a fixed-width design matrix.
pub(all) struct AftData {
  observations : Array[AftObservation]
  n_covariates : Int
} derive(Eq, Debug)

///|
/// Fitted AFT regression model.
pub(all) struct AftFit {
  distribution : AftDistribution
  coefficients : Array[Double]
  location : Double
  scale : Double
  log_likelihood : Double
  iterations : Int
  converged : Bool
  sample_size : Int
  event_count : Int
} derive(Eq, Debug)

///|
/// A prediction row produced by an AFT model.
pub(all) struct AftPrediction {
  time : Double
  survival : Double
  cumulative_hazard : Double
  hazard : Double
  density : Double
} derive(Eq, Debug)

///|
/// A coefficient summary with a time-ratio interpretation.
pub(all) struct AftCoefficientSummary {
  index : Int
  coefficient : Double
  time_ratio : Double
  standard_error : Double
  z_value : Double
  p_value : Double
} derive(Eq, Debug)

///|
/// Fit summaries for comparing candidate AFT distributions.
pub(all) struct AftModelScore {
  distribution : AftDistribution
  log_likelihood : Double
  aic : Double
  bic : Double
  sample_size : Int
  converged : Bool
} derive(Eq, Debug)

///|
/// Create and validate AFT data.
pub fn AftData::try_new(
  times : Array[Double],
  events : Array[Int],
  covariates : Array[Array[Double]],
) -> Result[AftData, AftError] {
  if times.length() == 0 {
    return Err(EmptyData)
  }
  if times.length() != events.length() || times.length() != covariates.length() {
    return Err(LengthMismatch)
  }
  let width = covariates[0].length()
  let observations : Array[AftObservation] = []
  for i in 0..<times.length() {
    let time = times[i]
    if time.is_nan() || time.is_inf() {
      return Err(InvalidTime)
    }
    if time <= 0.0 {
      return Err(NegativeTime)
    }
    if events[i] != 0 && events[i] != 1 {
      return Err(InvalidEvent)
    }
    if covariates[i].length() != width {
      return Err(CovariateLengthMismatch)
    }
    for value in covariates[i] {
      if value.is_nan() || value.is_inf() {
        return Err(InvalidCovariate)
      }
    }
    let row : Array[Double] = []
    for value in covariates[i] {
      row.push(value)
    }
    let observation : AftObservation = {
      time,
      event: events[i],
      covariates: row,
      aft_row: true,
    }
    observations.push(observation)
  }
  observations.sort_by(fn(a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      b.event.compare(a.event)
    }
  })
  Ok({ observations, n_covariates: width })
}

///|
/// Convenience constructor that preserves the historical panic-on-invalid
/// behavior used by small interactive examples.
pub fn AftData::new(
  times : Array[Double],
  events : Array[Int],
  covariates : Array[Array[Double]],
) -> AftData {
  match AftData::try_new(times, events, covariates) {
    Ok(data) => data
    Err(_) => panic()
  }
}

///|
/// Number of observations.
pub fn AftData::length(self : AftData) -> Int {
  self.observations.length()
}

///|
/// Number of covariates, excluding the implicit intercept/location.
pub fn AftData::covariate_count(self : AftData) -> Int {
  self.n_covariates
}

///|
/// Access an observation in sorted order.
pub fn AftData::get(self : AftData, index : Int) -> AftObservation {
  self.observations[index]
}

///|
/// Copy the response vector.
pub fn AftData::times(self : AftData) -> Array[Double] {
  let result : Array[Double] = []
  for observation in self.observations {
    result.push(observation.time)
  }
  result
}

///|
/// Copy event indicators.
pub fn AftData::events(self : AftData) -> Array[Int] {
  let result : Array[Int] = []
  for observation in self.observations {
    result.push(observation.event)
  }
  result
}

///|
/// Return a stable display name for a distribution.
pub fn aft_distribution_name(distribution : AftDistribution) -> String {
  match distribution {
    AftLogNormal => "lognormal"
    AftLogLogistic => "log-logistic"
    AftWeibull => "weibull"
  }
}

///|
/// Evaluate the log survival function for a standardized log time.
fn aft_log_survival(
  distribution : AftDistribution,
  standardized : Double,
) -> Double {
  match distribution {
    AftLogNormal => double_log(normal_cdf(-standardized))
    AftLogLogistic => -double_log(1.0 + double_exp(standardized))
    AftWeibull => -double_exp(standardized)
  }
}

///|
/// Evaluate the log density with respect to log time.
fn aft_log_density(
  distribution : AftDistribution,
  standardized : Double,
  scale : Double,
) -> Double {
  let safe_scale = clamp_positive(scale)
  match distribution {
    AftLogNormal => {
      let log_two_pi = 1.8378770664093453
      -0.5 * standardized * standardized -
      0.5 * log_two_pi -
      double_log(safe_scale)
    }
    AftLogLogistic => {
      let exp_z = double_exp(standardized)
      double_log(exp_z) - 2.0 * double_log(1.0 + exp_z) - double_log(safe_scale)
    }
    AftWeibull =>
      standardized - double_exp(standardized) - double_log(safe_scale)
  }
}

///|
/// Evaluate the log likelihood for an AFT parameter vector.
pub fn aft_log_likelihood(
  data : AftData,
  distribution : AftDistribution,
  coefficients : Array[Double],
  location : Double,
  scale : Double,
) -> Double {
  if coefficients.length() != data.covariate_count() || scale <= 0.0 {
    return -1.0e300
  }
  let mut value = 0.0
  for observation in data.observations {
    let linear = dot_product(observation.covariates, coefficients)
    let standardized = (double_log(observation.time) - location - linear) /
      scale
    if observation.event == 1 {
      value += aft_log_density(distribution, standardized, scale) -
        double_log(observation.time)
    } else {
      value += aft_log_survival(distribution, standardized)
    }
  }
  value
}

///|
/// Return the median log time used for a robust starting location.
fn aft_initial_location(data : AftData) -> Double {
  let values : Array[Double] = []
  for observation in data.observations {
    values.push(double_log(observation.time))
  }
  values.sort()
  let middle = values.length() / 2
  if values.length() % 2 == 1 {
    values[middle]
  } else {
    0.5 * (values[middle - 1] + values[middle])
  }
}

///|
/// Estimate a robust scale from log times around the initial location.
fn aft_initial_scale(data : AftData, location : Double) -> Double {
  let mut total = 0.0
  for observation in data.observations {
    let residual = double_log(observation.time) - location
    total += residual * residual
  }
  clamp_positive(double_sqrt(total / data.length().to_double()))
}

///|
/// Numerically differentiate an AFT likelihood component.
fn aft_numeric_score(
  data : AftData,
  distribution : AftDistribution,
  coefficients : Array[Double],
  location : Double,
  scale : Double,
  index : Int,
) -> Double {
  let step = 0.0001
  let plus = []
  let minus = []
  for coefficient in coefficients {
    plus.push(coefficient)
    minus.push(coefficient)
  }
  plus[index] += step
  minus[index] -= step
  (
    aft_log_likelihood(data, distribution, plus, location, scale) -
    aft_log_likelihood(data, distribution, minus, location, scale)
  ) /
  (2.0 * step)
}

///|
/// Fit a censored AFT model with deterministic bounded coordinate updates.
/// The routine is deliberately derivative-free at the boundary and is useful
/// for small clinical datasets where a failed Newton step should not abort a
/// pipeline.
pub fn fit_aft(
  data : AftData,
  distribution : AftDistribution,
  max_iterations? : Int = 80,
  tolerance? : Double = 0.000001,
) -> Result[AftFit, AftError] {
  if max_iterations <= 0 {
    return Err(InvalidIteration)
  }
  let coefficients = zeros(data.covariate_count())
  let mut location = aft_initial_location(data)
  let mut scale = aft_initial_scale(data, location)
  let mut current = aft_log_likelihood(
    data, distribution, coefficients, location, scale,
  )
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<max_iterations {
    let mut changed = 0.0
    for index in 0..<coefficients.length() {
      let score = aft_numeric_score(
        data, distribution, coefficients, location, scale, index,
      )
      let step = clamp_positive(
        0.02 / (1.0 + iteration.to_double() * 0.1),
        floor=0.0001,
      )
      let proposal = coefficients[index] +
        step * score / data.length().to_double()
      let old = coefficients[index]
      coefficients[index] = if proposal > 20.0 {
        20.0
      } else if proposal < -20.0 {
        -20.0
      } else {
        proposal
      }
      changed += (coefficients[index] - old).abs()
    }
    let location_step = 0.01 / (1.0 + iteration.to_double() * 0.1)
    let location_score = (
        aft_log_likelihood(
          data,
          distribution,
          coefficients,
          location + 0.0001,
          scale,
        ) -
        aft_log_likelihood(
          data,
          distribution,
          coefficients,
          location - 0.0001,
          scale,
        )
      ) /
      0.0002
    location += location_step * location_score / data.length().to_double()
    let scale_score = (
        aft_log_likelihood(
          data,
          distribution,
          coefficients,
          location,
          scale + 0.0001,
        ) -
        aft_log_likelihood(
          data,
          distribution,
          coefficients,
          location,
          clamp_positive(scale - 0.0001),
        )
      ) /
      0.0002
    scale = clamp_positive(
      scale + location_step * scale_score / data.length().to_double(),
    )
    let proposed = aft_log_likelihood(
      data, distribution, coefficients, location, scale,
    )
    if proposed.is_nan() || proposed.is_inf() {
      break
    }
    if (proposed - current).abs() <= tolerance && changed <= tolerance {
      converged = true
      iterations = iteration + 1
      current = proposed
      break
    }
    current = proposed
    iterations = iteration + 1
  }
  let mut event_count = 0
  for observation in data.observations {
    event_count += observation.event
  }
  Ok({
    distribution,
    coefficients,
    location,
    scale,
    log_likelihood: current,
    iterations,
    converged,
    sample_size: data.length(),
    event_count,
  })
}

///|
/// Fit all supported distributions for model selection.
pub fn fit_aft_families(
  data : AftData,
  max_iterations? : Int = 80,
) -> Array[AftFit] {
  let result : Array[AftFit] = []
  for distribution in [AftLogNormal, AftLogLogistic, AftWeibull] {
    match fit_aft(data, distribution, max_iterations~) {
      Ok(fit) => result.push(fit)
      Err(_) => ()
    }
  }
  result
}

///|
/// Return the linear predictor for a covariate vector.
pub fn AftFit::linear_predictor(
  self : AftFit,
  covariates : Array[Double],
) -> Double {
  if covariates.length() != self.coefficients.length() {
    return 0.0
  }
  dot_product(covariates, self.coefficients)
}

///|
/// Log time location conditional on covariates.
pub fn AftFit::conditional_location(
  self : AftFit,
  covariates : Array[Double],
) -> Double {
  self.location + self.linear_predictor(covariates)
}

///|
/// Predict survival at a time and covariate vector.
pub fn AftFit::survival_at(
  self : AftFit,
  time : Double,
  covariates : Array[Double],
) -> Double {
  if time <= 0.0 || covariates.length() != self.coefficients.length() {
    return 1.0
  }
  let z = (double_log(time) - self.conditional_location(covariates)) /
    clamp_positive(self.scale)
  clamp_probability(double_exp(aft_log_survival(self.distribution, z)))
}

///|
/// Predict the cumulative hazard.
pub fn AftFit::cumulative_hazard_at(
  self : AftFit,
  time : Double,
  covariates : Array[Double],
) -> Double {
  let survival = self.survival_at(time, covariates)
  -double_log(clamp_positive(survival))
}

///|
/// Predict a density on the original time scale.
pub fn AftFit::density_at(
  self : AftFit,
  time : Double,
  covariates : Array[Double],
) -> Double {
  if time <= 0.0 || covariates.length() != self.coefficients.length() {
    return 0.0
  }
  let z = (double_log(time) - self.conditional_location(covariates)) /
    clamp_positive(self.scale)
  let log_density = aft_log_density(self.distribution, z, self.scale) -
    double_log(time)
  double_exp(log_density)
}

///|
/// Predict the hazard on the original time scale.
pub fn AftFit::hazard_at(
  self : AftFit,
  time : Double,
  covariates : Array[Double],
) -> Double {
  let survival = self.survival_at(time, covariates)
  if survival <= survival_epsilon {
    0.0
  } else {
    self.density_at(time, covariates) / survival
  }
}

///|
/// Predict a quantile of event time.
pub fn AftFit::quantile(
  self : AftFit,
  probability : Double,
  covariates : Array[Double],
) -> Double {
  let p = clamp_probability(probability)
  if covariates.length() != self.coefficients.length() {
    return 0.0
  }
  let location = self.conditional_location(covariates)
  let time_scale = clamp_positive(self.scale)
  match self.distribution {
    AftLogNormal => double_exp(location + time_scale * normal_quantile(p))
    AftLogLogistic =>
      if p >= 1.0 {
        double_exp(location + 700.0)
      } else {
        double_exp(
          location + time_scale * double_log(p / clamp_positive(1.0 - p)),
        )
      }
    AftWeibull =>
      if p >= 1.0 {
        double_exp(location + 700.0)
      } else {
        double_exp(location + time_scale * double_log(-double_log(1.0 - p)))
      }
  }
}

///|
/// Produce predictions at a regular time grid.
pub fn AftFit::curve(
  self : AftFit,
  times : Array[Double],
  covariates : Array[Double],
) -> Array[AftPrediction] {
  let result : Array[AftPrediction] = []
  for time in times {
    let survival = self.survival_at(time, covariates)
    result.push({
      time,
      survival,
      cumulative_hazard: -double_log(clamp_positive(survival)),
      hazard: self.hazard_at(time, covariates),
      density: self.density_at(time, covariates),
    })
  }
  result
}

///|
/// Compute Cox-Snell-like residuals for censored observations.
pub fn AftFit::residuals(self : AftFit, data : AftData) -> Array[Double] {
  let result : Array[Double] = []
  for observation in data.observations {
    result.push(
      self.cumulative_hazard_at(observation.time, observation.covariates),
    )
  }
  result
}

///|
/// Estimate coefficient uncertainty from a diagonal finite-difference Hessian.
pub fn AftFit::coefficient_summary(
  self : AftFit,
  data : AftData,
) -> Array[AftCoefficientSummary] {
  let result : Array[AftCoefficientSummary] = []
  for index in 0..<self.coefficients.length() {
    let step = 0.001
    let plus = []
    let minus = []
    for coefficient in self.coefficients {
      plus.push(coefficient)
      minus.push(coefficient)
    }
    plus[index] += step
    minus[index] -= step
    let curvature = aft_log_likelihood(
        data,
        self.distribution,
        plus,
        self.location,
        self.scale,
      ) -
      2.0 * self.log_likelihood +
      aft_log_likelihood(
        data,
        self.distribution,
        minus,
        self.location,
        self.scale,
      )
    let variance = if curvature < -survival_epsilon {
      -1.0 / curvature
    } else {
      0.0
    }
    let standard_error = double_sqrt(variance)
    let z = if standard_error > survival_epsilon {
      self.coefficients[index] / standard_error
    } else {
      0.0
    }
    result.push({
      index,
      coefficient: self.coefficients[index],
      time_ratio: double_exp(self.coefficients[index]),
      standard_error,
      z_value: z,
      p_value: 2.0 * (1.0 - normal_cdf(z.abs())),
    })
  }
  result
}

///|
/// Compute AIC and BIC for a fitted AFT model.
pub fn aft_model_score(fit : AftFit) -> AftModelScore {
  let parameters = fit.coefficients.length() + 2
  {
    distribution: fit.distribution,
    log_likelihood: fit.log_likelihood,
    aic: 2.0 * parameters.to_double() - 2.0 * fit.log_likelihood,
    bic: parameters.to_double() * double_log(fit.sample_size.to_double()) -
    2.0 * fit.log_likelihood,
    sample_size: fit.sample_size,
    converged: fit.converged,
  }
}

///|
/// Score candidate fits in AIC order.
pub fn aft_model_scores(fits : Array[AftFit]) -> Array[AftModelScore] {
  let result : Array[AftModelScore] = []
  for fit in fits {
    result.push(aft_model_score(fit))
  }
  result.sort_by(fn(a, b) {
    if a.aic < b.aic {
      -1
    } else if a.aic > b.aic {
      1
    } else {
      0
    }
  })
  result
}

///|
/// Serialize model comparison output for experiment tracking.
pub fn aft_model_scores_to_csv(scores : Array[AftModelScore]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("distribution,log_likelihood,aic,bic,n,converged\n")
  for score in scores {
    builder.write_string(aft_distribution_name(score.distribution))
    builder.write_char(',')
    builder.write_string(score.log_likelihood.to_string())
    builder.write_char(',')
    builder.write_string(score.aic.to_string())
    builder.write_char(',')
    builder.write_string(score.bic.to_string())
    builder.write_char(',')
    builder.write_string(score.sample_size.to_string())
    builder.write_char(',')
    builder.write_string(score.converged.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Return the best converged fit by AIC, or None when no fit is available.
pub fn select_aft_model(fits : Array[AftFit]) -> AftFit? {
  let mut best : AftFit? = None
  for fit in fits {
    if !fit.converged {
      continue
    }
    match best {
      None => best = Some(fit)
      Some(current) =>
        if aft_model_score(fit).aic < aft_model_score(current).aic {
          best = Some(fit)
        }
    }
  }
  best
}

///|
/// Compute a censored integrated Brier score on an evaluation grid.
pub fn aft_integrated_brier_score(
  fit : AftFit,
  data : AftData,
  grid : Array[Double],
) -> Double {
  if grid.length() == 0 || data.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  let mut count = 0
  for time in grid {
    for observation in data.observations {
      let prediction = fit.survival_at(time, observation.covariates)
      let outcome = if observation.event == 1 && observation.time <= time {
        0.0
      } else {
        1.0
      }
      let error = prediction - outcome
      total += error * error
      count += 1
    }
  }
  total / count.to_double()
}

///|
/// Export coefficient estimates and model metadata.
pub fn AftFit::to_csv(self : AftFit) -> String {
  let builder = StringBuilder::new()
  builder.write_string("field,value\n")
  builder.write_string("distribution,")
  builder.write_string(aft_distribution_name(self.distribution))
  builder.write_char('\n')
  builder.write_string("location,")
  builder.write_string(self.location.to_string())
  builder.write_char('\n')
  builder.write_string("scale,")
  builder.write_string(self.scale.to_string())
  builder.write_char('\n')
  builder.write_string("log_likelihood,")
  builder.write_string(self.log_likelihood.to_string())
  builder.write_char('\n')
  for index in 0..<self.coefficients.length() {
    builder.write_string("beta_")
    builder.write_string(index.to_string())
    builder.write_char(',')
    builder.write_string(self.coefficients[index].to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
