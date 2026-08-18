///|
/// A single observation for the Cox proportional-hazards model.
pub(all) struct CoxObservation {
  time : Double
  event : Int
  covariates : Array[Double]
} derive(Eq, Debug)

///|
/// Validated and time-sorted Cox model input.
pub(all) struct CoxData {
  observations : Array[CoxObservation]
  n_covariates : Int
} derive(Eq, Debug)

///|
/// Create Cox data without panicking on external input.
pub fn CoxData::try_new(
  times : Array[Double],
  events : Array[Int],
  covs : Array[Array[Double]],
) -> Result[CoxData, DataError] {
  if times.length() != events.length() || times.length() != covs.length() {
    return Err(LengthMismatch)
  }
  let n_covariates = if covs.length() == 0 { 0 } else { covs[0].length() }
  let observations = []
  for i in 0..<times.length() {
    if times[i].is_nan() || times[i].is_inf() {
      return Err(InvalidTime)
    }
    if times[i] < 0.0 {
      return Err(NegativeTime)
    }
    if events[i] != 0 && events[i] != 1 {
      return Err(InvalidEvent(events[i]))
    }
    if covs[i].length() != n_covariates {
      return Err(CovariateLengthMismatch)
    }
    for value in covs[i] {
      if value.is_nan() || value.is_inf() {
        return Err(InvalidTime)
      }
    }
    observations.push({ time: times[i], event: events[i], covariates: covs[i] })
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
  Ok({ observations, n_covariates })
}

///|
/// Panicking convenience constructor retained for the original API.
pub fn CoxData::new(
  times : Array[Double],
  events : Array[Int],
  covs : Array[Array[Double]],
) -> CoxData {
  match CoxData::try_new(times, events, covs) {
    Ok(data) => data
    Err(_) => panic()
  }
}

///|
/// Number of observations in a Cox dataset.
pub fn CoxData::length(self : CoxData) -> Int {
  self.observations.length()
}

///|
/// Access a sorted Cox observation.
pub fn CoxData::get(self : CoxData, index : Int) -> CoxObservation {
  self.observations[index]
}

///|
/// Return a detached copy of covariate rows.
pub fn CoxData::covariates(self : CoxData) -> Array[Array[Double]] {
  let result = []
  for observation in self.observations {
    let row = []
    for value in observation.covariates {
      row.push(value)
    }
    result.push(row)
  }
  result
}

///|
/// Compute a linear predictor for one covariate vector.
pub fn cox_linear_predictor(
  covariates : Array[Double],
  beta : Array[Double],
) -> Double {
  dot_product(covariates, beta)
}

///|
/// Compute the Cox partial likelihood using Breslow handling for ties.
pub fn cox_partial_likelihood(data : CoxData, beta : Array[Double]) -> Double {
  if beta.length() != data.n_covariates {
    panic()
  }
  let mut likelihood = 1.0
  let mut index = 0
  while index < data.length() {
    let time = data.get(index).time
    let mut event_count = 0
    let mut numerator = 1.0
    while index < data.length() && data.get(index).time == time {
      let observation = data.get(index)
      if observation.event == 1 {
        event_count += 1
        numerator *= double_exp(
          cox_linear_predictor(observation.covariates, beta),
        )
      }
      index += 1
    }
    if event_count > 0 {
      let denominator = cox_risk_sum(data, beta, time)
      for _ in 0..<event_count {
        likelihood *= numerator / double_exp(double_log(denominator))
      }
    }
  }
  likelihood
}

///|
/// Log partial likelihood, which remains useful when the product underflows.
pub fn cox_log_partial_likelihood(
  data : CoxData,
  beta : Array[Double],
) -> Double {
  if beta.length() != data.n_covariates {
    panic()
  }
  let mut result = 0.0
  let mut index = 0
  while index < data.length() {
    let time = data.get(index).time
    let mut events = 0
    let mut numerator = 0.0
    while index < data.length() && data.get(index).time == time {
      let observation = data.get(index)
      if observation.event == 1 {
        events += 1
        numerator += dot_product(observation.covariates, beta)
      }
      index += 1
    }
    if events > 0 {
      result += numerator -
        events.to_double() * double_log(cox_risk_sum(data, beta, time))
    }
  }
  result
}

///|
/// Score vector of the Breslow partial likelihood.
pub fn cox_score_vector(data : CoxData, beta : Array[Double]) -> Array[Double] {
  let p = data.n_covariates
  if beta.length() != p {
    panic()
  }
  let score = zeros(p)
  let mut index = 0
  while index < data.length() {
    let time = data.get(index).time
    let events = []
    while index < data.length() && data.get(index).time == time {
      let observation = data.get(index)
      if observation.event == 1 {
        events.push(observation.covariates)
      }
      index += 1
    }
    if events.length() > 0 {
      let risk_sum = cox_risk_sum(data, beta, time)
      let risk_weighted = zeros(p)
      for observation in data.observations {
        if observation.time >= time {
          let weight = double_exp(dot_product(observation.covariates, beta))
          for j in 0..<p {
            risk_weighted[j] += weight * observation.covariates[j]
          }
        }
      }
      for event in events {
        for j in 0..<p {
          score[j] += event[j] - risk_weighted[j] / risk_sum
        }
      }
    }
  }
  score
}

///|
/// Observed information matrix (negative Hessian) for the Cox likelihood.
pub fn cox_information_matrix(
  data : CoxData,
  beta : Array[Double],
) -> Array[Array[Double]] {
  let p = data.n_covariates
  if beta.length() != p {
    panic()
  }
  let information = identity_matrix(p)
  for i in 0..<p {
    information[i][i] = 0.0
  }
  let mut index = 0
  while index < data.length() {
    let time = data.get(index).time
    let mut event_count = 0
    while index < data.length() && data.get(index).time == time {
      if data.get(index).event == 1 {
        event_count += 1
      }
      index += 1
    }
    if event_count > 0 {
      let risk_sum = cox_risk_sum(data, beta, time)
      let first = zeros(p)
      let second = []
      for _ in 0..<p {
        second.push(zeros(p))
      }
      for observation in data.observations {
        if observation.time >= time {
          let weight = double_exp(dot_product(observation.covariates, beta))
          for j in 0..<p {
            first[j] += weight * observation.covariates[j]
            for k in 0..<p {
              second[j][k] += weight *
                observation.covariates[j] *
                observation.covariates[k]
            }
          }
        }
      }
      for j in 0..<p {
        for k in 0..<p {
          let covariance = second[j][k] / risk_sum -
            first[j] * first[k] / (risk_sum * risk_sum)
          information[j][k] += event_count.to_double() * covariance
        }
      }
    }
  }
  matrix_symmetrize(information)
}

///|
/// Options controlling Newton-Raphson fitting.
pub(all) struct CoxFitOptions {
  max_iterations : Int
  tolerance : Double
  ridge : Double
  step_halving : Bool
} derive(Eq, Debug)

///|
/// Default fitting options for production-sized datasets.
pub fn cox_default_options() -> CoxFitOptions {
  {
    max_iterations: 50,
    tolerance: 0.000001,
    ridge: 0.000000001,
    step_halving: true,
  }
}

///|
/// One fitted coefficient with Wald diagnostics.
pub(all) struct CoxCoefficient {
  index : Int
  estimate : Double
  standard_error : Double
  z_value : Double
  p_value : Double
  hazard_ratio : Double
  lower_hr : Double
  upper_hr : Double
} derive(Eq, Debug)

///|
/// A step of the fitted baseline cumulative hazard.
pub(all) struct CoxBaselineStep {
  time : Double
  n_at_risk : Int
  n_events : Int
  hazard_increment : Double
  cumulative_hazard : Double
  baseline_survival : Double
} derive(Eq, Debug)

///|
/// Full Cox model fit and diagnostics.
pub(all) struct CoxFitResult {
  coefficients : Array[Double]
  standard_errors : Array[Double]
  covariance : Array[Array[Double]]
  log_likelihood : Double
  null_log_likelihood : Double
  iterations : Int
  converged : Bool
  coefficients_table : Array[CoxCoefficient]
  baseline : Array[CoxBaselineStep]
} derive(Eq, Debug)

///|
/// Fit a Cox model with convenient optional controls.
pub fn cox_fit(
  data : CoxData,
  max_iterations? : Int = 50,
  tolerance? : Double = 0.000001,
  ridge? : Double = 0.000000001,
) -> CoxFitResult {
  let options = { max_iterations, tolerance, ridge, step_halving: true }
  cox_fit_with_options(data, options)
}

///|
/// Fit using an explicit options record.
pub fn cox_fit_with_options(
  data : CoxData,
  options : CoxFitOptions,
) -> CoxFitResult {
  let mut beta = zeros(data.n_covariates)
  let null_likelihood = cox_log_partial_likelihood(data, beta)
  let mut current_likelihood = null_likelihood
  let mut converged = data.n_covariates == 0
  let mut iterations = 0
  for iteration in 0..<options.max_iterations {
    if converged {
      break
    }
    let score = cox_score_vector(data, beta)
    let information = matrix_ridge(
      cox_information_matrix(data, beta),
      options.ridge,
    )
    let step = match solve_linear_system(information, score) {
      Some(value) => value
      None => break
    }
    let mut factor = 1.0
    let candidate = vector_add(beta, step)
    let mut candidate_likelihood = cox_log_partial_likelihood(data, candidate)
    if options.step_halving {
      while candidate_likelihood < current_likelihood && factor > 0.0001 {
        factor *= 0.5
        candidate_likelihood = cox_log_partial_likelihood(
          data,
          vector_add(beta, vector_scale(step, factor)),
        )
      }
    }
    beta = vector_add(beta, vector_scale(step, factor))
    current_likelihood = candidate_likelihood
    iterations = iteration + 1
    if vector_norm(vector_scale(step, factor)) <= options.tolerance {
      converged = true
    }
  }
  let final_information = matrix_ridge(
    cox_information_matrix(data, beta),
    options.ridge,
  )
  let covariance = match matrix_inverse(final_information) {
    Some(value) => value
    None => identity_matrix(data.n_covariates)
  }
  let standard_errors = []
  for i in 0..<data.n_covariates {
    standard_errors.push(double_sqrt(covariance[i][i].abs()))
  }
  let coefficients_table = []
  let z = normal_quantile(0.975)
  for i in 0..<data.n_covariates {
    let standard_error = standard_errors[i]
    let z_value = if standard_error <= survival_epsilon {
      0.0
    } else {
      beta[i] / standard_error
    }
    let p_value = 2.0 * (1.0 - normal_cdf(z_value.abs()))
    let estimate = beta[i]
    coefficients_table.push({
      index: i,
      estimate,
      standard_error,
      z_value,
      p_value,
      hazard_ratio: double_exp(estimate),
      lower_hr: double_exp(estimate - z * standard_error),
      upper_hr: double_exp(estimate + z * standard_error),
    })
  }
  {
    coefficients: beta,
    standard_errors,
    covariance,
    log_likelihood: current_likelihood,
    null_log_likelihood: null_likelihood,
    iterations,
    converged,
    coefficients_table,
    baseline: cox_baseline(data, beta),
  }
}

///|
/// Predict a relative risk from a fitted coefficient vector.
pub fn cox_predict_risk(
  covariates : Array[Double],
  coefficients : Array[Double],
) -> Double {
  double_exp(cox_linear_predictor(covariates, coefficients))
}

///|
/// Predict an individual survival curve from a fitted model.
pub fn cox_predict_survival(
  fit : CoxFitResult,
  covariates : Array[Double],
  time : Double,
) -> Double {
  let linear_predictor = cox_linear_predictor(covariates, fit.coefficients)
  let baseline_hazard = cox_baseline_cumulative_hazard(fit.baseline, time)
  double_exp(-baseline_hazard * double_exp(linear_predictor))
}

///|
/// Return the baseline cumulative hazard at a time.
pub fn cox_baseline_cumulative_hazard(
  baseline : Array[CoxBaselineStep],
  time : Double,
) -> Double {
  let mut hazard = 0.0
  for step in baseline {
    if step.time <= time {
      hazard = step.cumulative_hazard
    } else {
      break
    }
  }
  hazard
}

///|
/// Return the fitted baseline survival at a time.
pub fn CoxFitResult::baseline_survival(
  self : CoxFitResult,
  time : Double,
) -> Double {
  double_exp(-cox_baseline_cumulative_hazard(self.baseline, time))
}

///|
/// Return deviance explained relative to the null partial likelihood.
pub fn CoxFitResult::likelihood_ratio(self : CoxFitResult) -> Double {
  2.0 * (self.log_likelihood - self.null_log_likelihood)
}

///|
/// Return the model's linear predictor.
pub fn CoxFitResult::linear_predictor(
  self : CoxFitResult,
  covariates : Array[Double],
) -> Double {
  cox_linear_predictor(covariates, self.coefficients)
}

///|
/// Return an individual relative risk.
pub fn CoxFitResult::risk(
  self : CoxFitResult,
  covariates : Array[Double],
) -> Double {
  cox_predict_risk(covariates, self.coefficients)
}

///|
/// Return an individual survival probability.
pub fn CoxFitResult::survival_at(
  self : CoxFitResult,
  covariates : Array[Double],
  time : Double,
) -> Double {
  cox_predict_survival(self, covariates, time)
}

///|
/// Compute score residuals for every observation.
pub fn cox_score_residuals(
  data : CoxData,
  beta : Array[Double],
) -> Array[Array[Double]] {
  let result = []
  let global_score = cox_score_vector(data, beta)
  for observation in data.observations {
    let value = vector_scale(
      observation.covariates,
      observation.event.to_double(),
    )
    result.push(
      vector_subtract(
        value,
        vector_scale(global_score, 1.0 / data.length().to_double()),
      ),
    )
  }
  result
}

///|
/// Construct the Breslow baseline cumulative hazard.
priv struct CoxEventRow {
  time : Double
  mut n_events : Int
}

///|
fn cox_baseline(data : CoxData, beta : Array[Double]) -> Array[CoxBaselineStep] {
  let result = []
  let rows : Array[CoxEventRow] = []
  for observation in data.observations {
    if rows.length() == 0 || rows[rows.length() - 1].time != observation.time {
      rows.push({ time: observation.time, n_events: 0 })
    }
    if observation.event == 1 {
      rows[rows.length() - 1].n_events += 1
    }
  }
  let mut cumulative = 0.0
  for row in rows {
    let denominator = cox_risk_sum(data, beta, row.time)
    let increment = row.n_events.to_double() / denominator
    cumulative += increment
    result.push({
      time: row.time,
      n_at_risk: cox_n_at_risk(data, row.time),
      n_events: row.n_events,
      hazard_increment: increment,
      cumulative_hazard: cumulative,
      baseline_survival: double_exp(-cumulative),
    })
  }
  result
}

///|
/// Compute a risk-set denominator using a max-shifted exponential sum.
fn cox_risk_sum(data : CoxData, beta : Array[Double], time : Double) -> Double {
  let mut maximum = -700.0
  for observation in data.observations {
    if observation.time >= time {
      let value = dot_product(observation.covariates, beta)
      if value > maximum {
        maximum = value
      }
    }
  }
  let mut scaled_sum = 0.0
  for observation in data.observations {
    if observation.time >= time {
      scaled_sum += double_exp(
        dot_product(observation.covariates, beta) - maximum,
      )
    }
  }
  double_exp(maximum) * scaled_sum
}

///|
/// Count members of a risk set.
fn cox_n_at_risk(data : CoxData, time : Double) -> Int {
  let mut count = 0
  for observation in data.observations {
    if observation.time >= time {
      count += 1
    }
  }
  count
}
