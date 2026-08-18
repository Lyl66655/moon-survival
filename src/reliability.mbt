///|
/// Exponential reliability fit under right censoring.
pub(all) struct ExponentialFit {
  rate : Double
  total_time : Double
  events : Int
  log_likelihood : Double
} derive(Eq, Debug)

///|
/// Fit a constant-hazard exponential model using all observed exposure time.
pub fn exponential_fit(data : SurvivalData) -> ExponentialFit {
  let mut total_time = 0.0
  for observation in data.to_array() {
    total_time += observation.time
  }
  let events = data.event_count()
  let rate = if events == 0 || total_time <= survival_epsilon {
    0.0
  } else {
    events.to_double() / total_time
  }
  {
    rate,
    total_time,
    events,
    log_likelihood: if rate <= 0.0 {
      0.0
    } else {
      events.to_double() * double_log(rate) - rate * total_time
    },
  }
}

///|
/// Exponential survival at a time.
pub fn ExponentialFit::survival_at(
  self : ExponentialFit,
  time : Double,
) -> Double {
  if self.rate <= 0.0 {
    1.0
  } else {
    double_exp(-self.rate * time)
  }
}

///|
/// Exponential cumulative hazard at a time.
pub fn ExponentialFit::cumulative_hazard_at(
  self : ExponentialFit,
  time : Double,
) -> Double {
  self.rate * time
}

///|
/// Exponential mean time to failure.
pub fn ExponentialFit::mean_time(self : ExponentialFit) -> Double {
  if self.rate <= 0.0 {
    0.0
  } else {
    1.0 / self.rate
  }
}

///|
/// Exponential quantile for a failure probability.
pub fn ExponentialFit::quantile(
  self : ExponentialFit,
  probability : Double,
) -> Double {
  if self.rate <= 0.0 {
    0.0
  } else {
    -double_log(1.0 - clamp_probability(probability)) / self.rate
  }
}

///|
/// Weibull reliability model parameters.
pub(all) struct WeibullFit {
  shape : Double
  scale : Double
  log_likelihood : Double
  iterations : Int
  converged : Bool
  events : Int
  total_time : Double
} derive(Eq, Debug)

///|
/// Fit a two-parameter Weibull model by profiling the censored likelihood over
/// the shape parameter and solving its one-dimensional score with Newton
/// steps.  This is a useful, dependency-free reliability baseline.
pub fn weibull_fit(data : SurvivalData) -> WeibullFit {
  let events = data.event_count()
  let mut total_time = 0.0
  for observation in data.to_array() {
    total_time += observation.time
  }
  if data.length() == 0 || events == 0 || total_time <= survival_epsilon {
    return {
      shape: 1.0,
      scale: 0.0,
      log_likelihood: 0.0,
      iterations: 0,
      converged: false,
      events,
      total_time,
    }
  }
  let mut shape = 1.0
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<40 {
    let step = 0.0001
    let left = weibull_profile_log_likelihood(data, shape - step)
    let center = weibull_profile_log_likelihood(data, shape)
    let right = weibull_profile_log_likelihood(data, shape + step)
    let first = (right - left) / (2.0 * step)
    let second = (right - 2.0 * center + left) / (step * step)
    if second.abs() <= survival_epsilon {
      break
    }
    let update = first / second
    let candidate = shape - update
    if candidate <= 0.05 || candidate > 20.0 {
      shape = if first > 0.0 { shape * 1.2 } else { shape * 0.8 }
    } else {
      shape = candidate
    }
    if update.abs() <= 0.000001 {
      converged = true
    }
    iterations = iteration + 1
    if converged {
      break
    }
  }
  let scale = weibull_profile_scale(data, shape)
  {
    shape,
    scale,
    log_likelihood: weibull_log_likelihood(data, shape, scale),
    iterations,
    converged,
    events,
    total_time,
  }
}

///|
/// Weibull survival function.
pub fn WeibullFit::survival_at(self : WeibullFit, time : Double) -> Double {
  if self.scale <= 0.0 || time <= 0.0 {
    1.0
  } else {
    double_exp(-@math.pow(time / self.scale, self.shape))
  }
}

///|
/// Weibull cumulative hazard.
pub fn WeibullFit::cumulative_hazard_at(
  self : WeibullFit,
  time : Double,
) -> Double {
  if self.scale <= 0.0 || time <= 0.0 {
    0.0
  } else {
    @math.pow(time / self.scale, self.shape)
  }
}

///|
/// Weibull hazard function.
pub fn WeibullFit::hazard_at(self : WeibullFit, time : Double) -> Double {
  if self.scale <= 0.0 || time <= 0.0 {
    0.0
  } else {
    self.shape / self.scale * @math.pow(time / self.scale, self.shape - 1.0)
  }
}

///|
/// Weibull failure-time quantile.
pub fn WeibullFit::quantile(self : WeibullFit, probability : Double) -> Double {
  if self.scale <= 0.0 {
    0.0
  } else {
    let p = clamp_probability(probability)
    self.scale * @math.pow(-double_log(1.0 - p), 1.0 / self.shape)
  }
}

///|
/// Weibull mean time to failure.
pub fn WeibullFit::mean_time(self : WeibullFit) -> Double {
  if self.scale <= 0.0 {
    0.0
  } else {
    self.scale * @math.exp(@math.ln(1.0 + 1.0 / self.shape))
  }
}

///|
/// Build a reliability report that compares empirical and parametric values.
pub(all) struct ReliabilityReport {
  sample_size : Int
  events : Int
  censoring_fraction : Double
  empirical_rmst : Double
  exponential_rate : Double
  weibull_shape : Double
  weibull_scale : Double
  weibull_log_likelihood : Double
} derive(Eq, Debug)

///|
/// Produce a compact report for maintenance or service-life dashboards.
pub fn reliability_report(
  data : SurvivalData,
  horizon : Double,
) -> ReliabilityReport {
  let exponential = exponential_fit(data)
  let weibull = weibull_fit(data)
  let censoring = data.censored_count()
  {
    sample_size: data.length(),
    events: data.event_count(),
    censoring_fraction: if data.length() == 0 {
      0.0
    } else {
      censoring.to_double() / data.length().to_double()
    },
    empirical_rmst: kaplan_meier(data).restricted_mean(horizon),
    exponential_rate: exponential.rate,
    weibull_shape: weibull.shape,
    weibull_scale: weibull.scale,
    weibull_log_likelihood: weibull.log_likelihood,
  }
}

///|
/// Fit the Weibull scale at a fixed shape under censoring.
fn weibull_profile_scale(data : SurvivalData, shape : Double) -> Double {
  let mut total = 0.0
  for observation in data.to_array() {
    total += @math.pow(observation.time, shape)
  }
  let events = data.event_count()
  if events == 0 {
    0.0
  } else {
    @math.pow(total / events.to_double(), 1.0 / shape)
  }
}

///|
/// Profiled Weibull log likelihood.
fn weibull_profile_log_likelihood(
  data : SurvivalData,
  shape : Double,
) -> Double {
  if shape <= 0.0 {
    -1.0e300
  } else {
    let scale = weibull_profile_scale(data, shape)
    weibull_log_likelihood(data, shape, scale)
  }
}

///|
/// Censored Weibull log likelihood.
fn weibull_log_likelihood(
  data : SurvivalData,
  shape : Double,
  scale : Double,
) -> Double {
  if scale <= 0.0 {
    -1.0e300
  } else {
    let mut result = 0.0
    let log_scale = double_log(scale)
    for observation in data.to_array() {
      let normalized = observation.time / scale
      let power = @math.pow(normalized, shape)
      result += -power
      if observation.event == 1 {
        result += double_log(shape) -
          log_scale +
          (shape - 1.0) * double_log(normalized)
      }
    }
    result
  }
}
