///|
/// A segment of a piecewise-constant hazard model.
pub(all) struct HazardSegment {
  start : Double
  end : Double
  exposure : Double
  events : Int
  rate : Double
  mut cumulative_hazard : Double
} derive(Eq, Debug)

///|
/// Piecewise exponential model fitted from right-censored observations.
pub(all) struct PiecewiseExponentialFit {
  segments : Array[HazardSegment]
  total_events : Int
  total_exposure : Double
  log_likelihood : Double
} derive(Eq, Debug)

///|
/// Fit a piecewise-constant hazard using exposure accounting. Cutpoints must
/// be sorted and the final endpoint must cover the analysis horizon.
pub fn piecewise_exponential_fit(
  data : SurvivalData,
  cutpoints : Array[Double],
) -> PiecewiseExponentialFit {
  if cutpoints.length() == 0 {
    panic()
  }
  for i in 1..<cutpoints.length() {
    if cutpoints[i] <= cutpoints[i - 1] {
      panic()
    }
  }
  let segments = []
  let mut previous = 0.0
  let mut total_exposure = 0.0
  let mut total_events = 0
  let mut log_likelihood = 0.0
  for endpoint in cutpoints {
    let mut exposure = 0.0
    let mut events = 0
    for observation in data.to_array() {
      let observed_end = if observation.time < endpoint {
        observation.time
      } else {
        endpoint
      }
      if observed_end > previous {
        let segment_exposure = observed_end - previous
        if observation.time >= previous {
          exposure += segment_exposure
        }
      }
      if observation.event == 1 &&
        observation.time > previous &&
        observation.time <= endpoint {
        events += 1
      }
    }
    let rate = if exposure <= survival_epsilon {
      0.0
    } else {
      events.to_double() / exposure
    }
    if events > 0 && rate > 0.0 {
      log_likelihood += events.to_double() * double_log(rate) - rate * exposure
    } else {
      log_likelihood += -rate * exposure
    }
    segments.push({
      start: previous,
      end: endpoint,
      exposure,
      events,
      rate,
      cumulative_hazard: 0.0,
    })
    total_exposure += exposure
    total_events += events
    previous = endpoint
  }
  let mut cumulative = 0.0
  for segment in segments {
    cumulative += segment.rate * (segment.end - segment.start)
    segment.cumulative_hazard = cumulative
  }
  { segments, total_events, total_exposure, log_likelihood }
}

///|
/// Cumulative hazard of a piecewise exponential fit.
pub fn PiecewiseExponentialFit::cumulative_hazard_at(
  self : PiecewiseExponentialFit,
  time : Double,
) -> Double {
  let mut cumulative = 0.0
  let mut previous = 0.0
  for segment in self.segments {
    if time <= segment.start {
      break
    }
    let endpoint = if time < segment.end { time } else { segment.end }
    if endpoint > previous {
      cumulative += segment.rate * (endpoint - previous)
    }
    previous = segment.end
    if time <= segment.end {
      break
    }
  }
  cumulative
}

///|
/// Piecewise exponential survival probability.
pub fn PiecewiseExponentialFit::survival_at(
  self : PiecewiseExponentialFit,
  time : Double,
) -> Double {
  double_exp(-self.cumulative_hazard_at(time))
}

///|
/// Hazard rate in the segment containing a time.
pub fn PiecewiseExponentialFit::hazard_at(
  self : PiecewiseExponentialFit,
  time : Double,
) -> Double {
  for segment in self.segments {
    if time >= segment.start && time <= segment.end {
      return segment.rate
    }
  }
  0.0
}

///|
/// Numerically locate a quantile of a piecewise exponential model.
pub fn PiecewiseExponentialFit::quantile(
  self : PiecewiseExponentialFit,
  probability : Double,
) -> Double {
  let target = -double_log(1.0 - clamp_probability(probability))
  let mut previous_hazard = 0.0
  for segment in self.segments {
    let segment_hazard = segment.rate * (segment.end - segment.start)
    if previous_hazard + segment_hazard >= target && segment.rate > 0.0 {
      return segment.start + (target - previous_hazard) / segment.rate
    }
    previous_hazard += segment_hazard
  }
  if self.segments.length() == 0 {
    0.0
  } else {
    self.segments[self.segments.length() - 1].end
  }
}

///|
/// Export a piecewise hazard model as CSV.
pub fn PiecewiseExponentialFit::to_csv(
  self : PiecewiseExponentialFit,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("start,end,exposure,events,rate,cumulative_hazard\n")
  for segment in self.segments {
    builder.write_string(
      "\{segment.start},\{segment.end},\{segment.exposure},\{segment.events},\{segment.rate},\{segment.cumulative_hazard}\n",
    )
  }
  builder.to_string()
}

///|
/// A grouped survival estimate for stratified reporting.
pub(all) struct StratifiedSurvival {
  stratum : String
  sample_size : Int
  events : Int
  survival : Double
  restricted_mean : Double
} derive(Eq, Debug)

///|
/// Fit one Kaplan-Meier curve per named stratum.
pub fn stratified_survival(
  strata : Array[String],
  groups : Array[SurvivalData],
  horizon : Double,
) -> Array[StratifiedSurvival] {
  if strata.length() != groups.length() {
    panic()
  }
  let result = []
  for i in 0..<groups.length() {
    let km = kaplan_meier(groups[i])
    result.push({
      stratum: strata[i],
      sample_size: groups[i].length(),
      events: groups[i].event_count(),
      survival: km.predict(horizon),
      restricted_mean: km.restricted_mean(horizon),
    })
  }
  result
}

///|
/// An actuarial life-table interval.
pub(all) struct LifeTableRow {
  start : Double
  end : Double
  n_entering : Int
  withdrawals : Int
  events : Int
  effective_at_risk : Double
  conditional_survival : Double
  cumulative_survival : Double
} derive(Eq, Debug)

///|
/// Create an actuarial life table with half-interval withdrawal correction.
pub fn actuarial_life_table(
  data : SurvivalData,
  cutpoints : Array[Double],
) -> Array[LifeTableRow] {
  if cutpoints.length() == 0 {
    return []
  }
  let rows = []
  let mut cumulative_survival = 1.0
  let mut entering = data.length()
  let mut previous = 0.0
  for endpoint in cutpoints {
    let mut withdrawals = 0
    let mut events = 0
    for observation in data.to_array() {
      if observation.time > previous && observation.time <= endpoint {
        if observation.event == 1 {
          events += 1
        } else {
          withdrawals += 1
        }
      }
    }
    let effective = entering.to_double() - withdrawals.to_double() / 2.0
    let conditional = if effective <= survival_epsilon {
      1.0
    } else {
      clamp_probability(1.0 - events.to_double() / effective)
    }
    cumulative_survival *= conditional
    rows.push({
      start: previous,
      end: endpoint,
      n_entering: entering,
      withdrawals,
      events,
      effective_at_risk: effective,
      conditional_survival: conditional,
      cumulative_survival,
    })
    entering -= withdrawals + events
    previous = endpoint
  }
  rows
}
