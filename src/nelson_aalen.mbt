///|
/// A step in the Nelson-Aalen cumulative hazard curve.
pub(all) struct NAStep {
  time : Double
  n_at_risk : Int
  n_events : Int
  n_censored : Int
  cum_hazard : Double
  variance : Double
  std_err : Double
  lower_ci : Double
  upper_ci : Double
} derive(Eq, Debug)

///|
/// Result of a Nelson-Aalen estimation.
pub(all) struct NAResult {
  steps : Array[NAStep]
  confidence_level : Double
  sample_size : Int
  event_count : Int
} derive(Eq, Debug)

///|
/// Estimates the Nelson-Aalen cumulative hazard curve from the data.
pub fn nelson_aalen(data : SurvivalData) -> NAResult {
  nelson_aalen_with_confidence(data, confidence=0.95)
}

///|
/// Estimate cumulative hazard with Greenwood-like variance and pointwise
/// normal confidence intervals.
pub fn nelson_aalen_with_confidence(
  data : SurvivalData,
  confidence? : Double = 0.95,
) -> NAResult {
  let level = if confidence <= 0.0 {
    0.95
  } else if confidence >= 1.0 {
    0.999999
  } else {
    confidence
  }
  let rows = event_table(data)
  let steps = []
  let mut hazard = 0.0
  let mut variance = 0.0
  let z = normal_quantile(0.5 + level / 2.0)
  for row in rows {
    if row.n_events > 0 {
      let risk = row.n_at_risk.to_double()
      let events = row.n_events.to_double()
      hazard += events / risk
      variance += events / (risk * risk)
    }
    let std_err = double_sqrt(variance)
    steps.push({
      time: row.time,
      n_at_risk: row.n_at_risk,
      n_events: row.n_events,
      n_censored: row.n_censored,
      cum_hazard: hazard,
      variance,
      std_err,
      lower_ci: if hazard - z * std_err < 0.0 {
        0.0
      } else {
        hazard - z * std_err
      },
      upper_ci: hazard + z * std_err,
    })
  }
  {
    steps,
    confidence_level: level,
    sample_size: data.length(),
    event_count: data.event_count(),
  }
}

///|
/// Predict cumulative hazard at a given time.
pub fn NAResult::predict(self : NAResult, t : Double) -> Double {
  let mut hazard = 0.0
  for step in self.steps {
    if step.time <= t {
      hazard = step.cum_hazard
    } else {
      break
    }
  }
  hazard
}

///|
/// Alias for reporting code.
pub fn NAResult::cumulative_hazard_at(self : NAResult, t : Double) -> Double {
  self.predict(t)
}

///|
/// Return the estimated instantaneous Nelson-Aalen increment at a time.
pub fn NAResult::increment_at(self : NAResult, t : Double) -> Double {
  match self.step_at(t) {
    Some(step) =>
      if step.n_at_risk == 0 {
        0.0
      } else {
        step.n_events.to_double() / step.n_at_risk.to_double()
      }
    None => 0.0
  }
}

///|
/// Find the latest hazard step before a time.
pub fn NAResult::step_at(self : NAResult, t : Double) -> NAStep? {
  let mut found : NAStep? = None
  for step in self.steps {
    if step.time <= t {
      found = Some(step)
    } else {
      break
    }
  }
  found
}

///|
/// Return confidence limits at a time.
pub fn NAResult::confidence_interval(
  self : NAResult,
  t : Double,
) -> (Double, Double) {
  match self.step_at(t) {
    Some(step) => (step.lower_ci, step.upper_ci)
    None => (0.0, 0.0)
  }
}

///|
/// Convert a cumulative hazard estimate to a survival estimate.
pub fn NAResult::survival_at(self : NAResult, t : Double) -> Double {
  double_exp(-self.predict(t))
}

///|
/// Return all cumulative hazard values at the requested times.
pub fn NAResult::predict_many(
  self : NAResult,
  times : Array[Double],
) -> Array[Double] {
  let result = []
  for time in times {
    result.push(self.predict(time))
  }
  result
}

///|
/// Export the cumulative hazard table as CSV.
pub fn NAResult::to_csv(self : NAResult) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,n_at_risk,n_events,n_censored,cum_hazard,variance,std_err,lower_ci,upper_ci\n",
  )
  for step in self.steps {
    builder.write_string(
      "\{step.time},\{step.n_at_risk},\{step.n_events},\{step.n_censored},\{step.cum_hazard},\{step.variance},\{step.std_err},\{step.lower_ci},\{step.upper_ci}\n",
    )
  }
  builder.to_string()
}

///|
/// A Nelson-Aalen estimate on weighted observations.
pub(all) struct WeightedNAStep {
  time : Double
  weight_at_risk : Double
  weight_events : Double
  weight_censored : Double
  cum_hazard : Double
  variance : Double
} derive(Eq, Debug)

///|
/// Weighted cumulative hazard result.
pub(all) struct WeightedNAResult {
  steps : Array[WeightedNAStep]
} derive(Eq, Debug)

///|
/// Estimate a weighted cumulative hazard curve.
pub fn weighted_nelson_aalen(data : WeightedSurvivalData) -> WeightedNAResult {
  let steps = []
  let mut at_risk = 0.0
  for i in 0..<data.length() {
    at_risk += data.get(i).weight
  }
  let mut index = 0
  let mut hazard = 0.0
  let mut variance = 0.0
  while index < data.length() {
    let time = data.get(index).time
    let mut event_weight = 0.0
    let mut censor_weight = 0.0
    while index < data.length() && data.get(index).time == time {
      let observation = data.get(index)
      if observation.event == 1 {
        event_weight += observation.weight
      } else {
        censor_weight += observation.weight
      }
      index += 1
    }
    if event_weight > 0.0 {
      hazard += event_weight / at_risk
      variance += event_weight / (at_risk * at_risk)
    }
    steps.push({
      time,
      weight_at_risk: at_risk,
      weight_events: event_weight,
      weight_censored: censor_weight,
      cum_hazard: hazard,
      variance,
    })
    at_risk -= event_weight + censor_weight
  }
  { steps, }
}

///|
/// Predict from a weighted cumulative hazard result.
pub fn WeightedNAResult::predict(self : WeightedNAResult, t : Double) -> Double {
  let mut hazard = 0.0
  for step in self.steps {
    if step.time <= t {
      hazard = step.cum_hazard
    } else {
      break
    }
  }
  hazard
}
