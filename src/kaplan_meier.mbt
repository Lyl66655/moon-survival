///|
/// A step in the Kaplan-Meier survival curve.
pub(all) struct KMStep {
  time : Double
  n_at_risk : Int
  n_events : Int
  n_censored : Int
  survival : Double
  std_err : Double
  cumulative_variance : Double
  lower_ci : Double
  upper_ci : Double
} derive(Eq, Debug)

///|
/// Result of a Kaplan-Meier estimation.
pub(all) struct KMResult {
  steps : Array[KMStep]
  confidence_level : Double
  sample_size : Int
  event_count : Int
} derive(Eq, Debug)

///|
/// Estimate a Kaplan-Meier curve with Greenwood standard errors and ordinary
/// pointwise confidence intervals.
pub fn kaplan_meier(data : SurvivalData) -> KMResult {
  kaplan_meier_with_confidence(data, confidence=0.95)
}

///|
/// Estimate a Kaplan-Meier curve at a caller-selected confidence level.
pub fn kaplan_meier_with_confidence(
  data : SurvivalData,
  confidence? : Double = 0.95,
) -> KMResult {
  let level = if confidence <= 0.0 {
    0.95
  } else if confidence >= 1.0 {
    0.999999
  } else {
    confidence
  }
  let rows = event_table(data)
  let steps = []
  let mut survival = 1.0
  let mut cumulative_variance = 0.0
  let z = normal_quantile(0.5 + level / 2.0)
  for row in rows {
    if row.n_events > 0 {
      let at_risk = row.n_at_risk.to_double()
      let events = row.n_events.to_double()
      survival *= 1.0 - events / at_risk
      if row.n_at_risk > row.n_events {
        cumulative_variance += events /
          (at_risk * (row.n_at_risk - row.n_events).to_double())
      }
    }
    let std_err = survival * double_sqrt(cumulative_variance)
    let margin = z * std_err
    steps.push({
      time: row.time,
      n_at_risk: row.n_at_risk,
      n_events: row.n_events,
      n_censored: row.n_censored,
      survival,
      std_err,
      cumulative_variance,
      lower_ci: clamp_probability(survival - margin),
      upper_ci: clamp_probability(survival + margin),
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
/// Predict survival at a time using the left-continuous step convention.
pub fn KMResult::predict(self : KMResult, t : Double) -> Double {
  let mut survival = 1.0
  for step in self.steps {
    if step.time <= t {
      survival = step.survival
    } else {
      break
    }
  }
  survival
}

///|
/// Alias that reads naturally in reporting code.
pub fn KMResult::survival_at(self : KMResult, t : Double) -> Double {
  self.predict(t)
}

///|
/// Return the step that contains a time, if the curve has any observations.
pub fn KMResult::step_at(self : KMResult, t : Double) -> KMStep? {
  let mut found : KMStep? = None
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
/// Estimate the discrete hazard at a curve step.
pub fn KMResult::hazard_at(self : KMResult, t : Double) -> Double {
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
/// Return a pointwise confidence interval at the requested time.
pub fn KMResult::confidence_interval(
  self : KMResult,
  t : Double,
) -> (Double, Double) {
  match self.step_at(t) {
    Some(step) => (step.lower_ci, step.upper_ci)
    None => (1.0, 1.0)
  }
}

///|
/// Return the first time at which the estimated survival is at most `p`.
pub fn KMResult::quantile(self : KMResult, p : Double) -> Double? {
  let target = 1.0 - clamp_probability(p)
  for step in self.steps {
    if step.survival <= target {
      return Some(step.time)
    }
  }
  None
}

///|
/// The median survival time, when the curve reaches one half.
pub fn KMResult::median(self : KMResult) -> Double? {
  self.quantile(0.5)
}

///|
/// Integrate the estimated survival curve from zero to a time horizon.
/// This is the restricted mean survival time (RMST).
pub fn KMResult::restricted_mean(self : KMResult, horizon : Double) -> Double {
  if horizon <= 0.0 {
    return 0.0
  }
  let mut area = 0.0
  let mut previous_time = 0.0
  let mut previous_survival = 1.0
  for step in self.steps {
    if step.time >= horizon {
      break
    }
    area += (step.time - previous_time) * previous_survival
    previous_time = step.time
    previous_survival = step.survival
  }
  area + (horizon - previous_time) * previous_survival
}

///|
/// Estimate the area after a landmark time until the supplied horizon.
pub fn KMResult::restricted_mean_after(
  self : KMResult,
  landmark : Double,
  horizon : Double,
) -> Double {
  if horizon <= landmark {
    0.0
  } else {
    self.restricted_mean(horizon) - self.restricted_mean(landmark)
  }
}

///|
/// Number of steps that contain observed events.
pub fn KMResult::observed_event_steps(self : KMResult) -> Int {
  let mut count = 0
  for step in self.steps {
    if step.n_events > 0 {
      count += 1
    }
  }
  count
}

///|
/// Survival probability at a vector of query times.
pub fn KMResult::predict_many(
  self : KMResult,
  times : Array[Double],
) -> Array[Double] {
  let result = []
  for time in times {
    result.push(self.predict(time))
  }
  result
}

///|
/// Export a compact, analysis-friendly event table as CSV.
pub fn KMResult::to_csv(self : KMResult) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,n_at_risk,n_events,n_censored,survival,std_err,lower_ci,upper_ci\n",
  )
  for step in self.steps {
    builder.write_string(
      "\{step.time},\{step.n_at_risk},\{step.n_events},\{step.n_censored},\{step.survival},\{step.std_err},\{step.lower_ci},\{step.upper_ci}\n",
    )
  }
  builder.to_string()
}

///|
/// Weighted Kaplan-Meier step for survey samples or inverse-probability
/// weighted observational studies.
pub(all) struct WeightedKMStep {
  time : Double
  weight_at_risk : Double
  weight_events : Double
  weight_censored : Double
  survival : Double
  std_err : Double
} derive(Eq, Debug)

///|
/// Weighted Kaplan-Meier result.
pub(all) struct WeightedKMResult {
  steps : Array[WeightedKMStep]
} derive(Eq, Debug)

///|
/// Estimate a product-limit curve using positive record weights.
pub fn weighted_kaplan_meier(data : WeightedSurvivalData) -> WeightedKMResult {
  let steps = []
  let mut index = 0
  let mut at_risk = 0.0
  for i in 0..<data.length() {
    at_risk += data.get(i).weight
  }
  let mut survival = 1.0
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
      survival *= 1.0 - event_weight / at_risk
      if at_risk > event_weight {
        variance += event_weight / (at_risk * (at_risk - event_weight))
      }
    }
    steps.push({
      time,
      weight_at_risk: at_risk,
      weight_events: event_weight,
      weight_censored: censor_weight,
      survival,
      std_err: survival * double_sqrt(variance),
    })
    at_risk -= event_weight + censor_weight
  }
  { steps, }
}

///|
/// Predict from a weighted curve.
pub fn WeightedKMResult::predict(self : WeightedKMResult, t : Double) -> Double {
  let mut value = 1.0
  for step in self.steps {
    if step.time <= t {
      value = step.survival
    } else {
      break
    }
  }
  value
}
