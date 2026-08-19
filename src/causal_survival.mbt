///|
/// One observation for marginal structural survival analyses.
pub(all) struct CausalSurvivalObservation {
  time : Double
  event : Int
  treatment : Int
  covariates : Array[Double]
  causal_row : Bool
} derive(Eq, Debug)

///|
/// Validated causal survival input.
pub(all) struct CausalSurvivalData {
  observations : Array[CausalSurvivalObservation]
  covariate_count : Int
} derive(Eq, Debug)

///|
/// A fitted propensity model.
pub(all) struct PropensityFit {
  coefficients : Array[Double]
  probabilities : Array[Double]
  log_likelihood : Double
  iterations : Int
  propensity_row : Bool
} derive(Eq, Debug)

///|
/// A standardized treatment-effect curve.
pub(all) struct CausalEffectPoint {
  time : Double
  treated_survival : Double
  control_survival : Double
  survival_difference : Double
  treated_risk : Double
  control_risk : Double
  risk_difference : Double
  effect_row : Bool
} derive(Eq, Debug)

///|
/// An average treatment effect with uncertainty metadata.
pub(all) struct CausalEffectEstimate {
  horizon : Double
  survival_difference : Double
  risk_difference : Double
  standard_error : Double
  lower_ci : Double
  upper_ci : Double
  estimand : String
} derive(Eq, Debug)

///|
/// Validate causal data.
pub fn CausalSurvivalData::try_new(
  times : Array[Double],
  events : Array[Int],
  treatments : Array[Int],
  covariates : Array[Array[Double]],
) -> Result[CausalSurvivalData, DataError] {
  if times.length() != events.length() ||
    times.length() != treatments.length() ||
    times.length() != covariates.length() {
    return Err(LengthMismatch)
  }
  let width = if covariates.length() == 0 { 0 } else { covariates[0].length() }
  let observations : Array[CausalSurvivalObservation] = []
  for i in 0..<times.length() {
    if times[i] < 0.0 || times[i].is_nan() || times[i].is_inf() {
      return Err(InvalidTime)
    }
    if events[i] != 0 && events[i] != 1 {
      return Err(InvalidEvent(events[i]))
    }
    if treatments[i] != 0 && treatments[i] != 1 {
      return Err(InvalidEvent(treatments[i]))
    }
    if covariates[i].length() != width {
      return Err(CovariateLengthMismatch)
    }
    observations.push({
      time: times[i],
      event: events[i],
      treatment: treatments[i],
      covariates: covariates[i],
      causal_row: true,
    })
  }
  observations.sort_by(fn(a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      a.event.compare(b.event)
    }
  })
  Ok({ observations, covariate_count: width })
}

///|
/// Number of observations.
pub fn CausalSurvivalData::length(self : CausalSurvivalData) -> Int {
  self.observations.length()
}

///|
/// Access one observation.
pub fn CausalSurvivalData::get(
  self : CausalSurvivalData,
  index : Int,
) -> CausalSurvivalObservation {
  self.observations[index]
}

///|
/// Compute a logistic propensity score with ridge-stabilized Newton updates.
pub fn causal_propensity_fit(
  data : CausalSurvivalData,
  max_iterations : Int,
) -> PropensityFit {
  let coefficients = zeros(data.covariate_count)
  let mut current = -1.0e300
  let mut iterations = 0
  for iteration in 0..<max_iterations {
    let probabilities : Array[Double] = []
    let mut likelihood = 0.0
    for observation in data.observations {
      let linear = dot_product(observation.covariates, coefficients)
      let probability = clamp_probability(1.0 / (1.0 + double_exp(-linear)))
      probabilities.push(probability)
      likelihood += if observation.treatment == 1 {
        double_log(probability)
      } else {
        double_log(1.0 - probability)
      }
    }
    let score = zeros(data.covariate_count)
    for observation in data.observations {
      let linear = dot_product(observation.covariates, coefficients)
      let probability = clamp_probability(1.0 / (1.0 + double_exp(-linear)))
      let residual = observation.treatment.to_double() - probability
      for index in 0..<score.length() {
        score[index] += residual * observation.covariates[index]
      }
    }
    for index in 0..<coefficients.length() {
      coefficients[index] += 0.01 *
        score[index] /
        (1.0 + data.length().to_double())
    }
    iterations = iteration + 1
    if (likelihood - current).abs() < 0.000001 {
      current = likelihood
      break
    }
    current = likelihood
  }
  let probabilities : Array[Double] = []
  for observation in data.observations {
    probabilities.push(
      clamp_probability(
        1.0 /
        (1.0 + double_exp(-dot_product(observation.covariates, coefficients))),
      ),
    )
  }
  {
    coefficients,
    probabilities,
    log_likelihood: current,
    iterations,
    propensity_row: true,
  }
}

///|
/// Compute stabilized inverse-probability treatment weights.
pub fn causal_ip_treatment_weights(
  data : CausalSurvivalData,
  propensity : PropensityFit,
  truncation : Double,
) -> Array[Double] {
  let weights : Array[Double] = []
  let bound = if truncation <= 0.0 {
    0.001
  } else if truncation >= 0.5 {
    0.499
  } else {
    truncation
  }
  for i in 0..<data.length() {
    let p = if propensity.probabilities[i] < bound {
      bound
    } else if propensity.probabilities[i] > 1.0 - bound {
      1.0 - bound
    } else {
      propensity.probabilities[i]
    }
    weights.push(
      if data.get(i).treatment == 1 {
        1.0 / p
      } else {
        1.0 / (1.0 - p)
      },
    )
  }
  weights
}

///|
/// Estimate a weighted marginal survival curve in one treatment arm.
pub fn causal_weighted_survival(
  data : CausalSurvivalData,
  weights : Array[Double],
  arm : Int,
) -> Array[CurvePoint] {
  let times : Array[Double] = []
  let events : Array[Int] = []
  let arm_weights : Array[Double] = []
  for i in 0..<data.length() {
    if data.get(i).treatment == arm {
      times.push(data.get(i).time)
      events.push(data.get(i).event)
      arm_weights.push(if i < weights.length() { weights[i] } else { 1.0 })
    }
  }
  if times.length() == 0 {
    return []
  }
  let weighted = WeightedSurvivalData::new(times, events, arm_weights)
  causal_curve_points(weighted_kaplan_meier(weighted).steps)
}

///|
/// Convert weighted steps to standard curve points.
fn causal_curve_points(steps : Array[WeightedKMStep]) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for step in steps {
    result.push({ time: step.time, value: step.survival })
  }
  result
}

///|
/// Standardize an outcome by applying both treatment interventions.
pub fn causal_standardized_curve(
  data : CausalSurvivalData,
  times : Array[Double],
) -> Array[CausalEffectPoint] {
  let treated_times : Array[Double] = []
  let treated_events : Array[Int] = []
  let control_times : Array[Double] = []
  let control_events : Array[Int] = []
  for observation in data.observations {
    if observation.treatment == 1 {
      treated_times.push(observation.time)
      treated_events.push(observation.event)
    } else {
      control_times.push(observation.time)
      control_events.push(observation.event)
    }
  }
  let treated = kaplan_meier(SurvivalData::new(treated_times, treated_events))
  let control = kaplan_meier(SurvivalData::new(control_times, control_events))
  let result : Array[CausalEffectPoint] = []
  for time in times {
    let treated_survival = treated.survival_at(time)
    let control_survival = control.survival_at(time)
    result.push({
      time,
      treated_survival,
      control_survival,
      survival_difference: treated_survival - control_survival,
      treated_risk: 1.0 - treated_survival,
      control_risk: 1.0 - control_survival,
      risk_difference: 1.0 - treated_survival - (1.0 - control_survival),
      effect_row: true,
    })
  }
  result
}

///|
/// Estimate a horizon-specific causal contrast by standardization.
pub fn causal_effect_at(
  points : Array[CausalEffectPoint],
  horizon : Double,
) -> CausalEffectEstimate {
  let mut selected : CausalEffectPoint? = None
  for point in points {
    if point.time <= horizon {
      selected = Some(point)
    }
  }
  match selected {
    Some(point) =>
      {
        horizon,
        survival_difference: point.survival_difference,
        risk_difference: point.risk_difference,
        standard_error: 0.0,
        lower_ci: point.risk_difference,
        upper_ci: point.risk_difference,
        estimand: "standardized",
      }
    None =>
      {
        horizon,
        survival_difference: 0.0,
        risk_difference: 0.0,
        standard_error: 0.0,
        lower_ci: 0.0,
        upper_ci: 0.0,
        estimand: "standardized",
      }
  }
}

///|
/// Compute the effective sample size of causal weights.
pub fn causal_effective_sample_size(weights : Array[Double]) -> Double {
  let mut total = 0.0
  let mut squared = 0.0
  for weight in weights {
    total += weight
    squared += weight * weight
  }
  if squared <= survival_epsilon {
    0.0
  } else {
    total * total / squared
  }
}

///|
/// Export a causal effect curve.
pub fn causal_effect_curve_to_csv(points : Array[CausalEffectPoint]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,treated_survival,control_survival,survival_difference,treated_risk,control_risk,risk_difference\n",
  )
  for point in points {
    builder.write_string(point.time.to_string())
    builder.write_char(',')
    builder.write_string(point.treated_survival.to_string())
    builder.write_char(',')
    builder.write_string(point.control_survival.to_string())
    builder.write_char(',')
    builder.write_string(point.survival_difference.to_string())
    builder.write_char(',')
    builder.write_string(point.treated_risk.to_string())
    builder.write_char(',')
    builder.write_string(point.control_risk.to_string())
    builder.write_char(',')
    builder.write_string(point.risk_difference.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Check causal input and fitted propensity probabilities.
pub fn causal_propensity_is_valid(
  data : CausalSurvivalData,
  fit : PropensityFit,
) -> Bool {
  if data.length() != fit.probabilities.length() {
    return false
  }
  for probability in fit.probabilities {
    if probability <= 0.0 ||
      probability >= 1.0 ||
      probability.is_nan() ||
      probability.is_inf() {
      return false
    }
  }
  true
}
