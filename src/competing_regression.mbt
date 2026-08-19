///|
/// A cause-specific hazard step from a competing-risk sample.
pub(all) struct CauseSpecificStep {
  time : Double
  cause : Int
  at_risk : Int
  events : Int
  hazard_increment : Double
  cumulative_hazard : Double
  cause_step : Bool
} derive(Eq, Debug)

///|
/// A fitted cause-specific cumulative hazard model.
pub(all) struct CauseSpecificFit {
  cause : Int
  steps : Array[CauseSpecificStep]
  sample_size : Int
  event_count : Int
  fit_row : Bool
} derive(Eq, Debug)

///|
/// Predicted cause-specific risk at a horizon.
pub(all) struct CauseSpecificPrediction {
  time : Double
  cause : Int
  cumulative_hazard : Double
  survival : Double
  cumulative_incidence : Double
  prediction_row : Bool
} derive(Eq, Debug)

///|
/// Fit a crude cause-specific cumulative hazard from event histories.
pub fn cause_specific_fit(
  data : CompetingRiskData,
  cause : Int,
) -> CauseSpecificFit {
  let steps : Array[CauseSpecificStep] = []
  if cause <= 0 || cause > data.n_causes {
    return {
      cause,
      steps,
      sample_size: data.length(),
      event_count: 0,
      fit_row: true,
    }
  }
  let times : Array[Double] = []
  for observation in data.observations {
    if observation.cause == cause && !times.contains(observation.time) {
      times.push(observation.time)
    }
  }
  times.sort()
  let mut cumulative = 0.0
  let mut events_total = 0
  for time in times {
    let mut at_risk = 0
    let mut events = 0
    for observation in data.observations {
      if observation.time >= time {
        at_risk += 1
      }
      if observation.time == time && observation.cause == cause {
        events += 1
      }
    }
    let increment = if at_risk == 0 {
      0.0
    } else {
      events.to_double() / at_risk.to_double()
    }
    cumulative += increment
    events_total += events
    steps.push({
      time,
      cause,
      at_risk,
      events,
      hazard_increment: increment,
      cumulative_hazard: cumulative,
      cause_step: true,
    })
  }
  {
    cause,
    steps,
    sample_size: data.length(),
    event_count: events_total,
    fit_row: true,
  }
}

///|
/// Fit all positive causes.
pub fn cause_specific_fit_all(
  data : CompetingRiskData,
) -> Array[CauseSpecificFit] {
  let result : Array[CauseSpecificFit] = []
  for cause in 1..<(data.n_causes + 1) {
    result.push(cause_specific_fit(data, cause))
  }
  result
}

///|
/// Return cumulative cause-specific hazard at a time.
pub fn cause_specific_hazard_at(
  fit : CauseSpecificFit,
  time : Double,
) -> Double {
  let mut value = 0.0
  for step in fit.steps {
    if step.time <= time {
      value = step.cumulative_hazard
    } else {
      break
    }
  }
  value
}

///|
/// Compute all-cause survival from fitted cause-specific hazards.
pub fn competing_survival_from_fits(
  fits : Array[CauseSpecificFit],
  time : Double,
) -> Double {
  let mut cumulative = 0.0
  for fit in fits {
    cumulative += cause_specific_hazard_at(fit, time)
  }
  double_exp(-cumulative)
}

///|
/// Predict a cumulative incidence using a discrete hazard decomposition.
pub fn cause_specific_predict(
  fit : CauseSpecificFit,
  fits : Array[CauseSpecificFit],
  time : Double,
) -> CauseSpecificPrediction {
  let mut incidence = 0.0
  let mut previous = 0.0
  let mut survival = 1.0
  for step in fit.steps {
    if step.time > time {
      break
    }
    let all_survival_before = competing_survival_from_fits(
      fits,
      step.time - survival_epsilon,
    )
    incidence += all_survival_before *
      (1.0 - double_exp(-step.hazard_increment))
    previous = step.time
    survival = competing_survival_from_fits(fits, previous)
  }
  {
    time,
    cause: fit.cause,
    cumulative_hazard: cause_specific_hazard_at(fit, time),
    survival,
    cumulative_incidence: clamp_probability(incidence),
    prediction_row: true,
  }
}

///|
/// Predict all causes on one time grid.
pub fn competing_regression_curve(
  fits : Array[CauseSpecificFit],
  times : Array[Double],
) -> Array[CauseSpecificPrediction] {
  let result : Array[CauseSpecificPrediction] = []
  for time in times {
    for fit in fits {
      result.push(cause_specific_predict(fit, fits, time))
    }
  }
  result
}

///|
/// Compute cause-specific hazard ratios from two fitted hazards.
pub fn cause_specific_rate_ratio(
  first : CauseSpecificFit,
  second : CauseSpecificFit,
  horizon : Double,
) -> Double {
  let first_hazard = cause_specific_hazard_at(first, horizon)
  let second_hazard = cause_specific_hazard_at(second, horizon)
  if second_hazard <= survival_epsilon {
    0.0
  } else {
    first_hazard / second_hazard
  }
}

///|
/// Compute the sum of cause-specific cumulative incidences.
pub fn competing_total_incidence(
  predictions : Array[CauseSpecificPrediction],
  time : Double,
) -> Double {
  let mut total = 0.0
  for prediction in predictions {
    if prediction.time == time {
      total += prediction.cumulative_incidence
    }
  }
  clamp_probability(total)
}

///|
/// Export cause-specific steps.
pub fn cause_specific_steps_to_csv(steps : Array[CauseSpecificStep]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,cause,at_risk,events,hazard_increment,cumulative_hazard\n",
  )
  for step in steps {
    builder.write_string(step.time.to_string())
    builder.write_char(',')
    builder.write_string(step.cause.to_string())
    builder.write_char(',')
    builder.write_string(step.at_risk.to_string())
    builder.write_char(',')
    builder.write_string(step.events.to_string())
    builder.write_char(',')
    builder.write_string(step.hazard_increment.to_string())
    builder.write_char(',')
    builder.write_string(step.cumulative_hazard.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export cumulative-incidence predictions.
pub fn cause_specific_predictions_to_csv(
  predictions : Array[CauseSpecificPrediction],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,cause,cumulative_hazard,survival,cumulative_incidence\n",
  )
  for prediction in predictions {
    builder.write_string(prediction.time.to_string())
    builder.write_char(',')
    builder.write_string(prediction.cause.to_string())
    builder.write_char(',')
    builder.write_string(prediction.cumulative_hazard.to_string())
    builder.write_char(',')
    builder.write_string(prediction.survival.to_string())
    builder.write_char(',')
    builder.write_string(prediction.cumulative_incidence.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate fitted cause-specific hazards.
pub fn cause_specific_fit_is_valid(fit : CauseSpecificFit) -> Bool {
  if fit.cause <= 0 || fit.sample_size < 0 || fit.event_count < 0 {
    return false
  }
  let mut previous = 0.0
  for step in fit.steps {
    if step.time < 0.0 ||
      step.hazard_increment < 0.0 ||
      step.cumulative_hazard < previous {
      return false
    }
    previous = step.cumulative_hazard
  }
  true
}
