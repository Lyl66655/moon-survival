///|
/// One observation augmented with an inverse-probability-of-censoring weight.
pub(all) struct IpcwObservation {
  time : Double
  event : Int
  weight : Double
  ipcw_row : Bool
} derive(Eq, Debug)

///|
/// A validated weighted survival sample.
pub(all) struct IpcwData {
  observations : Array[IpcwObservation]
  truncation : Double
  effective_sample_size : Double
} derive(Eq, Debug)

///|
/// A weighted survival step.
pub(all) struct IpcwStep {
  time : Double
  weighted_at_risk : Double
  weighted_events : Double
  survival : Double
  standard_error : Double
  lower_ci : Double
  upper_ci : Double
} derive(Eq, Debug)

///|
/// A treatment-group IPCW contrast at a prespecified horizon.
pub(all) struct IpcwTreatmentEffect {
  horizon : Double
  treated_survival : Double
  control_survival : Double
  survival_difference : Double
  treated_risk : Double
  control_risk : Double
  risk_difference : Double
  relative_risk : Double
  treated_n : Int
  control_n : Int
} derive(Eq, Debug)

///|
/// A weighted calibration bin for a predicted survival probability.
pub(all) struct IpcwCalibrationBin {
  lower : Double
  upper : Double
  count : Int
  weighted_count : Double
  mean_prediction : Double
  observed_survival : Double
  absolute_error : Double
} derive(Eq, Debug)

///|
/// Build the censoring distribution G(t) from the observed sample.
fn ipcw_censoring_km(data : SurvivalData) -> KMResult {
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in data.to_array() {
    times.push(observation.time)
    events.push(if observation.event == 0 { 1 } else { 0 })
  }
  kaplan_meier(SurvivalData::new(times, events))
}

///|
/// Estimate a left-limit censoring survival probability at an observation.
fn ipcw_censoring_probability(censoring : KMResult, time : Double) -> Double {
  clamp_positive(censoring.survival_at(time - survival_epsilon))
}

///|
/// Compute stabilized IPCW weights for a right-censored dataset.
pub fn ipcw_weights(
  data : SurvivalData,
  truncation? : Double = 0.05,
) -> IpcwData {
  let lower_bound = if truncation <= 0.0 {
    0.000001
  } else if truncation >= 1.0 {
    0.999999
  } else {
    truncation
  }
  let censoring = ipcw_censoring_km(data)
  let observations : Array[IpcwObservation] = []
  let mut sum_weight = 0.0
  let mut sum_weight_squared = 0.0
  for observation in data.to_array() {
    let probability = ipcw_censoring_probability(censoring, observation.time)
    let bounded = if probability < lower_bound {
      lower_bound
    } else {
      probability
    }
    let weight = 1.0 / bounded
    observations.push({
      time: observation.time,
      event: observation.event,
      weight,
      ipcw_row: true,
    })
    sum_weight += weight
    sum_weight_squared += weight * weight
  }
  let effective = if sum_weight_squared <= survival_epsilon {
    0.0
  } else {
    sum_weight * sum_weight / sum_weight_squared
  }
  { observations, truncation: lower_bound, effective_sample_size: effective }
}

///|
/// Number of weighted observations.
pub fn IpcwData::length(self : IpcwData) -> Int {
  self.observations.length()
}

///|
/// Return an observation by sorted position.
pub fn IpcwData::get(self : IpcwData, index : Int) -> IpcwObservation {
  self.observations[index]
}

///|
/// Return the Kish effective sample size.
pub fn IpcwData::effective_n(self : IpcwData) -> Double {
  self.effective_sample_size
}

///|
/// Copy all weighted rows.
pub fn IpcwData::to_array(self : IpcwData) -> Array[IpcwObservation] {
  let result : Array[IpcwObservation] = []
  for observation in self.observations {
    result.push(observation)
  }
  result
}

///|
/// Estimate a weighted Kaplan-Meier curve with Greenwood-like uncertainty.
pub fn ipcw_kaplan_meier(weighted : IpcwData) -> Array[IpcwStep] {
  let result : Array[IpcwStep] = []
  if weighted.length() == 0 {
    return result
  }
  let mut survival = 1.0
  let mut variance = 0.0
  let mut index = 0
  while index < weighted.length() {
    let time = weighted.get(index).time
    let mut at_risk = 0.0
    let mut events = 0.0
    while index < weighted.length() && weighted.get(index).time == time {
      index += 1
    }
    for observation in weighted.observations {
      if observation.time >= time {
        at_risk += observation.weight
      }
      if observation.time == time && observation.event == 1 {
        events += observation.weight
      }
    }
    if at_risk > survival_epsilon {
      let hazard = clamp_probability(events / at_risk)
      survival *= 1.0 - hazard
      if at_risk > events + survival_epsilon {
        variance += events / (at_risk * (at_risk - events))
      }
    }
    let standard_error = survival * double_sqrt(variance)
    result.push({
      time,
      weighted_at_risk: at_risk,
      weighted_events: events,
      survival,
      standard_error,
      lower_ci: clamp_probability(survival - 1.96 * standard_error),
      upper_ci: clamp_probability(survival + 1.96 * standard_error),
    })
  }
  result
}

///|
/// Predict a weighted survival step function at a horizon.
pub fn ipcw_survival_at(steps : Array[IpcwStep], horizon : Double) -> Double {
  let mut survival = 1.0
  for step in steps {
    if step.time <= horizon {
      survival = step.survival
    } else {
      break
    }
  }
  survival
}

///|
/// Compute a weighted risk at a horizon.
pub fn ipcw_risk_at(weighted : IpcwData, horizon : Double) -> Double {
  if weighted.length() == 0 {
    return 0.0
  }
  let mut numerator = 0.0
  let mut denominator = 0.0
  for observation in weighted.observations {
    let is_known = observation.event == 1 && observation.time <= horizon
    let is_at_risk = observation.time >= horizon || is_known
    if is_at_risk {
      denominator += observation.weight
      if is_known {
        numerator += observation.weight
      }
    }
  }
  if denominator <= survival_epsilon {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Compute a weighted restricted mean survival estimate from steps.
pub fn ipcw_restricted_mean(
  steps : Array[IpcwStep],
  horizon : Double,
) -> Double {
  if horizon <= 0.0 {
    return 0.0
  }
  let mut previous = 0.0
  let mut survival = 1.0
  let mut area = 0.0
  for step in steps {
    if step.time >= horizon {
      break
    }
    area += (step.time - previous) * survival
    previous = step.time
    survival = step.survival
  }
  area + (horizon - previous) * survival
}

///|
/// Estimate a treatment contrast using a binary group indicator.
pub fn ipcw_treatment_effect(
  data : SurvivalData,
  groups : Array[Int],
  horizon : Double,
  truncation? : Double = 0.05,
) -> IpcwTreatmentEffect {
  if groups.length() != data.length() || horizon < 0.0 {
    return {
      horizon,
      treated_survival: 0.0,
      control_survival: 0.0,
      survival_difference: 0.0,
      treated_risk: 0.0,
      control_risk: 0.0,
      risk_difference: 0.0,
      relative_risk: 0.0,
      treated_n: 0,
      control_n: 0,
    }
  }
  let all_rows = data.to_array()
  let treated_times : Array[Double] = []
  let treated_events : Array[Int] = []
  let control_times : Array[Double] = []
  let control_events : Array[Int] = []
  for i in 0..<all_rows.length() {
    if groups[i] == 1 {
      treated_times.push(all_rows[i].time)
      treated_events.push(all_rows[i].event)
    } else {
      control_times.push(all_rows[i].time)
      control_events.push(all_rows[i].event)
    }
  }
  let treated = ipcw_weights(
    SurvivalData::new(treated_times, treated_events),
    truncation~,
  )
  let control = ipcw_weights(
    SurvivalData::new(control_times, control_events),
    truncation~,
  )
  let treated_steps = ipcw_kaplan_meier(treated)
  let control_steps = ipcw_kaplan_meier(control)
  let treated_survival = ipcw_survival_at(treated_steps, horizon)
  let control_survival = ipcw_survival_at(control_steps, horizon)
  let treated_risk = 1.0 - treated_survival
  let control_risk = 1.0 - control_survival
  {
    horizon,
    treated_survival,
    control_survival,
    survival_difference: treated_survival - control_survival,
    treated_risk,
    control_risk,
    risk_difference: treated_risk - control_risk,
    relative_risk: if control_risk <= survival_epsilon {
      0.0
    } else {
      treated_risk / control_risk
    },
    treated_n: treated.length(),
    control_n: control.length(),
  }
}

///|
/// Calculate an IPCW Brier score for binary survival status at a horizon.
pub fn ipcw_brier_score(
  data : SurvivalData,
  predictions : Array[Double],
  horizon : Double,
  truncation? : Double = 0.05,
) -> Double {
  if predictions.length() != data.length() || data.length() == 0 {
    return 0.0
  }
  let weighted = ipcw_weights(data, truncation~)
  let mut total = 0.0
  let mut weight_total = 0.0
  for i in 0..<data.length() {
    let observation = data.get(i)
    let outcome = if observation.event == 1 && observation.time <= horizon {
      0.0
    } else {
      1.0
    }
    let error = predictions[i] - outcome
    total += weighted.get(i).weight * error * error
    weight_total += weighted.get(i).weight
  }
  if weight_total <= survival_epsilon {
    0.0
  } else {
    total / weight_total
  }
}

///|
/// Build equally spaced calibration bins from predicted survival.
pub fn ipcw_calibration(
  data : SurvivalData,
  predictions : Array[Double],
  horizon : Double,
  bins? : Int = 10,
  truncation? : Double = 0.05,
) -> Array[IpcwCalibrationBin] {
  if predictions.length() != data.length() || data.length() == 0 {
    return []
  }
  let bin_count = if bins <= 0 { 1 } else { bins }
  let weighted = ipcw_weights(data, truncation~)
  let result : Array[IpcwCalibrationBin] = []
  for bin in 0..<bin_count {
    let lower = bin.to_double() / bin_count.to_double()
    let upper = (bin + 1).to_double() / bin_count.to_double()
    let mut count = 0
    let mut weight_total = 0.0
    let mut prediction_total = 0.0
    let mut observed_total = 0.0
    for i in 0..<data.length() {
      let prediction = clamp_probability(predictions[i])
      let in_bin = if bin == bin_count - 1 {
        prediction >= lower && prediction <= upper
      } else {
        prediction >= lower && prediction < upper
      }
      if in_bin {
        let observation = data.get(i)
        let outcome = if observation.event == 1 && observation.time <= horizon {
          0.0
        } else {
          1.0
        }
        let weight = weighted.get(i).weight
        count += 1
        weight_total += weight
        prediction_total += weight * prediction
        observed_total += weight * outcome
      }
    }
    let mean_prediction = if weight_total <= survival_epsilon {
      0.0
    } else {
      prediction_total / weight_total
    }
    let observed_survival = if weight_total <= survival_epsilon {
      0.0
    } else {
      observed_total / weight_total
    }
    result.push({
      lower,
      upper,
      count,
      weighted_count: weight_total,
      mean_prediction,
      observed_survival,
      absolute_error: (mean_prediction - observed_survival).abs(),
    })
  }
  result
}

///|
/// Compute integrated calibration error over non-empty bins.
pub fn ipcw_calibration_error(bins : Array[IpcwCalibrationBin]) -> Double {
  let mut numerator = 0.0
  let mut denominator = 0.0
  for bin in bins {
    if bin.count > 0 {
      numerator += bin.weighted_count * bin.absolute_error
      denominator += bin.weighted_count
    }
  }
  if denominator <= survival_epsilon {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Compute the weighted event count and total person-time.
pub fn ipcw_event_rate(weighted : IpcwData, horizon : Double) -> Double {
  if horizon <= 0.0 {
    return 0.0
  }
  let mut event_count = 0.0
  let mut exposure = 0.0
  for observation in weighted.observations {
    let time = if observation.time < horizon {
      observation.time
    } else {
      horizon
    }
    exposure += observation.weight * time
    if observation.event == 1 && observation.time <= horizon {
      event_count += observation.weight
    }
  }
  if exposure <= survival_epsilon {
    0.0
  } else {
    event_count / exposure
  }
}

///|
/// Serialize weighted rows for audit and reproducibility.
pub fn IpcwData::to_csv(self : IpcwData) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,event,ipcw_weight\n")
  for observation in self.observations {
    builder.write_string(observation.time.to_string())
    builder.write_char(',')
    builder.write_string(observation.event.to_string())
    builder.write_char(',')
    builder.write_string(observation.weight.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Serialize a weighted survival curve.
pub fn ipcw_steps_to_csv(steps : Array[IpcwStep]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,weighted_at_risk,weighted_events,survival,standard_error,lower_ci,upper_ci\n",
  )
  for step in steps {
    builder.write_string(step.time.to_string())
    builder.write_char(',')
    builder.write_string(step.weighted_at_risk.to_string())
    builder.write_char(',')
    builder.write_string(step.weighted_events.to_string())
    builder.write_char(',')
    builder.write_string(step.survival.to_string())
    builder.write_char(',')
    builder.write_string(step.standard_error.to_string())
    builder.write_char(',')
    builder.write_string(step.lower_ci.to_string())
    builder.write_char(',')
    builder.write_string(step.upper_ci.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate a weighted sample before exporting it.
pub fn ipcw_is_valid(data : IpcwData) -> Bool {
  if data.truncation <= 0.0 || data.truncation >= 1.0 {
    return false
  }
  for observation in data.observations {
    if observation.time < 0.0 ||
      observation.weight <= 0.0 ||
      observation.weight.is_nan() ||
      observation.weight.is_inf() {
      return false
    }
    if observation.event != 0 && observation.event != 1 {
      return false
    }
  }
  true
}
