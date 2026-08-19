///|
/// A calibration regression summary for a survival endpoint.
pub(all) struct SurvivalCalibrationFit {
  intercept : Double
  slope : Double
  calibration_in_the_large : Double
  standard_error_intercept : Double
  standard_error_slope : Double
  horizon : Double
  calibration_row : Bool
} derive(Eq, Debug)

///|
/// A smooth calibration-curve point.
pub(all) struct SurvivalCalibrationPoint {
  predicted : Double
  observed : Double
  count : Int
  weight : Double
  lower : Double
  upper : Double
} derive(Eq, Debug)

///|
/// A global calibration score.
pub(all) struct SurvivalCalibrationScore {
  brier : Double
  calibration_error : Double
  slope_error : Double
  intercept_error : Double
  integrated_calibration_index : Double
  horizon : Double
} derive(Eq, Debug)

///|
/// Estimate calibration-in-the-large from weighted binary outcomes.
pub fn survival_calibration_fit(
  data : SurvivalData,
  predictions : Array[Double],
  horizon : Double,
) -> SurvivalCalibrationFit {
  if predictions.length() != data.length() || data.length() == 0 {
    return {
      intercept: 0.0,
      slope: 0.0,
      calibration_in_the_large: 0.0,
      standard_error_intercept: 0.0,
      standard_error_slope: 0.0,
      horizon,
      calibration_row: false,
    }
  }
  let observed : Array[Double] = []
  let logit : Array[Double] = []
  for i in 0..<data.length() {
    let outcome = if data.get(i).event == 1 && data.get(i).time <= horizon {
      1.0
    } else {
      0.0
    }
    observed.push(outcome)
    let p = clamp_probability(predictions[i])
    logit.push(double_log(p / clamp_positive(1.0 - p)))
  }
  let mean_observed = mean_doubles(observed)
  let mean_logit = mean_doubles(logit)
  let mut covariance = 0.0
  let mut variance = 0.0
  for i in 0..<data.length() {
    covariance += (observed[i] - mean_observed) * (logit[i] - mean_logit)
    variance += (logit[i] - mean_logit) * (logit[i] - mean_logit)
  }
  let slope = if variance <= survival_epsilon {
    0.0
  } else {
    covariance / variance
  }
  let intercept = mean_observed - slope * mean_logit
  let mut residual_sum = 0.0
  for value in observed {
    let delta = value - mean_observed
    residual_sum += delta * delta
  }
  let standard_error = double_sqrt(
    residual_sum / (data.length() * data.length()).to_double(),
  )
  {
    intercept,
    slope,
    calibration_in_the_large: mean_observed - mean_doubles(predictions),
    standard_error_intercept: standard_error,
    standard_error_slope: if variance <= survival_epsilon {
      0.0
    } else {
      standard_error / double_sqrt(variance)
    },
    horizon,
    calibration_row: true,
  }
}

///|
/// Construct quantile-based calibration groups.
pub fn survival_calibration_curve(
  data : SurvivalData,
  predictions : Array[Double],
  horizon : Double,
  groups : Int,
) -> Array[SurvivalCalibrationPoint] {
  if predictions.length() != data.length() || data.length() == 0 || groups <= 0 {
    return []
  }
  let indices : Array[Int] = []
  for i in 0..<data.length() {
    indices.push(i)
  }
  indices.sort_by(fn(a, b) {
    if predictions[a] < predictions[b] {
      -1
    } else if predictions[a] > predictions[b] {
      1
    } else {
      a.compare(b)
    }
  })
  let result : Array[SurvivalCalibrationPoint] = []
  for group in 0..<groups {
    let start = group * indices.length() / groups
    let end = (group + 1) * indices.length() / groups
    let mut prediction_total = 0.0
    let mut outcome_total = 0.0
    let mut count = 0
    for position in start..<end {
      if position < indices.length() {
        let index = indices[position]
        prediction_total += predictions[index]
        outcome_total += if data.get(index).event == 1 &&
          data.get(index).time <= horizon {
          1.0
        } else {
          0.0
        }
        count += 1
      }
    }
    let predicted = if count == 0 {
      0.0
    } else {
      prediction_total / count.to_double()
    }
    let observed = if count == 0 {
      0.0
    } else {
      outcome_total / count.to_double()
    }
    let standard_error = if count == 0 {
      0.0
    } else {
      double_sqrt(observed * (1.0 - observed) / count.to_double())
    }
    result.push({
      predicted,
      observed,
      count,
      weight: count.to_double(),
      lower: clamp_probability(observed - 1.96 * standard_error),
      upper: clamp_probability(observed + 1.96 * standard_error),
    })
  }
  result
}

///|
/// Compute the integrated calibration index.
pub fn survival_integrated_calibration_index(
  points : Array[SurvivalCalibrationPoint],
) -> Double {
  let mut total = 0.0
  let mut weight = 0.0
  for point in points {
    total += point.weight * (point.predicted - point.observed).abs()
    weight += point.weight
  }
  if weight <= survival_epsilon {
    0.0
  } else {
    total / weight
  }
}

///|
/// Compute a calibration score bundle.
pub fn survival_calibration_score(
  data : SurvivalData,
  predictions : Array[Double],
  horizon : Double,
  groups : Int,
) -> SurvivalCalibrationScore {
  let fit = survival_calibration_fit(data, predictions, horizon)
  let curve = survival_calibration_curve(data, predictions, horizon, groups)
  let brier = ipcw_brier_score(data, predictions, horizon)
  {
    brier,
    calibration_error: fit.calibration_in_the_large.abs(),
    slope_error: (fit.slope - 1.0).abs(),
    intercept_error: fit.intercept.abs(),
    integrated_calibration_index: survival_integrated_calibration_index(curve),
    horizon,
  }
}

///|
/// Compare calibration curves by integrated absolute difference.
pub fn survival_calibration_curve_distance(
  first : Array[SurvivalCalibrationPoint],
  second : Array[SurvivalCalibrationPoint],
) -> Double {
  let count = if first.length() < second.length() {
    first.length()
  } else {
    second.length()
  }
  if count == 0 {
    return 0.0
  }
  let mut total = 0.0
  for i in 0..<count {
    total += (first[i].observed - second[i].observed).abs()
  }
  total / count.to_double()
}

///|
/// Export calibration points.
pub fn survival_calibration_points_to_csv(
  points : Array[SurvivalCalibrationPoint],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("predicted,observed,count,weight,lower,upper\n")
  for point in points {
    builder.write_string(point.predicted.to_string())
    builder.write_char(',')
    builder.write_string(point.observed.to_string())
    builder.write_char(',')
    builder.write_string(point.count.to_string())
    builder.write_char(',')
    builder.write_string(point.weight.to_string())
    builder.write_char(',')
    builder.write_string(point.lower.to_string())
    builder.write_char(',')
    builder.write_string(point.upper.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export calibration fit diagnostics.
pub fn survival_calibration_fit_to_csv(fit : SurvivalCalibrationFit) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("intercept,")
  builder.write_string(fit.intercept.to_string())
  builder.write_char('\n')
  builder.write_string("slope,")
  builder.write_string(fit.slope.to_string())
  builder.write_char('\n')
  builder.write_string("calibration_in_the_large,")
  builder.write_string(fit.calibration_in_the_large.to_string())
  builder.write_char('\n')
  builder.write_string("standard_error_intercept,")
  builder.write_string(fit.standard_error_intercept.to_string())
  builder.write_char('\n')
  builder.write_string("standard_error_slope,")
  builder.write_string(fit.standard_error_slope.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Validate calibration points.
pub fn survival_calibration_is_valid(
  points : Array[SurvivalCalibrationPoint],
) -> Bool {
  for point in points {
    if point.predicted < 0.0 ||
      point.predicted > 1.0 ||
      point.observed < 0.0 ||
      point.observed > 1.0 ||
      point.count < 0 {
      return false
    }
  }
  true
}
