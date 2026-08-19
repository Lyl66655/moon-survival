///|
/// A hazard estimate on a query grid.
pub(all) struct SurvivalHazardPoint {
  time : Double
  hazard : Double
  cumulative_hazard : Double
  survival : Double
  hazard_row : Bool
} derive(Eq, Debug)

///|
/// Estimate a discrete hazard curve from core survival data.
pub fn survival_hazard_curve(data : SurvivalData) -> Array[SurvivalHazardPoint] {
  let result : Array[SurvivalHazardPoint] = []
  let fit = kaplan_meier(data)
  let mut cumulative = 0.0
  for step in fit.steps {
    let hazard = if step.n_at_risk == 0 {
      0.0
    } else {
      step.n_events.to_double() / step.n_at_risk.to_double()
    }
    cumulative += -double_log(clamp_positive(1.0 - hazard))
    result.push({
      time: step.time,
      hazard,
      cumulative_hazard: cumulative,
      survival: step.survival,
      hazard_row: true,
    })
  }
  result
}

///|
/// Return a hazard at a time using a step convention.
pub fn survival_hazard_at(
  points : Array[SurvivalHazardPoint],
  time : Double,
) -> Double {
  let mut hazard = 0.0
  for point in points {
    if point.time <= time {
      hazard = point.hazard
    } else {
      break
    }
  }
  hazard
}

///|
/// Return cumulative hazard at a time.
pub fn survival_cumulative_hazard_at(
  points : Array[SurvivalHazardPoint],
  time : Double,
) -> Double {
  let mut hazard = 0.0
  for point in points {
    if point.time <= time {
      hazard = point.cumulative_hazard
    } else {
      break
    }
  }
  hazard
}

///|
/// Smooth a hazard curve with a moving average.
pub fn survival_hazard_smooth(
  points : Array[SurvivalHazardPoint],
  window : Int,
) -> Array[SurvivalHazardPoint] {
  let result : Array[SurvivalHazardPoint] = []
  let values : Array[Double] = []
  for point in points {
    values.push(point.hazard)
  }
  let smoothed = survival_ext_rolling_mean(values, window)
  for i in 0..<points.length() {
    result.push({
      time: points[i].time,
      hazard: smoothed[i],
      cumulative_hazard: points[i].cumulative_hazard,
      survival: points[i].survival,
      hazard_row: true,
    })
  }
  result
}

///|
/// Convert cumulative hazard back to survival.
pub fn survival_from_cumulative_hazard(
  points : Array[CurvePoint],
) -> Array[CurvePoint] {
  survival_curve_map(points, fn(value) { double_exp(-value) })
}

///|
/// Compare hazard curves by integrated absolute difference.
pub fn survival_hazard_distance(
  first : Array[SurvivalHazardPoint],
  second : Array[SurvivalHazardPoint],
  times : Array[Double],
) -> Double {
  let mut total = 0.0
  for time in times {
    total += (survival_hazard_at(first, time) - survival_hazard_at(second, time)).abs()
  }
  if times.length() == 0 {
    0.0
  } else {
    total / times.length().to_double()
  }
}

///|
/// Compute a pointwise hazard ratio curve.
pub fn survival_hazard_ratio_curve(
  first : Array[SurvivalHazardPoint],
  second : Array[SurvivalHazardPoint],
  times : Array[Double],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for time in times {
    result.push({
      time,
      value: survival_ext_safe_ratio(
        survival_hazard_at(first, time),
        survival_hazard_at(second, time),
      ),
    })
  }
  result
}

///|
/// Integrate hazard over a horizon.
pub fn survival_hazard_area(
  points : Array[SurvivalHazardPoint],
  horizon : Double,
) -> Double {
  survival_curve_area(survival_hazard_to_curve(points), horizon)
}

///|
/// Convert hazard points to a generic curve.
pub fn survival_hazard_to_curve(
  points : Array[SurvivalHazardPoint],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push({ time: point.time, value: point.hazard })
  }
  result
}

///|
/// Check hazard monotonicity of cumulative hazard.
pub fn survival_hazard_is_valid(points : Array[SurvivalHazardPoint]) -> Bool {
  let mut previous = 0.0
  for point in points {
    if point.time < 0.0 ||
      point.hazard < 0.0 ||
      point.cumulative_hazard < previous ||
      point.survival < 0.0 ||
      point.survival > 1.0 {
      return false
    }
    previous = point.cumulative_hazard
  }
  true
}

///|
/// Export hazard points.
pub fn survival_hazard_to_csv(points : Array[SurvivalHazardPoint]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,hazard,cumulative_hazard,survival\n")
  for point in points {
    builder.write_string(point.time.to_string())
    builder.write_char(',')
    builder.write_string(point.hazard.to_string())
    builder.write_char(',')
    builder.write_string(point.cumulative_hazard.to_string())
    builder.write_char(',')
    builder.write_string(point.survival.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
