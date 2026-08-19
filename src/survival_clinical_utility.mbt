///|
/// A utility-weighted survival point.
pub(all) struct SurvivalUtilityPoint {
  time : Double
  survival : Double
  utility : Double
  qaly_density : Double
  cumulative_qaly : Double
  utility_row : Bool
} derive(Eq, Debug)

///|
/// An expected-value summary over a horizon.
pub(all) struct SurvivalUtilitySummary {
  horizon : Double
  restricted_mean_survival : Double
  restricted_mean_qaly : Double
  discounted_qaly : Double
  average_utility : Double
  summary_row : Bool
} derive(Eq, Debug)

///|
/// Build a utility-weighted curve.
pub fn survival_utility_curve(
  points : Array[CurvePoint],
  utilities : Array[Double],
) -> Array[SurvivalUtilityPoint] {
  let result : Array[SurvivalUtilityPoint] = []
  let mut cumulative = 0.0
  for i in 0..<points.length() {
    let utility = if i < utilities.length() {
      clamp_probability(utilities[i])
    } else {
      1.0
    }
    let density = points[i].value * utility
    if i > 0 {
      cumulative += (points[i].time - points[i - 1].time) * density
    }
    result.push({
      time: points[i].time,
      survival: points[i].value,
      utility,
      qaly_density: density,
      cumulative_qaly: cumulative,
      utility_row: true,
    })
  }
  result
}

///|
/// Integrate utility-weighted survival to a horizon.
pub fn survival_utility_rmst(
  points : Array[SurvivalUtilityPoint],
  horizon : Double,
) -> Double {
  if points.length() == 0 || horizon <= 0.0 {
    return 0.0
  }
  let curve : Array[CurvePoint] = []
  for point in points {
    curve.push({ time: point.time, value: point.qaly_density })
  }
  survival_curve_area(curve, horizon)
}

///|
/// Apply continuous discounting to utility points.
pub fn survival_utility_discount(
  points : Array[SurvivalUtilityPoint],
  rate : Double,
) -> Array[SurvivalUtilityPoint] {
  let result : Array[SurvivalUtilityPoint] = []
  let mut cumulative = 0.0
  for point in points {
    let discounted = point.qaly_density *
      double_exp(-clamp_positive(rate, floor=0.0) * point.time)
    cumulative += discounted
    result.push({
      time: point.time,
      survival: point.survival,
      utility: point.utility,
      qaly_density: discounted,
      cumulative_qaly: cumulative,
      utility_row: true,
    })
  }
  result
}

///|
/// Summarize QALYs over a survival curve.
pub fn survival_utility_summary(
  points : Array[SurvivalUtilityPoint],
  horizon : Double,
  discount_rate : Double,
) -> SurvivalUtilitySummary {
  let discounted = survival_utility_discount(points, discount_rate)
  let rmst = survival_curve_area(survival_utility_to_curve(points), horizon)
  let qaly = survival_curve_area(survival_utility_to_curve(discounted), horizon)
  {
    horizon,
    restricted_mean_survival: survival_curve_area(
      survival_utility_survival_curve(points),
      horizon,
    ),
    restricted_mean_qaly: rmst,
    discounted_qaly: qaly,
    average_utility: if rmst <= survival_epsilon {
      0.0
    } else {
      qaly / rmst
    },
    summary_row: true,
  }
}

///|
/// Convert utility points to a generic QALY-density curve.
pub fn survival_utility_to_curve(
  points : Array[SurvivalUtilityPoint],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push({ time: point.time, value: point.qaly_density })
  }
  result
}

///|
/// Convert utility points to a survival curve.
pub fn survival_utility_survival_curve(
  points : Array[SurvivalUtilityPoint],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push({ time: point.time, value: point.survival })
  }
  result
}

///|
/// Compute QALY difference between two utility curves.
pub fn survival_utility_difference(
  first : Array[SurvivalUtilityPoint],
  second : Array[SurvivalUtilityPoint],
  horizon : Double,
) -> Double {
  survival_utility_rmst(first, horizon) - survival_utility_rmst(second, horizon)
}

///|
/// Create a utility profile from scalar utility.
pub fn survival_scalar_utility_curve(
  points : Array[CurvePoint],
  utility : Double,
) -> Array[SurvivalUtilityPoint] {
  survival_utility_curve(points, Array::make(points.length(), utility))
}

///|
/// Export utility points.
pub fn survival_utility_to_csv(points : Array[SurvivalUtilityPoint]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,survival,utility,qaly_density,cumulative_qaly\n")
  for point in points {
    builder.write_string(point.time.to_string())
    builder.write_char(',')
    builder.write_string(point.survival.to_string())
    builder.write_char(',')
    builder.write_string(point.utility.to_string())
    builder.write_char(',')
    builder.write_string(point.qaly_density.to_string())
    builder.write_char(',')
    builder.write_string(point.cumulative_qaly.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Check utility range and time ordering.
pub fn survival_utility_is_valid(points : Array[SurvivalUtilityPoint]) -> Bool {
  for i in 0..<points.length() {
    if points[i].time < 0.0 ||
      points[i].utility < 0.0 ||
      points[i].utility > 1.0 ||
      (i > 0 && points[i].time < points[i - 1].time) {
      return false
    }
  }
  true
}
