///|
/// A curve sample used to interchange estimator output with plotting tools.
pub(all) struct CurvePoint {
  time : Double
  value : Double
} derive(Eq, Debug)

///|
/// Sample a KM curve on a regular grid.
pub fn sample_km_curve(
  fit : KMResult,
  start : Double,
  end : Double,
  points : Int,
) -> Array[CurvePoint] {
  let result = []
  for time in linspace(start, end, points) {
    result.push({ time, value: fit.predict(time) })
  }
  result
}

///|
/// Sample a Nelson-Aalen cumulative hazard curve.
pub fn sample_na_curve(
  fit : NAResult,
  start : Double,
  end : Double,
  points : Int,
) -> Array[CurvePoint] {
  let result = []
  for time in linspace(start, end, points) {
    result.push({ time, value: fit.predict(time) })
  }
  result
}

///|
/// Convert cumulative hazard values into survival values.
pub fn hazard_to_survival(values : Array[Double]) -> Array[Double] {
  let result = []
  for value in values {
    result.push(double_exp(-clamp_positive(value, floor=0.0)))
  }
  result
}

///|
/// Convert survival values into cumulative hazards.
pub fn survival_to_hazard(values : Array[Double]) -> Array[Double] {
  let result = []
  for value in values {
    result.push(
      -double_log(
        clamp_positive(
          1.0 - (1.0 - clamp_probability(value)),
          floor=survival_epsilon,
        ),
      ),
    )
  }
  result
}

///|
/// Interpolate linearly between two curve points.
pub fn interpolate_curve(
  left : CurvePoint,
  right : CurvePoint,
  time : Double,
) -> Double {
  if right.time <= left.time {
    left.value
  } else if time <= left.time {
    left.value
  } else if time >= right.time {
    right.value
  } else {
    let ratio = (time - left.time) / (right.time - left.time)
    left.value + ratio * (right.value - left.value)
  }
}

///|
/// Integrate sampled points with the trapezoid rule.
pub fn integrate_curve(points : Array[CurvePoint]) -> Double {
  if points.length() < 2 {
    0.0
  } else {
    let mut area = 0.0
    for i in 1..<points.length() {
      area += (points[i].time - points[i - 1].time) *
        (points[i].value + points[i - 1].value) /
        2.0
    }
    area
  }
}

///|
/// Resample a cumulative hazard curve into discrete interval hazards.
pub fn interval_hazards(points : Array[CurvePoint]) -> Array[CurvePoint] {
  if points.length() < 2 {
    []
  } else {
    let result = []
    for i in 1..<points.length() {
      result.push({
        time: points[i].time,
        value: points[i].value - points[i - 1].value,
      })
    }
    result
  }
}

///|
/// Locate a time at which a monotone curve reaches a target value.
pub fn inverse_curve(points : Array[CurvePoint], target : Double) -> Double? {
  if points.length() == 0 {
    None
  } else {
    for point in points {
      if point.value >= target {
        return Some(point.time)
      }
    }
    None
  }
}

///|
/// Serialize sampled points as CSV.
pub fn curve_points_to_csv(points : Array[CurvePoint]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,value\n")
  for point in points {
    builder.write_string("\{point.time},\{point.value}\n")
  }
  builder.to_string()
}
