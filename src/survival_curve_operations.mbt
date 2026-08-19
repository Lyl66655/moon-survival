///|
/// Validate an ordered curve.
pub fn survival_curve_is_valid(points : Array[CurvePoint]) -> Bool {
  for index in 0..<points.length() {
    if points[index].time < 0.0 ||
      points[index].value.is_nan() ||
      points[index].value.is_inf() ||
      (index > 0 && points[index].time < points[index - 1].time) {
      return false
    }
  }
  true
}

///|
/// Return the time grid of a curve.
pub fn survival_curve_times(points : Array[CurvePoint]) -> Array[Double] {
  let result : Array[Double] = []
  for point in points {
    result.push(point.time)
  }
  result
}

///|
/// Return the values of a curve.
pub fn survival_curve_values(points : Array[CurvePoint]) -> Array[Double] {
  let result : Array[Double] = []
  for point in points {
    result.push(point.value)
  }
  result
}

///|
/// Copy a curve.
pub fn survival_curve_copy(points : Array[CurvePoint]) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push(point)
  }
  result
}

///|
/// Evaluate a curve with a left-continuous step convention.
pub fn survival_curve_step_at(
  points : Array[CurvePoint],
  time : Double,
) -> Double {
  let mut value = if points.length() == 0 { 0.0 } else { points[0].value }
  for point in points {
    if point.time <= time {
      value = point.value
    } else {
      break
    }
  }
  value
}

///|
/// Evaluate a curve by linear interpolation.
pub fn survival_curve_linear_at(
  points : Array[CurvePoint],
  time : Double,
) -> Double {
  survival_ext_interpolate(
    survival_curve_times(points),
    survival_curve_values(points),
    time,
  )
}

///|
/// Sample a curve on a new grid.
pub fn survival_curve_sample(
  points : Array[CurvePoint],
  times : Array[Double],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for time in times {
    result.push({ time, value: survival_curve_step_at(points, time) })
  }
  result
}

///|
/// Add two curves on a shared grid.
pub fn survival_curve_add(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  let times = survival_curve_times(first)
  for time in times {
    result.push({
      time,
      value: survival_curve_step_at(first, time) +
      survival_curve_step_at(second, time),
    })
  }
  result
}

///|
/// Subtract two curves on a shared grid.
pub fn survival_curve_subtract(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  let times = survival_curve_times(first)
  for time in times {
    result.push({
      time,
      value: survival_curve_step_at(first, time) -
      survival_curve_step_at(second, time),
    })
  }
  result
}

///|
/// Scale a curve.
pub fn survival_curve_scale(
  points : Array[CurvePoint],
  factor : Double,
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push({ time: point.time, value: point.value * factor })
  }
  result
}

///|
/// Apply a value transform to a curve.
pub fn survival_curve_map(
  points : Array[CurvePoint],
  transform : (Double) -> Double,
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push({ time: point.time, value: transform(point.value) })
  }
  result
}

///|
/// Compute the area under a step curve.
pub fn survival_curve_area(
  points : Array[CurvePoint],
  horizon : Double,
) -> Double {
  survival_ext_step_area(
    survival_curve_times(points),
    survival_curve_values(points),
    horizon,
  )
}

///|
/// Return a curve's first threshold crossing.
pub fn survival_curve_crossing(
  points : Array[CurvePoint],
  threshold : Double,
) -> Double? {
  for point in points {
    if point.value <= threshold {
      return Some(point.time)
    }
  }
  None
}

///|
/// Return the maximum absolute difference between curves.
pub fn survival_curve_distance(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
  times : Array[Double],
) -> Double {
  let mut distance = 0.0
  for time in times {
    let difference = (survival_curve_step_at(first, time) -
    survival_curve_step_at(second, time)).abs()
    if difference > distance {
      distance = difference
    }
  }
  distance
}

///|
/// Compute an integrated squared curve error.
pub fn survival_curve_squared_error(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
  times : Array[Double],
) -> Double {
  let result : Array[Double] = []
  for time in times {
    let difference = survival_curve_step_at(first, time) -
      survival_curve_step_at(second, time)
    result.push(difference * difference)
  }
  survival_ext_mean(result)
}

///|
/// Return a curve restricted to a time interval.
pub fn survival_curve_window(
  points : Array[CurvePoint],
  start : Double,
  end : Double,
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    if point.time >= start && point.time <= end {
      result.push(point)
    }
  }
  result
}

///|
/// Return a curve with values clamped to a probability interval.
pub fn survival_curve_probability_clamp(
  points : Array[CurvePoint],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in points {
    result.push({ time: point.time, value: clamp_probability(point.value) })
  }
  result
}

///|
/// Make a curve monotone non-increasing by cumulative minimum.
pub fn survival_curve_monotone_decreasing(
  points : Array[CurvePoint],
) -> Array[CurvePoint] {
  let values = survival_ext_cumulative_min(survival_curve_values(points))
  let result : Array[CurvePoint] = []
  for index in 0..<points.length() {
    result.push({ time: points[index].time, value: values[index] })
  }
  result
}

///|
/// Make a curve monotone non-decreasing by cumulative maximum.
pub fn survival_curve_monotone_increasing(
  points : Array[CurvePoint],
) -> Array[CurvePoint] {
  let values = survival_ext_cumulative_max(survival_curve_values(points))
  let result : Array[CurvePoint] = []
  for index in 0..<points.length() {
    result.push({ time: points[index].time, value: values[index] })
  }
  result
}

///|
/// Compute an RMST contrast from two survival curves.
pub fn survival_curve_rmst_contrast(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
  horizon : Double,
) -> Double {
  survival_curve_area(first, horizon) - survival_curve_area(second, horizon)
}

///|
/// Return the pointwise ratio of two curves.
pub fn survival_curve_ratio(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
  times : Array[Double],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for time in times {
    result.push({
      time,
      value: survival_ext_safe_ratio(
        survival_curve_step_at(first, time),
        survival_curve_step_at(second, time),
      ),
    })
  }
  result
}

///|
/// Return the pointwise log ratio of two curves.
pub fn survival_curve_log_ratio(
  first : Array[CurvePoint],
  second : Array[CurvePoint],
  times : Array[Double],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for time in times {
    result.push({
      time,
      value: double_log(
        survival_ext_safe_ratio(
          survival_curve_step_at(first, time),
          survival_curve_step_at(second, time),
        ),
      ),
    })
  }
  result
}

///|
/// Return curve values at their own grid as a matrix.
pub fn survival_curve_matrix(
  curves : Array[Array[CurvePoint]],
  times : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for curve in curves {
    result.push(survival_curve_values(survival_curve_sample(curve, times)))
  }
  result
}

///|
/// Compute a pointwise mean curve.
pub fn survival_curve_mean(
  curves : Array[Array[CurvePoint]],
  times : Array[Double],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for time in times {
    let values : Array[Double] = []
    for curve in curves {
      values.push(survival_curve_step_at(curve, time))
    }
    result.push({ time, value: survival_ext_mean(values) })
  }
  result
}

///|
/// Compute pointwise standard deviation.
pub fn survival_curve_standard_deviation(
  curves : Array[Array[CurvePoint]],
  times : Array[Double],
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for time in times {
    let values : Array[Double] = []
    for curve in curves {
      values.push(survival_curve_step_at(curve, time))
    }
    result.push({ time, value: survival_ext_standard_deviation(values) })
  }
  result
}

///|
/// Export curve operations.
pub fn survival_curve_operation_to_csv(points : Array[CurvePoint]) -> String {
  curve_points_to_csv(points)
}
