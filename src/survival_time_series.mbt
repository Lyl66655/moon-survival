///|
/// A time-indexed survival monitoring observation.
pub(all) struct SurvivalTimeSeriesPoint {
  time : Double
  value : Double
  exposure : Double
  event_count : Int
  series_point : Bool
} derive(Eq, Debug)

///|
/// A time-series trend summary.
pub(all) struct SurvivalTrendSummary {
  slope : Double
  intercept : Double
  correlation : Double
  first_value : Double
  last_value : Double
  trend_row : Bool
} derive(Eq, Debug)

///|
/// Aggregate event data into fixed time bins.
pub fn survival_time_series_aggregate(
  data : SurvivalData,
  width : Double,
  bins : Int,
) -> Array[SurvivalTimeSeriesPoint] {
  let result : Array[SurvivalTimeSeriesPoint] = []
  if width <= 0.0 || bins <= 0 {
    return result
  }
  for bin in 0..<bins {
    let start = bin.to_double() * width
    let end = start + width
    let mut exposure = 0.0
    let mut events = 0
    let mut values = 0
    for observation in data.to_array() {
      if observation.time >= start && observation.time < end {
        values += 1
        exposure += end - observation.time
        if observation.event == 1 {
          events += 1
        }
      }
    }
    result.push({
      time: start + width / 2.0,
      value: if values == 0 {
        0.0
      } else {
        events.to_double() / values.to_double()
      },
      exposure,
      event_count: events,
      series_point: true,
    })
  }
  result
}

///|
/// Estimate a linear trend by least squares.
pub fn survival_time_series_trend(
  points : Array[SurvivalTimeSeriesPoint],
) -> SurvivalTrendSummary {
  if points.length() < 2 {
    return {
      slope: 0.0,
      intercept: if points.length() == 0 {
        0.0
      } else {
        points[0].value
      },
      correlation: 0.0,
      first_value: if points.length() == 0 {
        0.0
      } else {
        points[0].value
      },
      last_value: if points.length() == 0 {
        0.0
      } else {
        points[points.length() - 1].value
      },
      trend_row: true,
    }
  }
  let times : Array[Double] = []
  let values : Array[Double] = []
  for point in points {
    times.push(point.time)
    values.push(point.value)
  }
  let time_mean = survival_ext_mean(times)
  let value_mean = survival_ext_mean(values)
  let mut covariance = 0.0
  let mut time_variance = 0.0
  let mut value_variance = 0.0
  for i in 0..<times.length() {
    let dt = times[i] - time_mean
    let dv = values[i] - value_mean
    covariance += dt * dv
    time_variance += dt * dt
    value_variance += dv * dv
  }
  let slope = survival_ext_safe_ratio(covariance, time_variance)
  {
    slope,
    intercept: value_mean - slope * time_mean,
    correlation: survival_ext_safe_ratio(
      covariance,
      double_sqrt(time_variance * value_variance),
    ),
    first_value: values[0],
    last_value: values[values.length() - 1],
    trend_row: true,
  }
}

///|
/// Forecast a value from a trend.
pub fn survival_time_series_forecast(
  trend : SurvivalTrendSummary,
  time : Double,
) -> Double {
  trend.intercept + trend.slope * time
}

///|
/// Detect a trend threshold crossing.
pub fn survival_time_series_crossing(
  points : Array[SurvivalTimeSeriesPoint],
  threshold : Double,
) -> Double? {
  for point in points {
    if point.value >= threshold {
      return Some(point.time)
    }
  }
  None
}

///|
/// Smooth a time series with a moving average.
pub fn survival_time_series_smooth(
  points : Array[SurvivalTimeSeriesPoint],
  window : Int,
) -> Array[SurvivalTimeSeriesPoint] {
  let result : Array[SurvivalTimeSeriesPoint] = []
  let values : Array[Double] = []
  for point in points {
    values.push(point.value)
  }
  let smoothed = survival_ext_rolling_mean(values, window)
  for i in 0..<points.length() {
    result.push({
      time: points[i].time,
      value: smoothed[i],
      exposure: points[i].exposure,
      event_count: points[i].event_count,
      series_point: true,
    })
  }
  result
}

///|
/// Compute a cumulative event count series.
pub fn survival_time_series_cumulative(
  points : Array[SurvivalTimeSeriesPoint],
) -> Array[SurvivalTimeSeriesPoint] {
  let result : Array[SurvivalTimeSeriesPoint] = []
  let mut events = 0
  let mut exposure = 0.0
  for point in points {
    events += point.event_count
    exposure += point.exposure
    result.push({
      time: point.time,
      value: events.to_double(),
      exposure,
      event_count: events,
      series_point: true,
    })
  }
  result
}

///|
/// Normalize event rates by exposure.
pub fn survival_time_series_rate(
  points : Array[SurvivalTimeSeriesPoint],
) -> Array[SurvivalTimeSeriesPoint] {
  let result : Array[SurvivalTimeSeriesPoint] = []
  for point in points {
    result.push({
      time: point.time,
      value: survival_ext_safe_ratio(
        point.event_count.to_double(),
        point.exposure,
      ),
      exposure: point.exposure,
      event_count: point.event_count,
      series_point: true,
    })
  }
  result
}

///|
/// Export time-series points.
pub fn survival_time_series_to_csv(
  points : Array[SurvivalTimeSeriesPoint],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,value,exposure,event_count\n")
  for point in points {
    builder.write_string(point.time.to_string())
    builder.write_char(',')
    builder.write_string(point.value.to_string())
    builder.write_char(',')
    builder.write_string(point.exposure.to_string())
    builder.write_char(',')
    builder.write_string(point.event_count.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate a time series.
pub fn survival_time_series_is_valid(
  points : Array[SurvivalTimeSeriesPoint],
) -> Bool {
  for i in 0..<points.length() {
    if points[i].time < 0.0 ||
      points[i].exposure < 0.0 ||
      points[i].event_count < 0 ||
      (i > 0 && points[i].time < points[i - 1].time) {
      return false
    }
  }
  true
}
