///|
/// Dataset quality report used before fitting any estimator.
pub(all) struct DataQualityReport {
  sample_size : Int
  event_count : Int
  censored_count : Int
  duplicate_time_count : Int
  zero_time_count : Int
  maximum_time : Double
  event_fraction : Double
  warnings : Array[String]
  valid : Bool
} derive(Eq, Debug)

///|
/// Assess ties, event balance, and empty/degenerate follow-up before analysis.
pub fn assess_data_quality(data : SurvivalData) -> DataQualityReport {
  let warnings : Array[String] = []
  let duplicates = event_table(data)
  let mut duplicate_count = 0
  for row in duplicates {
    if row.n_events + row.n_censored > 1 {
      duplicate_count += row.n_events + row.n_censored - 1
    }
  }
  let mut zero_count = 0
  for observation in data.to_array() {
    if observation.time == 0.0 {
      zero_count += 1
    }
  }
  let summary = data.summary()
  if data.length() == 0 {
    warnings.push("empty_dataset")
  }
  if data.length() > 0 && data.event_count() == 0 {
    warnings.push("no_events")
  }
  if data.length() > 0 && data.event_count() == data.length() {
    warnings.push("no_censoring")
  }
  if duplicate_count > 0 {
    warnings.push("tied_times_present")
  }
  if zero_count > 0 {
    warnings.push("zero_time_present")
  }
  {
    sample_size: data.length(),
    event_count: data.event_count(),
    censored_count: data.censored_count(),
    duplicate_time_count: duplicate_count,
    zero_time_count: zero_count,
    maximum_time: summary.max_time,
    event_fraction: summary.event_fraction,
    warnings,
    valid: data.length() > 0 && data.event_count() > 0,
  }
}

///|
/// Validate a predicted-survival vector for a requested horizon.
pub fn validate_predictions(predictions : Array[Double]) -> Bool {
  for prediction in predictions {
    if prediction.is_nan() ||
      prediction.is_inf() ||
      prediction < 0.0 ||
      prediction > 1.0 {
      return false
    }
  }
  true
}

///|
/// A common time grid for reporting and integration.
pub fn analysis_time_grid(data : SurvivalData, points : Int) -> Array[Double] {
  let summary = data.summary()
  if points <= 0 || data.length() == 0 {
    []
  } else {
    linspace(summary.min_time, summary.max_time, points)
  }
}

///|
/// Maximum absolute difference between two step curves at query times.
pub fn curve_max_difference(
  first : KMResult,
  second : KMResult,
  times : Array[Double],
) -> Double {
  let mut maximum = 0.0
  for time in times {
    let difference = (first.predict(time) - second.predict(time)).abs()
    if difference > maximum {
      maximum = difference
    }
  }
  maximum
}

///|
/// Integrated absolute difference between two curves.
pub fn curve_integrated_difference(
  first : KMResult,
  second : KMResult,
  horizon : Double,
  points : Int,
) -> Double {
  if horizon <= 0.0 || points < 2 {
    0.0
  } else {
    let grid = linspace(0.0, horizon, points)
    let mut total = 0.0
    for i in 1..<grid.length() {
      let left = (first.predict(grid[i - 1]) - second.predict(grid[i - 1])).abs()
      let right = (first.predict(grid[i]) - second.predict(grid[i])).abs()
      total += (right + left) * (grid[i] - grid[i - 1]) / 2.0
    }
    total
  }
}

///|
/// Summarize a vector of quality flags as one line.
pub fn quality_flags_to_string(flags : Array[String]) -> String {
  let builder = StringBuilder::new()
  for i in 0..<flags.length() {
    if i > 0 {
      builder.write_string(",")
    }
    builder.write_string(flags[i])
  }
  builder.to_string()
}

///|
/// Return whether a curve is monotone non-increasing on a query grid.
pub fn curve_is_monotone(
  curve : KMResult,
  times : Array[Double],
  tolerance? : Double = 0.000000001,
) -> Bool {
  for i in 1..<times.length() {
    if curve.predict(times[i]) > curve.predict(times[i - 1]) + tolerance {
      return false
    }
  }
  true
}
