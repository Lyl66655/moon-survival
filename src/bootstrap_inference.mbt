///|
/// A bootstrap estimate with percentile and normal intervals.
pub(all) struct SurvivalBootstrapEstimate {
  estimate : Double
  standard_error : Double
  percentile_lower : Double
  percentile_upper : Double
  normal_lower : Double
  normal_upper : Double
  replicates : Int
  bootstrap_row : Bool
} derive(Eq, Debug)

///|
/// One bootstrap survival curve summary.
pub(all) struct SurvivalBootstrapCurvePoint {
  time : Double
  estimate : Double
  lower : Double
  upper : Double
  replicate_count : Int
} derive(Eq, Debug)

///|
/// Deterministically resample core survival data.
pub fn survival_bootstrap_sample(
  data : SurvivalData,
  seed : UInt,
) -> SurvivalData {
  if data.length() == 0 {
    return SurvivalData::new([], [])
  }
  let rng = DeterministicRng::new(seed)
  let times : Array[Double] = []
  let events : Array[Int] = []
  for _ in 0..<data.length() {
    let index = (rng.next_u32() % data.length().reinterpret_as_uint()).reinterpret_as_int()
    times.push(data.get(index).time)
    events.push(data.get(index).event)
  }
  SurvivalData::new(times, events)
}

///|
/// Calculate a sorted percentile from bootstrap values.
pub fn survival_bootstrap_percentile(
  values : Array[Double],
  probability : Double,
) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted : Array[Double] = []
  for value in values {
    sorted.push(value)
  }
  sorted.sort()
  sorted[(clamp_probability(probability) * (sorted.length() - 1).to_double()).to_int()]
}

///|
/// Bootstrap a scalar functional represented by a fitted KM curve at a horizon.
pub fn survival_bootstrap_km(
  data : SurvivalData,
  horizon : Double,
  replicates : Int,
  seed : UInt,
) -> SurvivalBootstrapEstimate {
  if data.length() == 0 || replicates <= 0 {
    return {
      estimate: 0.0,
      standard_error: 0.0,
      percentile_lower: 0.0,
      percentile_upper: 0.0,
      normal_lower: 0.0,
      normal_upper: 0.0,
      replicates: 0,
      bootstrap_row: true,
    }
  }
  let estimate = kaplan_meier(data).survival_at(horizon)
  let values : Array[Double] = []
  for replicate in 0..<replicates {
    let sample = survival_bootstrap_sample(
      data,
      seed + replicate.reinterpret_as_uint(),
    )
    values.push(kaplan_meier(sample).survival_at(horizon))
  }
  let mean = mean_doubles(values)
  let mut variance = 0.0
  for value in values {
    let delta = value - mean
    variance += delta * delta
  }
  let standard_error = double_sqrt(variance / values.length().to_double())
  {
    estimate,
    standard_error,
    percentile_lower: survival_bootstrap_percentile(values, 0.025),
    percentile_upper: survival_bootstrap_percentile(values, 0.975),
    normal_lower: clamp_probability(estimate - 1.96 * standard_error),
    normal_upper: clamp_probability(estimate + 1.96 * standard_error),
    replicates: values.length(),
    bootstrap_row: true,
  }
}

///|
/// Bootstrap restricted mean survival time.
pub fn survival_bootstrap_rmst(
  data : SurvivalData,
  horizon : Double,
  replicates : Int,
  seed : UInt,
) -> SurvivalBootstrapEstimate {
  if data.length() == 0 || replicates <= 0 {
    return {
      estimate: 0.0,
      standard_error: 0.0,
      percentile_lower: 0.0,
      percentile_upper: 0.0,
      normal_lower: 0.0,
      normal_upper: 0.0,
      replicates: 0,
      bootstrap_row: true,
    }
  }
  let estimate = kaplan_meier(data).restricted_mean(horizon)
  let values : Array[Double] = []
  for replicate in 0..<replicates {
    values.push(
      kaplan_meier(
        survival_bootstrap_sample(data, seed + replicate.reinterpret_as_uint()),
      ).restricted_mean(horizon),
    )
  }
  let mean = mean_doubles(values)
  let mut variance = 0.0
  for value in values {
    let delta = value - mean
    variance += delta * delta
  }
  let standard_error = double_sqrt(variance / values.length().to_double())
  {
    estimate,
    standard_error,
    percentile_lower: survival_bootstrap_percentile(values, 0.025),
    percentile_upper: survival_bootstrap_percentile(values, 0.975),
    normal_lower: estimate - 1.96 * standard_error,
    normal_upper: estimate + 1.96 * standard_error,
    replicates: values.length(),
    bootstrap_row: true,
  }
}

///|
/// Build bootstrap KM bands on a query grid.
pub fn survival_bootstrap_curve(
  data : SurvivalData,
  times : Array[Double],
  replicates : Int,
  seed : UInt,
) -> Array[SurvivalBootstrapCurvePoint] {
  let result : Array[SurvivalBootstrapCurvePoint] = []
  let central = kaplan_meier(data)
  for time in times {
    let values : Array[Double] = []
    for replicate in 0..<replicates {
      values.push(
        kaplan_meier(
          survival_bootstrap_sample(
            data,
            seed + replicate.reinterpret_as_uint(),
          ),
        ).survival_at(time),
      )
    }
    result.push({
      time,
      estimate: central.survival_at(time),
      lower: survival_bootstrap_percentile(values, 0.025),
      upper: survival_bootstrap_percentile(values, 0.975),
      replicate_count: values.length(),
    })
  }
  result
}

///|
/// Calculate leave-one-out KM horizon estimates.
pub fn survival_jackknife_km(
  data : SurvivalData,
  horizon : Double,
) -> Array[Double] {
  let result : Array[Double] = []
  if data.length() <= 1 {
    return result
  }
  for omitted in 0..<data.length() {
    let times : Array[Double] = []
    let events : Array[Int] = []
    for index in 0..<data.length() {
      if index != omitted {
        times.push(data.get(index).time)
        events.push(data.get(index).event)
      }
    }
    result.push(
      kaplan_meier(SurvivalData::new(times, events)).survival_at(horizon),
    )
  }
  result
}

///|
/// Compute a jackknife acceleration estimate.
pub fn survival_jackknife_acceleration(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let mean = mean_doubles(values)
  let mut numerator = 0.0
  let mut denominator = 0.0
  for value in values {
    let delta = mean - value
    numerator += delta * delta * delta
    denominator += delta * delta
  }
  if denominator <= survival_epsilon {
    0.0
  } else {
    numerator / (6.0 * double_exp(1.5 * double_log(denominator)))
  }
}

///|
/// Export bootstrap scalar summaries.
pub fn survival_bootstrap_estimate_to_csv(
  estimate : SurvivalBootstrapEstimate,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("estimate,")
  builder.write_string(estimate.estimate.to_string())
  builder.write_char('\n')
  builder.write_string("standard_error,")
  builder.write_string(estimate.standard_error.to_string())
  builder.write_char('\n')
  builder.write_string("percentile_lower,")
  builder.write_string(estimate.percentile_lower.to_string())
  builder.write_char('\n')
  builder.write_string("percentile_upper,")
  builder.write_string(estimate.percentile_upper.to_string())
  builder.write_char('\n')
  builder.write_string("normal_lower,")
  builder.write_string(estimate.normal_lower.to_string())
  builder.write_char('\n')
  builder.write_string("normal_upper,")
  builder.write_string(estimate.normal_upper.to_string())
  builder.write_char('\n')
  builder.write_string("replicates,")
  builder.write_string(estimate.replicates.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Export bootstrap curve bands.
pub fn survival_bootstrap_curve_to_csv(
  points : Array[SurvivalBootstrapCurvePoint],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,estimate,lower,upper,replicate_count\n")
  for point in points {
    builder.write_string(point.time.to_string())
    builder.write_char(',')
    builder.write_string(point.estimate.to_string())
    builder.write_char(',')
    builder.write_string(point.lower.to_string())
    builder.write_char(',')
    builder.write_string(point.upper.to_string())
    builder.write_char(',')
    builder.write_string(point.replicate_count.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
