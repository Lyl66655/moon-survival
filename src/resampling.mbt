///|
/// Percentile interval for a scalar bootstrap estimate.
pub(all) struct PercentileInterval {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  replicates : Int
} derive(Eq, Debug)

///|
/// Deterministic bootstrap for a scalar statistic on a vector.
pub fn bootstrap_scalar(
  values : Array[Double],
  statistic : (Array[Double]) -> Double,
  replicates : Int,
  seed : UInt,
) -> PercentileInterval {
  if values.length() == 0 || replicates <= 0 {
    return {
      estimate: 0.0,
      standard_error: 0.0,
      lower: 0.0,
      upper: 0.0,
      replicates: 0,
    }
  }
  let estimate = statistic(values)
  let rng = DeterministicRng::new(seed)
  let samples = []
  for _ in 0..<replicates {
    let sample = []
    for _ in 0..<values.length() {
      let index = (rng.next_u32() % values.length().reinterpret_as_uint()).reinterpret_as_int()
      sample.push(values[index])
    }
    samples.push(statistic(sample))
  }
  let standard_error = double_sqrt(population_variance(samples))
  {
    estimate,
    standard_error,
    lower: quantile_doubles(samples, 0.025).unwrap_or(estimate),
    upper: quantile_doubles(samples, 0.975).unwrap_or(estimate),
    replicates,
  }
}

///|
/// Bootstrap the sample mean.
pub fn bootstrap_mean(
  values : Array[Double],
  replicates : Int,
  seed : UInt,
) -> PercentileInterval {
  bootstrap_scalar(values, mean_doubles, replicates, seed)
}

///|
/// Jackknife estimate and standard error for a scalar statistic.
pub fn jackknife_scalar(
  values : Array[Double],
  statistic : (Array[Double]) -> Double,
) -> PercentileInterval {
  if values.length() < 2 {
    return {
      estimate: statistic(values),
      standard_error: 0.0,
      lower: 0.0,
      upper: 0.0,
      replicates: 0,
    }
  }
  let estimates = []
  for omitted in 0..<values.length() {
    let sample = []
    for i in 0..<values.length() {
      if i != omitted {
        sample.push(values[i])
      }
    }
    estimates.push(statistic(sample))
  }
  let mean = mean_doubles(estimates)
  let mut total = 0.0
  for estimate in estimates {
    let difference = estimate - mean
    total += difference * difference
  }
  let standard_error = double_sqrt(
    (values.length() - 1).to_double() / values.length().to_double() * total,
  )
  {
    estimate: statistic(values),
    standard_error,
    lower: statistic(values) - 1.96 * standard_error,
    upper: statistic(values) + 1.96 * standard_error,
    replicates: values.length(),
  }
}

///|
/// Bootstrap RMST at a fixed horizon.
pub fn bootstrap_rmst(
  data : SurvivalData,
  horizon : Double,
  replicates : Int,
  seed : UInt,
) -> PercentileInterval {
  let observations = data.to_array()
  if data.length() == 0 {
    {
      estimate: 0.0,
      standard_error: 0.0,
      lower: 0.0,
      upper: 0.0,
      replicates: 0,
    }
  } else {
    let estimate = kaplan_meier(data).restricted_mean(horizon)
    let rng = DeterministicRng::new(seed)
    let samples = []
    for _ in 0..<replicates {
      let sample_times = []
      let sample_events = []
      for _ in 0..<data.length() {
        let index = (rng.next_u32() % data.length().reinterpret_as_uint()).reinterpret_as_int()
        sample_times.push(observations[index].time)
        sample_events.push(observations[index].event)
      }
      samples.push(
        kaplan_meier(SurvivalData::new(sample_times, sample_events)).restricted_mean(
          horizon,
        ),
      )
    }
    {
      estimate,
      standard_error: double_sqrt(population_variance(samples)),
      lower: quantile_doubles(samples, 0.025).unwrap_or(estimate),
      upper: quantile_doubles(samples, 0.975).unwrap_or(estimate),
      replicates,
    }
  }
}

///|
/// Difference in two bootstrap mean estimates.
pub fn bootstrap_mean_difference(
  left : Array[Double],
  right : Array[Double],
  replicates : Int,
  seed : UInt,
) -> PercentileInterval {
  let combined = []
  for value in left {
    combined.push(value)
  }
  for value in right {
    combined.push(value)
  }
  let estimate = mean_doubles(left) - mean_doubles(right)
  let rng = DeterministicRng::new(seed)
  let samples = []
  for _ in 0..<replicates {
    let left_sample = []
    let right_sample = []
    for _ in 0..<left.length() {
      let index = (rng.next_u32() % combined.length().reinterpret_as_uint()).reinterpret_as_int()
      left_sample.push(combined[index])
    }
    for _ in 0..<right.length() {
      let index = (rng.next_u32() % combined.length().reinterpret_as_uint()).reinterpret_as_int()
      right_sample.push(combined[index])
    }
    samples.push(mean_doubles(left_sample) - mean_doubles(right_sample))
  }
  {
    estimate,
    standard_error: double_sqrt(population_variance(samples)),
    lower: quantile_doubles(samples, 0.025).unwrap_or(estimate),
    upper: quantile_doubles(samples, 0.975).unwrap_or(estimate),
    replicates,
  }
}
