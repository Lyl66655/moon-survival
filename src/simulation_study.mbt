///|
/// A simulation-study design for survival-method validation.
pub(all) struct SurvivalSimulationStudy {
  name : String
  replicates : Int
  sample_size : Int
  baseline_rate : Double
  censor_rate : Double
  hazard_ratio : Double
  alpha : Double
  seed : UInt
  study_row : Bool
} derive(Eq, Debug)

///|
/// One replicate-level operating-characteristic record.
pub(all) struct SimulationReplicateResult {
  replicate : Int
  sample_size : Int
  events : Int
  estimate : Double
  standard_error : Double
  p_value : Double
  rejected : Bool
  simulation_row : Bool
} derive(Eq, Debug)

///|
/// Aggregated operating characteristics.
pub(all) struct SimulationOperatingCharacteristics {
  replicates : Int
  mean_estimate : Double
  bias : Double
  monte_carlo_se : Double
  empirical_power : Double
  mean_events : Double
  coverage : Double
  target : Double
} derive(Eq, Debug)

///|
/// Create a validated study design.
pub fn survival_simulation_study(
  name : String,
  replicates : Int,
  sample_size : Int,
  baseline_rate : Double,
  censor_rate : Double,
  hazard_ratio : Double,
  alpha : Double,
  seed : UInt,
) -> SurvivalSimulationStudy {
  {
    name,
    replicates: if replicates <= 0 {
      1
    } else {
      replicates
    },
    sample_size: if sample_size <= 1 {
      2
    } else {
      sample_size
    },
    baseline_rate: clamp_positive(baseline_rate),
    censor_rate: if censor_rate < 0.0 {
      0.0
    } else {
      censor_rate
    },
    hazard_ratio: clamp_positive(hazard_ratio),
    alpha: if alpha <= 0.0 {
      0.000001
    } else if alpha >= 1.0 {
      0.999999
    } else {
      alpha
    },
    seed,
    study_row: true,
  }
}

///|
/// Simulate one two-arm proportional-hazards study.
pub fn simulate_two_arm_study(
  study : SurvivalSimulationStudy,
  replicate : Int,
) -> (SurvivalData, SurvivalData) {
  let rng = DeterministicRng::new(study.seed + replicate.reinterpret_as_uint())
  let first_times : Array[Double] = []
  let first_events : Array[Int] = []
  let second_times : Array[Double] = []
  let second_events : Array[Int] = []
  let first_size = study.sample_size / 2
  for arm in 0..<2 {
    let size = if arm == 0 {
      first_size
    } else {
      study.sample_size - first_size
    }
    for _ in 0..<size {
      let arm_rate = if arm == 0 { 1.0 } else { study.hazard_ratio }
      let failure = -double_log(clamp_positive(rng.uniform(), floor=0.0000001)) /
        (study.baseline_rate * arm_rate)
      let censor = if study.censor_rate <= 0.0 {
        failure + 1000000.0
      } else {
        -double_log(clamp_positive(rng.uniform(), floor=0.0000001)) /
        study.censor_rate
      }
      let time = if failure <= censor { failure } else { censor }
      let event = if failure <= censor { 1 } else { 0 }
      if arm == 0 {
        first_times.push(time)
        first_events.push(event)
      } else {
        second_times.push(time)
        second_events.push(event)
      }
    }
  }
  (
    SurvivalData::new(first_times, first_events),
    SurvivalData::new(second_times, second_events),
  )
}

///|
/// Run a log-rank replicate and retain a binary decision.
pub fn survival_simulation_replicate(
  study : SurvivalSimulationStudy,
  replicate : Int,
) -> SimulationReplicateResult {
  let (first, second) = simulate_two_arm_study(study, replicate)
  let result = logrank_test(first, second)
  let events = first.event_count() + second.event_count()
  let estimate = if result.observed.length() < 2 ||
    result.expected.length() < 2 ||
    result.expected[1] <= survival_epsilon {
    0.0
  } else {
    double_log(result.observed[1] / result.expected[1])
  }
  let standard_error = double_sqrt(
    if result.variance <= 0.0 {
      0.0
    } else {
      1.0 / result.variance
    },
  )
  {
    replicate,
    sample_size: first.length() + second.length(),
    events,
    estimate,
    standard_error,
    p_value: result.p_value,
    rejected: result.p_value < study.alpha,
    simulation_row: true,
  }
}

///|
/// Run all replicates in a reproducible study.
pub fn survival_simulation_run(
  study : SurvivalSimulationStudy,
) -> Array[SimulationReplicateResult] {
  let result : Array[SimulationReplicateResult] = []
  for replicate in 0..<study.replicates {
    result.push(survival_simulation_replicate(study, replicate))
  }
  result
}

///|
/// Aggregate bias, power, coverage and Monte Carlo uncertainty.
pub fn survival_simulation_summarize(
  study : SurvivalSimulationStudy,
  results : Array[SimulationReplicateResult],
) -> SimulationOperatingCharacteristics {
  if results.length() == 0 {
    return {
      replicates: 0,
      mean_estimate: 0.0,
      bias: 0.0,
      monte_carlo_se: 0.0,
      empirical_power: 0.0,
      mean_events: 0.0,
      coverage: 0.0,
      target: double_log(study.hazard_ratio),
    }
  }
  let target = double_log(study.hazard_ratio)
  let estimates : Array[Double] = []
  let mut rejected = 0
  let mut events = 0.0
  let mut covered = 0
  for result in results {
    estimates.push(result.estimate)
    if result.rejected {
      rejected += 1
    }
    events += result.events.to_double()
    let lower = result.estimate - 1.96 * result.standard_error
    let upper = result.estimate + 1.96 * result.standard_error
    if target >= lower && target <= upper {
      covered += 1
    }
  }
  let mean = mean_doubles(estimates)
  let mut variance = 0.0
  for estimate in estimates {
    let delta = estimate - mean
    variance += delta * delta
  }
  let mc_se = double_sqrt(
    variance / (estimates.length() * estimates.length()).to_double(),
  )
  {
    replicates: results.length(),
    mean_estimate: mean,
    bias: mean - target,
    monte_carlo_se: mc_se,
    empirical_power: rejected.to_double() / results.length().to_double(),
    mean_events: events / results.length().to_double(),
    coverage: covered.to_double() / results.length().to_double(),
    target,
  }
}

///|
/// Estimate an empirical quantile of replicate estimates.
pub fn survival_simulation_quantile(
  results : Array[SimulationReplicateResult],
  probability : Double,
) -> Double {
  if results.length() == 0 {
    return 0.0
  }
  let values : Array[Double] = []
  for result in results {
    values.push(result.estimate)
  }
  values.sort()
  let index = (clamp_probability(probability) *
  (values.length() - 1).to_double()).to_int()
  values[index]
}

///|
/// Estimate the probability of a sign-correct effect.
pub fn survival_simulation_sign_accuracy(
  results : Array[SimulationReplicateResult],
  target : Double,
) -> Double {
  if results.length() == 0 {
    return 0.0
  }
  let mut correct = 0
  for result in results {
    if result.estimate * target >= 0.0 {
      correct += 1
    }
  }
  correct.to_double() / results.length().to_double()
}

///|
/// Serialize replicate results.
pub fn simulation_replicates_to_csv(
  results : Array[SimulationReplicateResult],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "replicate,sample_size,events,estimate,standard_error,p_value,rejected\n",
  )
  for result in results {
    builder.write_string(result.replicate.to_string())
    builder.write_char(',')
    builder.write_string(result.sample_size.to_string())
    builder.write_char(',')
    builder.write_string(result.events.to_string())
    builder.write_char(',')
    builder.write_string(result.estimate.to_string())
    builder.write_char(',')
    builder.write_string(result.standard_error.to_string())
    builder.write_char(',')
    builder.write_string(result.p_value.to_string())
    builder.write_char(',')
    builder.write_string(result.rejected.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Serialize operating characteristics.
pub fn simulation_summary_to_csv(
  summary : SimulationOperatingCharacteristics,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("replicates,")
  builder.write_string(summary.replicates.to_string())
  builder.write_char('\n')
  builder.write_string("mean_estimate,")
  builder.write_string(summary.mean_estimate.to_string())
  builder.write_char('\n')
  builder.write_string("bias,")
  builder.write_string(summary.bias.to_string())
  builder.write_char('\n')
  builder.write_string("monte_carlo_se,")
  builder.write_string(summary.monte_carlo_se.to_string())
  builder.write_char('\n')
  builder.write_string("empirical_power,")
  builder.write_string(summary.empirical_power.to_string())
  builder.write_char('\n')
  builder.write_string("mean_events,")
  builder.write_string(summary.mean_events.to_string())
  builder.write_char('\n')
  builder.write_string("coverage,")
  builder.write_string(summary.coverage.to_string())
  builder.write_char('\n')
  builder.write_string("target,")
  builder.write_string(summary.target.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Validate an operating-characteristics summary.
pub fn simulation_summary_is_valid(
  summary : SimulationOperatingCharacteristics,
) -> Bool {
  summary.replicates >= 0 &&
  summary.empirical_power >= 0.0 &&
  summary.empirical_power <= 1.0 &&
  summary.coverage >= 0.0 &&
  summary.coverage <= 1.0 &&
  !summary.mean_estimate.is_nan()
}
