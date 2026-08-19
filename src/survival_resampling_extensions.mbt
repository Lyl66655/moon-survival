///|
/// Reproducible resampling utilities for survival cohorts.
pub(all) struct SurvivalResamplePlan {
  sample_size : Int
  draws : Int
  seed : UInt
  replicate_row : Int
  replacement : Bool
} derive(Eq, Debug)

///|
pub(all) struct SurvivalResampleSummary {
  replicate_row : Int
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  successful_replicates : Int
  requested_replicates : Int
} derive(Eq, Debug)

///|
/// Sample a survival cohort with replacement using the project RNG.
pub fn survival_resample_with_replacement(
  data : SurvivalData,
  draws : Int,
  seed : UInt,
) -> SurvivalData {
  if data.length() == 0 || draws <= 0 {
    return SurvivalData::new([], [])
  }
  let observations = data.to_array()
  let rng = DeterministicRng::new(seed)
  let times : Array[Double] = []
  let events : Array[Int] = []
  for _ in 0..<draws {
    let index = (rng.next_u32() % observations.length().reinterpret_as_uint()).reinterpret_as_int()
    times.push(observations[index].time)
    events.push(observations[index].event)
  }
  SurvivalData::new(times, events)
}

///|
/// Build independent replicate seeds from one study seed.
pub fn survival_resample_seeds(seed : UInt, replicates : Int) -> Array[UInt] {
  let output : Array[UInt] = []
  if replicates <= 0 {
    return output
  }
  let rng = DeterministicRng::new(seed)
  for _ in 0..<replicates {
    output.push(rng.next_u32())
  }
  output
}

///|
/// Bootstrap the difference in RMST between two arms.
pub fn survival_resample_rmst_difference(
  left : SurvivalData,
  right : SurvivalData,
  horizon : Double,
  replicates : Int,
  seed : UInt,
) -> SurvivalResampleSummary {
  let observed = kaplan_meier(left).restricted_mean(horizon) -
    kaplan_meier(right).restricted_mean(horizon)
  if replicates <= 0 || left.length() == 0 || right.length() == 0 {
    return {
      replicate_row: 0,
      estimate: observed,
      standard_error: 0.0,
      lower: observed,
      upper: observed,
      successful_replicates: 0,
      requested_replicates: replicates.max(0),
    }
  }
  let seeds = survival_resample_seeds(seed, replicates)
  let estimates : Array[Double] = []
  for replicate_seed in seeds {
    let left_sample = survival_resample_with_replacement(
      left,
      left.length(),
      replicate_seed,
    )
    let right_sample = survival_resample_with_replacement(
      right,
      right.length(),
      replicate_seed + 17,
    )
    estimates.push(
      kaplan_meier(left_sample).restricted_mean(horizon) -
      kaplan_meier(right_sample).restricted_mean(horizon),
    )
  }
  let standard_error = double_sqrt(population_variance(estimates))
  {
    replicate_row: 0,
    estimate: observed,
    standard_error,
    lower: quantile_doubles(estimates, 0.025).unwrap_or(observed),
    upper: quantile_doubles(estimates, 0.975).unwrap_or(observed),
    successful_replicates: estimates.length(),
    requested_replicates: replicates,
  }
}

///|
/// Resample every group independently while keeping its original size.
pub fn survival_resample_groups(
  groups : Array[SurvivalData],
  seed : UInt,
) -> Array[SurvivalData] {
  let output : Array[SurvivalData] = []
  let seeds = survival_resample_seeds(seed, groups.length())
  for i = 0; i < groups.length(); i = i + 1 {
    output.push(
      survival_resample_with_replacement(
        groups[i],
        groups[i].length(),
        seeds[i],
      ),
    )
  }
  output
}

///|
/// Return a deterministic plan describing a resampling job.
pub fn survival_resample_plan(
  sample_size : Int,
  draws : Int,
  seed : UInt,
  replacement : Bool,
) -> SurvivalResamplePlan {
  {
    sample_size: sample_size.max(0),
    draws: draws.max(0),
    seed,
    replicate_row: 0,
    replacement,
  }
}

///|
pub fn survival_resample_summary_to_csv(
  summary : SurvivalResampleSummary,
) -> String {
  "replicate_row,estimate,standard_error,lower,upper,successful_replicates,requested_replicates\n" +
  summary.replicate_row.to_string() +
  "," +
  summary.estimate.to_string() +
  "," +
  summary.standard_error.to_string() +
  "," +
  summary.lower.to_string() +
  "," +
  summary.upper.to_string() +
  "," +
  summary.successful_replicates.to_string() +
  "," +
  summary.requested_replicates.to_string() +
  "\n"
}

///|
pub fn survival_resample_plan_is_valid(plan : SurvivalResamplePlan) -> Bool {
  plan.sample_size >= 0 && plan.draws >= 0
}

///|
pub fn survival_resample_summary_is_valid(
  summary : SurvivalResampleSummary,
) -> Bool {
  summary.successful_replicates >= 0 &&
  summary.requested_replicates >= summary.successful_replicates &&
  summary.standard_error >= 0.0 &&
  summary.lower <= summary.upper
}
