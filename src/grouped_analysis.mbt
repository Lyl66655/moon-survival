///|
/// A grouped row combining a label and a survival observation.
pub(all) struct GroupedObservation {
  group : String
  time : Double
  event : Int
} derive(Eq, Debug)

///|
/// Summary of one group in a stratified analysis.
pub(all) struct GroupSurvivalSummary {
  group : String
  sample_size : Int
  events : Int
  censoring_fraction : Double
  median : Double?
  restricted_mean : Double
  survival_at_horizon : Double
} derive(Eq, Debug)

///|
/// Fit one curve per group label while retaining the original observations.
pub fn grouped_survival_summary(
  groups : Array[String],
  times : Array[Double],
  events : Array[Int],
  horizon : Double,
) -> Array[GroupSurvivalSummary] {
  if groups.length() != times.length() || groups.length() != events.length() {
    panic()
  }
  let unique : Array[String] = []
  for group in groups {
    if !unique.contains(group) {
      unique.push(group)
    }
  }
  let result = []
  for group in unique {
    let group_times = []
    let group_events = []
    for i in 0..<groups.length() {
      if groups[i] == group {
        group_times.push(times[i])
        group_events.push(events[i])
      }
    }
    let data = SurvivalData::new(group_times, group_events)
    let km = kaplan_meier(data)
    let censoring = data.censored_count()
    result.push({
      group,
      sample_size: data.length(),
      events: data.event_count(),
      censoring_fraction: if data.length() == 0 {
        0.0
      } else {
        censoring.to_double() / data.length().to_double()
      },
      median: km.median(),
      restricted_mean: km.restricted_mean(horizon),
      survival_at_horizon: km.predict(horizon),
    })
  }
  result
}

///|
/// Fixed-effect inverse-variance meta-analysis result.
pub(all) struct MetaAnalysisResult {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  z_value : Double
  p_value : Double
  heterogeneity_q : Double
  studies : Int
} derive(Eq, Debug)

///|
/// Combine log hazard ratios with inverse-variance weights.
pub fn fixed_effect_meta_analysis(
  estimates : Array[Double],
  standard_errors : Array[Double],
) -> MetaAnalysisResult {
  if estimates.length() != standard_errors.length() || estimates.length() == 0 {
    panic()
  }
  let mut weight_sum = 0.0
  let mut weighted_estimate = 0.0
  for i in 0..<estimates.length() {
    let variance = standard_errors[i] * standard_errors[i]
    if variance > survival_epsilon {
      let weight = 1.0 / variance
      weight_sum += weight
      weighted_estimate += weight * estimates[i]
    }
  }
  let estimate = if weight_sum <= survival_epsilon {
    0.0
  } else {
    weighted_estimate / weight_sum
  }
  let standard_error = if weight_sum <= survival_epsilon {
    0.0
  } else {
    double_sqrt(1.0 / weight_sum)
  }
  let mut heterogeneity = 0.0
  for i in 0..<estimates.length() {
    let variance = standard_errors[i] * standard_errors[i]
    if variance > survival_epsilon {
      heterogeneity += (estimates[i] - estimate) *
        (estimates[i] - estimate) /
        variance
    }
  }
  let z_value = if standard_error <= survival_epsilon {
    0.0
  } else {
    estimate / standard_error
  }
  {
    estimate,
    standard_error,
    lower: estimate - 1.96 * standard_error,
    upper: estimate + 1.96 * standard_error,
    z_value,
    p_value: 2.0 * (1.0 - normal_cdf(z_value.abs())),
    heterogeneity_q: heterogeneity,
    studies: estimates.length(),
  }
}

///|
/// Convert grouped summaries into CSV for cohort reporting.
pub fn grouped_summary_to_csv(
  summaries : Array[GroupSurvivalSummary],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "group,sample_size,events,censoring_fraction,restricted_mean,survival_at_horizon\n",
  )
  for summary in summaries {
    builder.write_string(
      "\{summary.group},\{summary.sample_size},\{summary.events},\{summary.censoring_fraction},\{summary.restricted_mean},\{summary.survival_at_horizon}\n",
    )
  }
  builder.to_string()
}

///|
/// Convert a log hazard ratio to a hazard ratio with confidence limits.
pub fn meta_hazard_ratio(
  result : MetaAnalysisResult,
) -> (Double, Double, Double) {
  (
    double_exp(result.estimate),
    double_exp(result.lower),
    double_exp(result.upper),
  )
}
