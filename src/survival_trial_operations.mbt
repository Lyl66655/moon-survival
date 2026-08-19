///|
/// Operational summaries for two-arm time-to-event studies.
pub(all) struct TrialOperationsArmSummary {
  arm_row : Int
  label : String
  sample_size : Int
  events : Int
  censored : Int
  event_fraction : Double
  median_follow_up : Double
  mean_follow_up : Double
  trial_row : Int
} derive(Eq, Debug)

///|
pub(all) struct TrialEnrollmentRow {
  period_row : Int
  period_start : Double
  period_end : Double
  planned : Int
  enrolled : Int
  cumulative : Int
  completion_fraction : Double
  enrollment_rate : Double
} derive(Eq, Debug)

///|
pub(all) struct TrialEventTarget {
  target_row : Int
  information_fraction : Double
  target_events : Int
  target_sample_size : Int
  expected_power : Double
  achieved : Bool
} derive(Eq, Debug)

///|
/// Summarize one arm of a survival trial.
pub fn survival_trial_arm_summary(
  arm_row : Int,
  label : String,
  data : SurvivalData,
) -> TrialOperationsArmSummary {
  let times : Array[Double] = []
  for observation in data.to_array() {
    times.push(observation.time)
  }
  times.sort()
  let mean = if times.length() == 0 { 0.0 } else { survival_ext_mean(times) }
  let median = if times.length() == 0 { 0.0 } else { times[times.length() / 2] }
  let sample = data.length()
  {
    arm_row,
    label,
    sample_size: sample,
    events: data.event_count(),
    censored: data.censored_count(),
    event_fraction: if sample == 0 {
      0.0
    } else {
      data.event_count().to_double() / sample.to_double()
    },
    median_follow_up: median,
    mean_follow_up: mean,
    trial_row: arm_row,
  }
}

///|
/// Build cumulative enrollment counts from planned and observed period counts.
pub fn survival_trial_enrollment(
  period_lengths : Array[Double],
  planned : Array[Int],
  enrolled : Array[Int],
) -> Array[TrialEnrollmentRow] {
  let output : Array[TrialEnrollmentRow] = []
  let count = period_lengths
    .length()
    .min(planned.length())
    .min(enrolled.length())
  let mut cumulative = 0
  let mut elapsed = 0.0
  for i = 0; i < count; i = i + 1 {
    let start = elapsed
    let length = period_lengths[i].max(0.0)
    elapsed += length
    cumulative += enrolled[i].max(0)
    let target = planned[i].max(0)
    let completion = if target == 0 {
      0.0
    } else {
      enrolled[i].max(0).to_double() / target.to_double()
    }
    let rate = if length == 0.0 {
      0.0
    } else {
      enrolled[i].max(0).to_double() / length
    }
    output.push({
      period_row: i,
      period_start: start,
      period_end: elapsed,
      planned: target,
      enrolled: enrolled[i].max(0),
      cumulative,
      completion_fraction: completion,
      enrollment_rate: rate,
    })
  }
  output
}

///|
/// Construct event targets for a staged interim monitoring plan.
pub fn survival_trial_event_targets(
  total_events : Int,
  total_sample_size : Int,
  information_fractions : Array[Double],
  baseline_power : Double,
) -> Array[TrialEventTarget] {
  let output : Array[TrialEventTarget] = []
  if total_events <= 0 || total_sample_size <= 0 {
    return output
  }
  for i = 0; i < information_fractions.length(); i = i + 1 {
    let fraction = clamp_probability(information_fractions[i])
    let target_events = (total_events.to_double() * fraction).ceil().to_int()
    let target_sample = (total_sample_size.to_double() * fraction)
      .ceil()
      .to_int()
    let power = clamp_probability(baseline_power * (0.5 + 0.5 * fraction))
    output.push({
      target_row: i,
      information_fraction: fraction,
      target_events,
      target_sample_size: target_sample,
      expected_power: power,
      achieved: false,
    })
  }
  output
}

///|
/// Mark target rows reached by the current observed event and sample counts.
pub fn survival_trial_update_targets(
  targets : Array[TrialEventTarget],
  observed_events : Int,
  observed_sample_size : Int,
) -> Array[TrialEventTarget] {
  let output : Array[TrialEventTarget] = []
  for target in targets {
    output.push({
      ..target,
      achieved: observed_events >= target.target_events &&
      observed_sample_size >= target.target_sample_size,
    })
  }
  output
}

///|
/// Compute the treatment allocation ratio for completed enrollment.
pub fn survival_trial_allocation_ratio(treated : Int, control : Int) -> Double {
  if control <= 0 {
    0.0
  } else {
    treated.to_double() / control.to_double()
  }
}

///|
pub fn survival_trial_enrollment_to_csv(
  rows : Array[TrialEnrollmentRow],
) -> String {
  let mut out = "period_row,period_start,period_end,planned,enrolled,cumulative,completion_fraction,enrollment_rate\n"
  for row in rows {
    out = out +
      row.period_row.to_string() +
      "," +
      row.period_start.to_string() +
      "," +
      row.period_end.to_string() +
      "," +
      row.planned.to_string() +
      "," +
      row.enrolled.to_string() +
      "," +
      row.cumulative.to_string() +
      "," +
      row.completion_fraction.to_string() +
      "," +
      row.enrollment_rate.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_trial_targets_to_csv(
  targets : Array[TrialEventTarget],
) -> String {
  let mut out = "target_row,information_fraction,target_events,target_sample_size,expected_power,achieved\n"
  for target in targets {
    out = out +
      target.target_row.to_string() +
      "," +
      target.information_fraction.to_string() +
      "," +
      target.target_events.to_string() +
      "," +
      target.target_sample_size.to_string() +
      "," +
      target.expected_power.to_string() +
      "," +
      target.achieved.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_trial_enrollment_is_valid(
  rows : Array[TrialEnrollmentRow],
) -> Bool {
  for row in rows {
    if row.period_end < row.period_start ||
      row.planned < 0 ||
      row.enrolled < 0 ||
      row.cumulative < row.enrolled ||
      row.completion_fraction < 0.0 ||
      row.enrollment_rate < 0.0 {
      return false
    }
  }
  true
}

///|
pub fn survival_trial_targets_are_valid(
  targets : Array[TrialEventTarget],
) -> Bool {
  for target in targets {
    if target.information_fraction < 0.0 ||
      target.information_fraction > 1.0 ||
      target.target_events < 0 ||
      target.target_sample_size < 0 ||
      target.expected_power < 0.0 ||
      target.expected_power > 1.0 {
      return false
    }
  }
  true
}
