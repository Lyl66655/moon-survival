///|
/// Diagnostics for the censoring process in right-censored cohorts.
pub(all) struct CensoringTimeRow {
  row : Int
  time : Double
  at_risk : Int
  events : Int
  censored : Int
  censoring_rate : Double
  event_rate : Double
} derive(Eq, Debug)

///|
pub(all) struct CensoringSummary {
  sample_size : Int
  events : Int
  censored : Int
  censoring_fraction : Double
  median_censoring_time : Double
  last_observed_time : Double
  censoring_curve : KMResult
} derive(Eq, Debug)

///|
/// Build an event table for the censoring process itself.
pub fn survival_censoring_table(data : SurvivalData) -> Array[CensoringTimeRow] {
  let output : Array[CensoringTimeRow] = []
  let observations = data.to_array()
  let mut index = 0
  let mut row = 0
  while index < observations.length() {
    let time = observations[index].time
    let mut events = 0
    let mut censored = 0
    while index < observations.length() && observations[index].time == time {
      if observations[index].event == 1 {
        events += 1
      } else {
        censored += 1
      }
      index += 1
    }
    let mut at_risk = 0
    for observation in observations {
      if observation.time >= time {
        at_risk += 1
      }
    }
    let censoring_rate = if at_risk == 0 {
      0.0
    } else {
      censored.to_double() / at_risk.to_double()
    }
    let event_rate = if at_risk == 0 {
      0.0
    } else {
      events.to_double() / at_risk.to_double()
    }
    output.push({
      row,
      time,
      at_risk,
      events,
      censored,
      censoring_rate,
      event_rate,
    })
    row += 1
  }
  output
}

///|
/// Estimate the censoring survival distribution G(t).
pub fn survival_censoring_curve(data : SurvivalData) -> KMResult {
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in data.to_array() {
    times.push(observation.time)
    events.push(if observation.event == 0 { 1 } else { 0 })
  }
  kaplan_meier(SurvivalData::new(times, events))
}

///|
/// Summarize the censoring burden and its estimated distribution.
pub fn survival_censoring_summary(data : SurvivalData) -> CensoringSummary {
  let curve = survival_censoring_curve(data)
  let censoring_times : Array[Double] = []
  for observation in data.to_array() {
    if observation.event == 0 {
      censoring_times.push(observation.time)
    }
  }
  censoring_times.sort()
  let median = if censoring_times.length() == 0 {
    -1.0
  } else {
    censoring_times[censoring_times.length() / 2]
  }
  let sample_size = data.length()
  let censored = data.censored_count()
  let fraction = if sample_size == 0 {
    0.0
  } else {
    censored.to_double() / sample_size.to_double()
  }
  let last_time = if sample_size == 0 {
    0.0
  } else {
    data.get(sample_size - 1).time
  }
  {
    sample_size,
    events: data.event_count(),
    censored,
    censoring_fraction: fraction,
    median_censoring_time: median,
    last_observed_time: last_time,
    censoring_curve: curve,
  }
}

///|
/// Estimate the censoring probability at a time using a left limit.
pub fn survival_censoring_probability(
  data : SurvivalData,
  time : Double,
) -> Double {
  let curve = survival_censoring_curve(data)
  curve.survival_at(time - survival_epsilon)
}

///|
/// Return the largest time at which the censoring curve remains above a limit.
pub fn survival_censoring_support(
  data : SurvivalData,
  probability_limit : Double,
) -> Double {
  let curve = survival_censoring_curve(data)
  let limit = clamp_probability(probability_limit)
  let mut support = 0.0
  for step in curve.steps {
    if step.survival >= limit {
      support = step.time
    }
  }
  support
}

///|
pub fn survival_censoring_table_to_csv(
  rows : Array[CensoringTimeRow],
) -> String {
  let mut out = "row,time,at_risk,events,censored,censoring_rate,event_rate\n"
  for row in rows {
    out = out +
      row.row.to_string() +
      "," +
      row.time.to_string() +
      "," +
      row.at_risk.to_string() +
      "," +
      row.events.to_string() +
      "," +
      row.censored.to_string() +
      "," +
      row.censoring_rate.to_string() +
      "," +
      row.event_rate.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_censoring_table_is_valid(
  rows : Array[CensoringTimeRow],
) -> Bool {
  for row in rows {
    if row.at_risk < 0 ||
      row.events < 0 ||
      row.censored < 0 ||
      row.events + row.censored > row.at_risk ||
      row.censoring_rate < 0.0 ||
      row.censoring_rate > 1.0 ||
      row.event_rate < 0.0 ||
      row.event_rate > 1.0 {
      return false
    }
  }
  true
}

///|
pub fn survival_censoring_summary_is_valid(summary : CensoringSummary) -> Bool {
  summary.sample_size >= 0 &&
  summary.events >= 0 &&
  summary.censored >= 0 &&
  summary.events + summary.censored == summary.sample_size &&
  summary.censoring_fraction >= 0.0 &&
  summary.censoring_fraction <= 1.0
}
