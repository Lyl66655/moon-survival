///|
/// Follow-up and administrative censoring utilities.
pub(all) struct SurvivalFollowupRow {
  row : Int
  observed_time : Double
  analysis_time : Double
  event : Int
  administratively_censored : Bool
  truncated : Bool
  exposure : Double
  followup_row : Int
} derive(Eq, Debug)

///|
pub(all) struct SurvivalFollowupSummary {
  sample_size : Int
  events : Int
  administrative_censorings : Int
  truncated_rows : Int
  total_observed_time : Double
  total_analysis_time : Double
  mean_followup : Double
  maximum_followup : Double
  horizon : Double
  followup_row : Int
} derive(Eq, Debug)

///|
/// Apply an administrative horizon to a cohort.
pub fn survival_apply_followup_horizon(
  data : SurvivalData,
  horizon : Double,
) -> Array[SurvivalFollowupRow] {
  let output : Array[SurvivalFollowupRow] = []
  let safe_horizon = horizon.max(0.0)
  for i = 0; i < data.length(); i = i + 1 {
    let observation = data.get(i)
    let truncated = observation.time > safe_horizon
    let analysis_time = if truncated { safe_horizon } else { observation.time }
    let event = if truncated { 0 } else { observation.event }
    output.push({
      row: i,
      observed_time: observation.time,
      analysis_time,
      event,
      administratively_censored: truncated,
      truncated,
      exposure: analysis_time,
      followup_row: i,
    })
  }
  output
}

///|
/// Convert a follow-up table back into validated SurvivalData.
pub fn survival_followup_to_data(
  rows : Array[SurvivalFollowupRow],
) -> SurvivalData {
  let times : Array[Double] = []
  let events : Array[Int] = []
  for row in rows {
    times.push(row.analysis_time)
    events.push(row.event)
  }
  SurvivalData::new(times, events)
}

///|
pub fn survival_followup_summary(
  rows : Array[SurvivalFollowupRow],
  horizon : Double,
) -> SurvivalFollowupSummary {
  let mut events = 0
  let mut administrative = 0
  let mut truncated = 0
  let mut observed_time = 0.0
  let mut analysis_time = 0.0
  let mut maximum = 0.0
  for row in rows {
    if row.event == 1 {
      events += 1
    }
    if row.administratively_censored {
      administrative += 1
    }
    if row.truncated {
      truncated += 1
    }
    observed_time += row.observed_time
    analysis_time += row.analysis_time
    if row.analysis_time > maximum {
      maximum = row.analysis_time
    }
  }
  let mean = if rows.length() == 0 {
    0.0
  } else {
    analysis_time / rows.length().to_double()
  }
  {
    sample_size: rows.length(),
    events,
    administrative_censorings: administrative,
    truncated_rows: truncated,
    total_observed_time: observed_time,
    total_analysis_time: analysis_time,
    mean_followup: mean,
    maximum_followup: maximum,
    horizon,
    followup_row: rows.length(),
  }
}

///|
/// Compute person-time in a closed-open interval.
pub fn survival_followup_exposure_between(
  rows : Array[SurvivalFollowupRow],
  start : Double,
  end : Double,
) -> Double {
  let mut total = 0.0
  for row in rows {
    let left = row.analysis_time.min(end).max(start)
    if row.analysis_time > start && end > start {
      total += (row.analysis_time.min(end) - start).max(0.0)
    }
    if left < start {
      total += 0.0
    }
  }
  total
}

///|
/// Count events observed before an administrative horizon.
pub fn survival_followup_event_count(
  rows : Array[SurvivalFollowupRow],
  horizon : Double,
) -> Int {
  let mut count = 0
  for row in rows {
    if row.event == 1 && row.analysis_time <= horizon {
      count += 1
    }
  }
  count
}

///|
pub fn survival_followup_to_csv(rows : Array[SurvivalFollowupRow]) -> String {
  let mut out = "row,observed_time,analysis_time,event,administratively_censored,truncated,exposure,followup_row\n"
  for row in rows {
    out = out +
      row.row.to_string() +
      "," +
      row.observed_time.to_string() +
      "," +
      row.analysis_time.to_string() +
      "," +
      row.event.to_string() +
      "," +
      row.administratively_censored.to_string() +
      "," +
      row.truncated.to_string() +
      "," +
      row.exposure.to_string() +
      "," +
      row.followup_row.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_followup_is_valid(rows : Array[SurvivalFollowupRow]) -> Bool {
  for row in rows {
    if row.observed_time < 0.0 ||
      row.analysis_time < 0.0 ||
      row.analysis_time > row.observed_time ||
      row.event < 0 ||
      row.event > 1 ||
      row.exposure < 0.0 {
      return false
    }
  }
  true
}

///|
pub fn survival_followup_summary_is_valid(
  summary : SurvivalFollowupSummary,
) -> Bool {
  summary.sample_size >= 0 &&
  summary.events >= 0 &&
  summary.events <= summary.sample_size &&
  summary.administrative_censorings >= 0 &&
  summary.truncated_rows >= 0 &&
  summary.mean_followup >= 0.0
}
