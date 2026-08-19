///|
/// Cohort summaries for operational time windows.
///
/// These helpers make it possible to turn a survival cohort into a set of
/// reproducible calendar-like windows without losing the event/censor counts.
pub(all) struct SurvivalTimeWindow {
  window_row : Int
  start : Double
  end : Double
  at_risk : Int
  events : Int
  censored : Int
  event_rate : Double
  exposure : Double
} derive(Eq, Debug)

///|
pub(all) struct SurvivalSplitResult {
  before : SurvivalData
  after : SurvivalData
  boundary : Double
  before_events : Int
  after_events : Int
  split_row : Int
} derive(Eq, Debug)

///|
fn survival_time_window_empty(
  row : Int,
  start : Double,
  end : Double,
) -> SurvivalTimeWindow {
  {
    window_row: row,
    start,
    end,
    at_risk: 0,
    events: 0,
    censored: 0,
    event_rate: 0.0,
    exposure: 0.0,
  }
}

///|
/// Aggregate observations into equal-width windows.
pub fn survival_time_windows(
  data : SurvivalData,
  start : Double,
  end : Double,
  width : Double,
) -> Array[SurvivalTimeWindow] {
  if end <= start || width <= 0.0 {
    return []
  }
  let count = ((end - start) / width).ceil().to_int()
  let output = Array::make(count, survival_time_window_empty(0, 0.0, 0.0))
  for row = 0; row < count; row = row + 1 {
    let left = start + row.to_double() * width
    let right = if left + width < end { left + width } else { end }
    output[row] = survival_time_window_empty(row, left, right)
  }
  for observation in data.to_array() {
    if observation.time >= start && observation.time <= end {
      let raw = ((observation.time - start) / width).floor().to_int()
      let row = if raw >= count { count - 1 } else { raw }
      if row >= 0 && row < count {
        let current = output[row]
        let event_increment = if observation.event == 1 { 1 } else { 0 }
        let censor_increment = if observation.event == 0 { 1 } else { 0 }
        let at_risk = current.at_risk + 1
        let events = current.events + event_increment
        let censored = current.censored + censor_increment
        let exposure = current.exposure +
          (observation.time - current.start).max(0.0)
        let event_rate = if exposure == 0.0 {
          0.0
        } else {
          events.to_double() / exposure
        }
        output[row] = {
          ..current,
          at_risk,
          events,
          censored,
          event_rate,
          exposure,
        }
      } else {
        ()
      }
    }
  }
  output
}

///|
/// Split a cohort at a time boundary while preserving valid SurvivalData.
pub fn survival_time_split(
  data : SurvivalData,
  boundary : Double,
) -> SurvivalSplitResult {
  let before_times = []
  let before_events = []
  let after_times = []
  let after_events = []
  let mut before_event_count = 0
  let mut after_event_count = 0
  for observation in data.to_array() {
    if observation.time <= boundary {
      before_times.push(observation.time)
      before_events.push(observation.event)
      if observation.event == 1 {
        before_event_count = before_event_count + 1
      }
    } else {
      after_times.push(observation.time - boundary)
      after_events.push(observation.event)
      if observation.event == 1 {
        after_event_count = after_event_count + 1
      }
    }
  }
  {
    before: SurvivalData::new(before_times, before_events),
    after: SurvivalData::new(after_times, after_events),
    boundary,
    before_events: before_event_count,
    after_events: after_event_count,
    split_row: data.length(),
  }
}

///|
/// Return the number of observations that are still at risk at a landmark.
pub fn survival_at_risk_at(data : SurvivalData, landmark : Double) -> Int {
  let mut count = 0
  for observation in data.to_array() {
    if observation.time >= landmark {
      count = count + 1
    }
  }
  count
}

///|
/// Count events in a closed-open interval.
pub fn survival_events_between(
  data : SurvivalData,
  start : Double,
  end : Double,
) -> Int {
  let mut count = 0
  for observation in data.to_array() {
    if observation.event == 1 &&
      observation.time >= start &&
      observation.time < end {
      count = count + 1
    }
  }
  count
}

///|
/// Count censorings in a closed-open interval.
pub fn survival_censors_between(
  data : SurvivalData,
  start : Double,
  end : Double,
) -> Int {
  let mut count = 0
  for observation in data.to_array() {
    if observation.event == 0 &&
      observation.time >= start &&
      observation.time < end {
      count = count + 1
    }
  }
  count
}

///|
pub fn survival_time_windows_to_csv(
  windows : Array[SurvivalTimeWindow],
) -> String {
  let mut out = "window_row,start,end,at_risk,events,censored,event_rate,exposure\n"
  for row in windows {
    out = out +
      row.window_row.to_string() +
      "," +
      row.start.to_string() +
      "," +
      row.end.to_string() +
      "," +
      row.at_risk.to_string() +
      "," +
      row.events.to_string() +
      "," +
      row.censored.to_string() +
      "," +
      row.event_rate.to_string() +
      "," +
      row.exposure.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_time_windows_is_valid(
  windows : Array[SurvivalTimeWindow],
) -> Bool {
  for row in windows {
    if row.start >= row.end ||
      row.at_risk < 0 ||
      row.events < 0 ||
      row.censored < 0 ||
      row.events + row.censored > row.at_risk ||
      row.event_rate < 0.0 ||
      row.exposure < 0.0 {
      return false
    }
  }
  true
}

///|
/// Produce a compact table of event and censor counts at requested landmarks.
pub fn survival_landmark_counts(
  data : SurvivalData,
  landmarks : Array[Double],
) -> Array[SurvivalTimeWindow] {
  let output = []
  for i = 0; i < landmarks.length(); i = i + 1 {
    let landmark = landmarks[i]
    let next = if i + 1 < landmarks.length() {
      landmarks[i + 1]
    } else {
      landmark + 1.0
    }
    let events = survival_events_between(data, landmark, next)
    let censored = survival_censors_between(data, landmark, next)
    let at_risk = survival_at_risk_at(data, landmark)
    let exposure = (next - landmark).max(0.0) * at_risk.to_double()
    let rate = if exposure == 0.0 { 0.0 } else { events.to_double() / exposure }
    output.push({
      window_row: i,
      start: landmark,
      end: next,
      at_risk,
      events,
      censored,
      event_rate: rate,
      exposure,
    })
  }
  output
}
