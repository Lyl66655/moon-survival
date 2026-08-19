///|
/// Patient-level event timelines used to prepare longitudinal survival data.
pub(all) struct SurvivalTimelineEntry {
  patient_row : Int
  entry_row : Int
  start : Double
  end : Double
  state : String
  event : Int
  exposure : Double
  timeline_row : Int
} derive(Eq, Debug)

///|
pub(all) struct SurvivalPatientTimeline {
  patient_row : Int
  entries : Array[SurvivalTimelineEntry]
  total_exposure : Double
  event_count : Int
  last_time : Double
  timeline_row : Int
} derive(Eq, Debug)

///|
/// Create a timeline from start/end/state arrays for one patient.
pub fn survival_patient_timeline(
  patient_row : Int,
  starts : Array[Double],
  ends : Array[Double],
  states : Array[String],
  events : Array[Int],
) -> SurvivalPatientTimeline {
  let entries : Array[SurvivalTimelineEntry] = []
  let count = starts
    .length()
    .min(ends.length())
    .min(states.length())
    .min(events.length())
  let mut exposure = 0.0
  let mut event_count = 0
  let mut last_time = 0.0
  for i = 0; i < count; i = i + 1 {
    let start = starts[i]
    let end = ends[i].max(start)
    let event = if events[i] == 1 { 1 } else { 0 }
    let duration = end - start
    exposure += duration
    if event == 1 {
      event_count += 1
    }
    if end > last_time {
      last_time = end
    }
    entries.push({
      patient_row,
      entry_row: i,
      start,
      end,
      state: states[i],
      event,
      exposure: duration,
      timeline_row: i,
    })
  }
  {
    patient_row,
    entries,
    total_exposure: exposure,
    event_count,
    last_time,
    timeline_row: patient_row,
  }
}

///|
/// Merge adjacent entries when their state and event status agree.
pub fn survival_timeline_merge(
  entries : Array[SurvivalTimelineEntry],
) -> Array[SurvivalTimelineEntry] {
  if entries.length() == 0 {
    return []
  }
  let output : Array[SurvivalTimelineEntry] = []
  for entry in entries {
    if output.length() == 0 {
      output.push(entry)
    } else {
      let last = output[output.length() - 1]
      if last.state == entry.state &&
        last.event == entry.event &&
        last.end <= entry.start + survival_epsilon {
        let merged_end = if entry.end > last.end { entry.end } else { last.end }
        let merged = {
          ..last,
          end: merged_end,
          exposure: merged_end - last.start,
        }
        output[output.length() - 1] = merged
      } else {
        output.push(entry)
      }
    }
  }
  output
}

///|
/// Restrict a patient timeline to a follow-up window.
pub fn survival_timeline_clip(
  entries : Array[SurvivalTimelineEntry],
  start : Double,
  end : Double,
) -> Array[SurvivalTimelineEntry] {
  let output : Array[SurvivalTimelineEntry] = []
  if end <= start {
    return output
  }
  for entry in entries {
    let left = if entry.start > start { entry.start } else { start }
    let right = if entry.end < end { entry.end } else { end }
    if right > left {
      output.push({ ..entry, start: left, end: right, exposure: right - left })
    }
  }
  output
}

///|
/// Convert a timeline into a right-censored endpoint at the supplied horizon.
pub fn survival_timeline_to_endpoint(
  entries : Array[SurvivalTimelineEntry],
  horizon : Double,
) -> Observation {
  let mut time = 0.0
  let mut event = 0
  for entry in entries {
    if entry.end <= horizon && entry.end > time {
      time = entry.end
    }
    if entry.event == 1 && entry.end <= horizon {
      event = 1
      time = entry.end
    }
  }
  let endpoint_time = if time <= 0.0 { horizon } else { time }
  { time: endpoint_time, event }
}

///|
pub fn survival_timeline_exposure(
  entries : Array[SurvivalTimelineEntry],
  state : String,
) -> Double {
  let mut exposure = 0.0
  for entry in entries {
    if entry.state == state {
      exposure += entry.exposure
    }
  }
  exposure
}

///|
pub fn survival_timeline_events(
  entries : Array[SurvivalTimelineEntry],
  state : String,
) -> Int {
  let mut count = 0
  for entry in entries {
    if entry.state == state && entry.event == 1 {
      count += 1
    }
  }
  count
}

///|
pub fn survival_timeline_to_csv(
  entries : Array[SurvivalTimelineEntry],
) -> String {
  let mut out = "patient_row,entry_row,start,end,state,event,exposure,timeline_row\n"
  for entry in entries {
    out = out +
      entry.patient_row.to_string() +
      "," +
      entry.entry_row.to_string() +
      "," +
      entry.start.to_string() +
      "," +
      entry.end.to_string() +
      "," +
      entry.state +
      "," +
      entry.event.to_string() +
      "," +
      entry.exposure.to_string() +
      "," +
      entry.timeline_row.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_timeline_is_valid(
  entries : Array[SurvivalTimelineEntry],
) -> Bool {
  for entry in entries {
    if entry.start < 0.0 ||
      entry.end < entry.start ||
      entry.event < 0 ||
      entry.event > 1 ||
      entry.exposure < 0.0 {
      return false
    }
  }
  true
}

///|
pub fn survival_patient_timeline_is_valid(
  timeline : SurvivalPatientTimeline,
) -> Bool {
  timeline.total_exposure >= 0.0 &&
  timeline.event_count >= 0 &&
  timeline.event_count <= timeline.entries.length() &&
  survival_timeline_is_valid(timeline.entries)
}
