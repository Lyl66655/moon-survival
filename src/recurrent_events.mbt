///|
/// A start-stop record for recurrent-event analysis.
pub(all) struct RecurrentRecord {
  subject : Int
  start : Double
  stop : Double
  event : Int
  recurrent_row : Bool
} derive(Eq, Debug)

///|
/// Input validation failures for recurrent-event data.
pub(all) enum RecurrentError {
  EmptyData
  LengthMismatch
  InvalidSubject
  InvalidInterval
  InvalidEvent
  UnsortedSubject
  InvalidWindow
} derive(Eq, Debug)

///|
/// Validated start-stop recurrent-event data.
pub(all) struct RecurrentData {
  records : Array[RecurrentRecord]
  subject_count : Int
} derive(Eq, Debug)

///|
/// A recurrent-event risk-table row.
pub(all) struct RecurrentRiskStep {
  time : Double
  subjects_at_risk : Int
  intervals_at_risk : Int
  events : Int
  event_rate : Double
  cumulative_mean : Double
  standard_error : Double
} derive(Eq, Debug)

///|
/// Per-subject event and exposure summary.
pub(all) struct RecurrentSubjectSummary {
  subject : Int
  first_entry : Double
  last_exit : Double
  event_count : Int
  exposure : Double
  event_rate : Double
  recurrent_summary : Bool
} derive(Eq, Debug)

///|
/// A gap-time record generated from recurrent episodes.
pub(all) struct GapTimeRecord {
  subject : Int
  episode : Int
  gap_start : Double
  gap_stop : Double
  event : Int
  gap_row : Bool
} derive(Eq, Debug)

///|
/// A simple robust rate-ratio estimate between recurrent-event groups.
pub(all) struct RecurrentRateContrast {
  exposed_events : Int
  exposed_person_time : Double
  reference_events : Int
  reference_person_time : Double
  rate_ratio : Double
  log_rate_ratio : Double
  standard_error : Double
  lower_ci : Double
  upper_ci : Double
} derive(Eq, Debug)

///|
/// Construct recurrent data from parallel start-stop arrays.
pub fn RecurrentData::try_new(
  subjects : Array[Int],
  starts : Array[Double],
  stops : Array[Double],
  events : Array[Int],
) -> Result[RecurrentData, RecurrentError] {
  if subjects.length() == 0 {
    return Err(EmptyData)
  }
  if subjects.length() != starts.length() ||
    subjects.length() != stops.length() ||
    subjects.length() != events.length() {
    return Err(LengthMismatch)
  }
  let records : Array[RecurrentRecord] = []
  let mut maximum_subject = -1
  for i in 0..<subjects.length() {
    if subjects[i] < 0 {
      return Err(InvalidSubject)
    }
    if starts[i].is_nan() ||
      starts[i].is_inf() ||
      stops[i].is_nan() ||
      stops[i].is_inf() ||
      starts[i] < 0.0 ||
      stops[i] <= starts[i] {
      return Err(InvalidInterval)
    }
    if events[i] != 0 && events[i] != 1 {
      return Err(InvalidEvent)
    }
    if subjects[i] > maximum_subject {
      maximum_subject = subjects[i]
    }
    records.push({
      subject: subjects[i],
      start: starts[i],
      stop: stops[i],
      event: events[i],
      recurrent_row: true,
    })
  }
  records.sort_by(fn(a, b) {
    if a.subject < b.subject {
      -1
    } else if a.subject > b.subject {
      1
    } else if a.start < b.start {
      -1
    } else if a.start > b.start {
      1
    } else {
      a.stop.compare(b.stop)
    }
  })
  let mut previous_subject = -1
  let mut previous_stop = 0.0
  for record in records {
    if record.subject != previous_subject {
      previous_subject = record.subject
      previous_stop = record.start
    } else if record.start < previous_stop {
      return Err(UnsortedSubject)
    }
    previous_stop = record.stop
  }
  Ok({ records, subject_count: maximum_subject + 1 })
}

///|
/// Construct recurrent data and panic only for programmer-supplied invalid data.
pub fn RecurrentData::new(
  subjects : Array[Int],
  starts : Array[Double],
  stops : Array[Double],
  events : Array[Int],
) -> RecurrentData {
  match RecurrentData::try_new(subjects, starts, stops, events) {
    Ok(data) => data
    Err(_) => panic()
  }
}

///|
/// Number of interval records.
pub fn RecurrentData::length(self : RecurrentData) -> Int {
  self.records.length()
}

///|
/// Number of distinct subject identifiers represented by the maximum id.
pub fn RecurrentData::subject_count(self : RecurrentData) -> Int {
  self.subject_count
}

///|
/// Access a start-stop record.
pub fn RecurrentData::get(self : RecurrentData, index : Int) -> RecurrentRecord {
  self.records[index]
}

///|
/// Copy all records.
pub fn RecurrentData::to_array(self : RecurrentData) -> Array[RecurrentRecord] {
  let result : Array[RecurrentRecord] = []
  for record in self.records {
    result.push(record)
  }
  result
}

///|
/// Return the distinct event times.
pub fn recurrent_event_times(data : RecurrentData) -> Array[Double] {
  let result : Array[Double] = []
  for record in data.records {
    if record.event == 1 && !result.contains(record.stop) {
      result.push(record.stop)
    }
  }
  result.sort()
  result
}

///|
/// Calculate all-subject exposure and event totals.
pub fn recurrent_totals(data : RecurrentData) -> (Int, Double) {
  let mut events = 0
  let mut exposure = 0.0
  for record in data.records {
    events += record.event
    exposure += record.stop - record.start
  }
  (events, exposure)
}

///|
/// Create a subject-level summary table.
pub fn recurrent_subject_summaries(
  data : RecurrentData,
) -> Array[RecurrentSubjectSummary] {
  let result : Array[RecurrentSubjectSummary] = []
  for subject in 0..<data.subject_count {
    let mut seen = false
    let mut first_entry = 0.0
    let mut last_exit = 0.0
    let mut events = 0
    let mut exposure = 0.0
    for record in data.records {
      if record.subject == subject {
        if !seen {
          first_entry = record.start
          seen = true
        }
        if record.stop > last_exit {
          last_exit = record.stop
        }
        events += record.event
        exposure += record.stop - record.start
      }
    }
    if seen {
      result.push({
        subject,
        first_entry,
        last_exit,
        event_count: events,
        exposure,
        event_rate: if exposure <= survival_epsilon {
          0.0
        } else {
          events.to_double() / exposure
        },
        recurrent_summary: true,
      })
    }
  }
  result
}

///|
/// Estimate the mean cumulative function using the event-rate identity.
pub fn recurrent_mean_cumulative_function(
  data : RecurrentData,
) -> Array[RecurrentRiskStep] {
  let result : Array[RecurrentRiskStep] = []
  let times = recurrent_event_times(data)
  let mut cumulative = 0.0
  let mut variance = 0.0
  for time in times {
    let mut intervals = 0
    let mut events = 0
    let seen_subjects : Array[Int] = []
    for record in data.records {
      if record.start < time && record.stop >= time {
        intervals += 1
        if !seen_subjects.contains(record.subject) {
          seen_subjects.push(record.subject)
        }
      }
      if record.stop == time && record.event == 1 {
        events += 1
      }
    }
    let rate = if intervals == 0 {
      0.0
    } else {
      events.to_double() / intervals.to_double()
    }
    cumulative += rate
    if intervals > 0 {
      variance += events.to_double() /
        (intervals.to_double() * intervals.to_double())
    }
    result.push({
      time,
      subjects_at_risk: seen_subjects.length(),
      intervals_at_risk: intervals,
      events,
      event_rate: rate,
      cumulative_mean: cumulative,
      standard_error: double_sqrt(variance),
    })
  }
  result
}

///|
/// Estimate the mean cumulative count at a time horizon.
pub fn recurrent_cumulative_at(
  steps : Array[RecurrentRiskStep],
  horizon : Double,
) -> Double {
  let mut value = 0.0
  for step in steps {
    if step.time <= horizon {
      value = step.cumulative_mean
    } else {
      break
    }
  }
  value
}

///|
/// Convert records to gap-time episodes within each subject.
pub fn recurrent_to_gap_time(data : RecurrentData) -> Array[GapTimeRecord] {
  let result : Array[GapTimeRecord] = []
  let summaries = recurrent_subject_summaries(data)
  for summary in summaries {
    let mut episode = 0
    let mut origin = summary.first_entry
    for record in data.records {
      if record.subject == summary.subject {
        result.push({
          subject: record.subject,
          episode,
          gap_start: record.start - origin,
          gap_stop: record.stop - origin,
          event: record.event,
          gap_row: true,
        })
        if record.event == 1 {
          episode += 1
          origin = record.stop
        }
      }
    }
  }
  result
}

///|
/// Estimate a constant event rate for a recurrent dataset.
pub fn recurrent_event_rate(data : RecurrentData) -> Double {
  let (events, exposure) = recurrent_totals(data)
  if exposure <= survival_epsilon {
    0.0
  } else {
    events.to_double() / exposure
  }
}

///|
/// Estimate a rate ratio from two independent recurrent datasets.
pub fn recurrent_rate_contrast(
  exposed : RecurrentData,
  reference : RecurrentData,
) -> RecurrentRateContrast {
  let (exposed_events, exposed_time) = recurrent_totals(exposed)
  let (reference_events, reference_time) = recurrent_totals(reference)
  let exposed_rate = if exposed_time <= survival_epsilon {
    0.0
  } else {
    exposed_events.to_double() / exposed_time
  }
  let reference_rate = if reference_time <= survival_epsilon {
    0.0
  } else {
    reference_events.to_double() / reference_time
  }
  let ratio = if reference_rate <= survival_epsilon {
    0.0
  } else {
    exposed_rate / reference_rate
  }
  let standard_error = double_sqrt(
    if exposed_events > 0 && reference_events > 0 {
      1.0 / exposed_events.to_double() + 1.0 / reference_events.to_double()
    } else {
      0.0
    },
  )
  let log_ratio = if ratio <= survival_epsilon {
    0.0
  } else {
    double_log(ratio)
  }
  {
    exposed_events,
    exposed_person_time: exposed_time,
    reference_events,
    reference_person_time: reference_time,
    rate_ratio: ratio,
    log_rate_ratio: log_ratio,
    standard_error,
    lower_ci: double_exp(log_ratio - 1.96 * standard_error),
    upper_ci: double_exp(log_ratio + 1.96 * standard_error),
  }
}

///|
/// Count event episodes by subject.
pub fn recurrent_event_counts(data : RecurrentData) -> Array[Int] {
  let counts = Array::make(data.subject_count, 0)
  for record in data.records {
    counts[record.subject] += record.event
  }
  counts
}

///|
/// Return a truncated dataset ending at a follow-up horizon.
pub fn recurrent_truncate(
  data : RecurrentData,
  horizon : Double,
) -> RecurrentData {
  if horizon <= 0.0 {
    return RecurrentData::new([], [], [], [])
  }
  let subjects : Array[Int] = []
  let starts : Array[Double] = []
  let stops : Array[Double] = []
  let events : Array[Int] = []
  for record in data.records {
    if record.start >= horizon {
      continue
    }
    subjects.push(record.subject)
    stops.push(if record.stop > horizon { horizon } else { record.stop })
    starts.push(record.start)
    events.push(if record.stop <= horizon { record.event } else { 0 })
  }
  RecurrentData::new(subjects, starts, stops, events)
}

///|
/// Return the cumulative event count for each subject up to a horizon.
pub fn recurrent_counts_at(
  data : RecurrentData,
  horizon : Double,
) -> Array[Int] {
  let counts = Array::make(data.subject_count, 0)
  for record in data.records {
    if record.event == 1 && record.stop <= horizon {
      counts[record.subject] += 1
    }
  }
  counts
}

///|
/// Export the start-stop table.
pub fn RecurrentData::to_csv(self : RecurrentData) -> String {
  let builder = StringBuilder::new()
  builder.write_string("subject,start,stop,event\n")
  for record in self.records {
    builder.write_string(record.subject.to_string())
    builder.write_char(',')
    builder.write_string(record.start.to_string())
    builder.write_char(',')
    builder.write_string(record.stop.to_string())
    builder.write_char(',')
    builder.write_string(record.event.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export the mean cumulative function.
pub fn recurrent_steps_to_csv(steps : Array[RecurrentRiskStep]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,subjects_at_risk,intervals_at_risk,events,event_rate,cumulative_mean,standard_error\n",
  )
  for step in steps {
    builder.write_string(step.time.to_string())
    builder.write_char(',')
    builder.write_string(step.subjects_at_risk.to_string())
    builder.write_char(',')
    builder.write_string(step.intervals_at_risk.to_string())
    builder.write_char(',')
    builder.write_string(step.events.to_string())
    builder.write_char(',')
    builder.write_string(step.event_rate.to_string())
    builder.write_char(',')
    builder.write_string(step.cumulative_mean.to_string())
    builder.write_char(',')
    builder.write_string(step.standard_error.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate interval continuity and event coding.
pub fn recurrent_is_valid(data : RecurrentData) -> Bool {
  if data.records.length() == 0 || data.subject_count <= 0 {
    return false
  }
  for record in data.records {
    if record.subject < 0 ||
      record.start < 0.0 ||
      record.stop <= record.start ||
      (record.event != 0 && record.event != 1) {
      return false
    }
  }
  true
}
