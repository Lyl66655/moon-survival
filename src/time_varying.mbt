///|
/// One counting-process interval for time-varying covariates.
pub(all) struct TimeVaryingObservation {
  start : Double
  stop : Double
  event : Int
  covariates : Array[Double]
  subject : Int
} derive(Eq, Debug)

///|
/// Validated counting-process input.
pub(all) struct TimeVaryingData {
  observations : Array[TimeVaryingObservation]
  n_covariates : Int
  n_subjects : Int
} derive(Eq, Debug)

///|
/// Construct counting-process data.  Every interval must satisfy
/// `0 <= start < stop`, and at most one event indicator is allowed.
pub fn TimeVaryingData::new(
  starts : Array[Double],
  stops : Array[Double],
  events : Array[Int],
  covs : Array[Array[Double]],
  subjects : Array[Int],
) -> TimeVaryingData {
  if starts.length() != stops.length() ||
    starts.length() != events.length() ||
    starts.length() != covs.length() ||
    starts.length() != subjects.length() {
    panic()
  }
  let n_covariates = if covs.length() == 0 { 0 } else { covs[0].length() }
  let observations = []
  let mut n_subjects = 0
  for i in 0..<starts.length() {
    if starts[i].is_nan() ||
      stops[i].is_nan() ||
      starts[i] < 0.0 ||
      stops[i] <= starts[i] {
      panic()
    }
    if events[i] != 0 && events[i] != 1 {
      panic()
    }
    if subjects[i] < 0 || covs[i].length() != n_covariates {
      panic()
    }
    if subjects[i] + 1 > n_subjects {
      n_subjects = subjects[i] + 1
    }
    observations.push({
      start: starts[i],
      stop: stops[i],
      event: events[i],
      covariates: covs[i],
      subject: subjects[i],
    })
  }
  observations.sort_by(fn(a, b) {
    if a.stop < b.stop {
      -1
    } else if a.stop > b.stop {
      1
    } else {
      a.subject.compare(b.subject)
    }
  })
  { observations, n_covariates, n_subjects }
}

///|
/// Number of counting-process intervals.
pub fn TimeVaryingData::length(self : TimeVaryingData) -> Int {
  self.observations.length()
}

///|
/// Access an interval.
pub fn TimeVaryingData::get(
  self : TimeVaryingData,
  index : Int,
) -> TimeVaryingObservation {
  self.observations[index]
}

///|
/// Return all intervals for one subject.
pub fn TimeVaryingData::subject_intervals(
  self : TimeVaryingData,
  subject : Int,
) -> Array[TimeVaryingObservation] {
  let result = []
  for observation in self.observations {
    if observation.subject == subject {
      result.push(observation)
    }
  }
  result
}

///|
/// Check that intervals for a subject are non-overlapping and ordered.
pub fn TimeVaryingData::is_valid_subject(
  self : TimeVaryingData,
  subject : Int,
) -> Bool {
  let intervals = self.subject_intervals(subject)
  if intervals.length() < 2 {
    true
  } else {
    let mut valid = true
    for i in 1..<intervals.length() {
      if intervals[i].start < intervals[i - 1].stop {
        valid = false
      }
    }
    valid
  }
}

///|
/// Count-process event table row.
pub(all) struct TimeVaryingRiskRow {
  stop : Double
  n_at_risk : Int
  n_events : Int
} derive(Eq, Debug)

///|
/// Build risk-set counts using interval membership at each event time.
pub fn time_varying_risk_table(
  data : TimeVaryingData,
) -> Array[TimeVaryingRiskRow] {
  let event_times = []
  for observation in data.observations {
    if observation.event == 1 {
      event_times.push(observation.stop)
    }
  }
  event_times.sort()
  let unique = []
  for time in event_times {
    if unique.length() == 0 || unique[unique.length() - 1] != time {
      unique.push(time)
    }
  }
  let result = []
  for time in unique {
    let mut risk = 0
    let mut events = 0
    for observation in data.observations {
      if observation.start < time && observation.stop >= time {
        risk += 1
      }
      if observation.event == 1 && observation.stop == time {
        events += 1
      }
    }
    result.push({ stop: time, n_at_risk: risk, n_events: events })
  }
  result
}

///|
/// Convert the final interval per subject to ordinary right-censored data.
pub fn TimeVaryingData::to_survival_data(
  self : TimeVaryingData,
) -> SurvivalData {
  let times = []
  let events = []
  for subject in 0..<self.n_subjects {
    let intervals = self.subject_intervals(subject)
    if intervals.length() > 0 {
      let last = intervals[intervals.length() - 1]
      times.push(last.stop)
      events.push(last.event)
    }
  }
  SurvivalData::new(times, events)
}

///|
/// A landmark snapshot of time-varying covariates.
pub(all) struct LandmarkRow {
  subject : Int
  time : Double
  covariates : Array[Double]
  event_after_landmark : Int
} derive(Eq, Debug)

///|
/// Extract the active covariate row for every subject at a landmark time.
pub fn landmark_snapshot(
  data : TimeVaryingData,
  landmark : Double,
) -> Array[LandmarkRow] {
  let result = []
  for subject in 0..<data.n_subjects {
    let intervals = data.subject_intervals(subject)
    let mut found : TimeVaryingObservation? = None
    for observation in intervals {
      if observation.start <= landmark && landmark < observation.stop {
        found = Some(observation)
      }
    }
    match found {
      Some(observation) => {
        let covariates = []
        for value in observation.covariates {
          covariates.push(value)
        }
        result.push({
          subject,
          time: observation.stop,
          covariates,
          event_after_landmark: observation.event,
        })
      }
      None => ()
    }
  }
  result
}

///|
/// A start-stop record formatted for interchange with common survival tools.
pub fn TimeVaryingData::to_csv(self : TimeVaryingData) -> String {
  let builder = StringBuilder::new()
  builder.write_string("subject,start,stop,event,covariates\n")
  for observation in self.observations {
    builder.write_string(
      "\{observation.subject},\{observation.start},\{observation.stop},\{observation.event},",
    )
    for i in 0..<observation.covariates.length() {
      if i > 0 {
        builder.write_string(";")
      }
      builder.write_string("\{observation.covariates[i]}")
    }
    builder.write_string("\n")
  }
  builder.to_string()
}
