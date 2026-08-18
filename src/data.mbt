///|
/// A single right-censored survival observation.
///
/// `event` is one when the event of interest was observed and zero when the
/// observation ended in right censoring.  Time is measured in the unit chosen
/// by the caller (days, hours, cycles, or releases).
pub(all) struct Observation {
  /// The event or censoring time.
  time : Double
  /// One for an event and zero for right censoring.
  event : Int
} derive(Eq, Debug, Compare)

///|
/// Structured validation failures returned by the non-panicking constructors.
pub(all) enum DataError {
  LengthMismatch
  InvalidEvent(Int)
  InvalidTime
  NegativeTime
  EmptyCovariates
  CovariateLengthMismatch
} derive(Eq, Debug)

///|
/// SurvivalData owns a sorted collection of observations.
pub(all) struct SurvivalData {
  observations : Array[Observation]
} derive(Eq, Debug)

///|
/// Create a validated dataset without panicking.
pub fn SurvivalData::try_new(
  times : Array[Double],
  events : Array[Int],
) -> Result[SurvivalData, DataError] {
  if times.length() != events.length() {
    return Err(LengthMismatch)
  }
  let observations = []
  for i in 0..<times.length() {
    let time = times[i]
    let event = events[i]
    if time.is_nan() || time.is_inf() {
      return Err(InvalidTime)
    }
    if time < 0.0 {
      return Err(NegativeTime)
    }
    if event != 0 && event != 1 {
      return Err(InvalidEvent(event))
    }
    observations.push({ time, event })
  }
  observations.sort_by(observation_compare)
  Ok({ observations, })
}

///|
/// Creates a SurvivalData instance from parallel arrays.
///
/// # Panics
/// Panics when array lengths, times, or event indicators are invalid. Use
/// `SurvivalData::try_new` when input comes from a file or user interface.
pub fn SurvivalData::new(
  times : Array[Double],
  events : Array[Int],
) -> SurvivalData {
  match SurvivalData::try_new(times, events) {
    Ok(data) => data
    Err(_) => panic()
  }
}

///|
/// Create a dataset from already materialized observations.
pub fn SurvivalData::from_observations(
  observations : Array[Observation],
) -> SurvivalData {
  let times = []
  let events = []
  for observation in observations {
    times.push(observation.time)
    events.push(observation.event)
  }
  SurvivalData::new(times, events)
}

///|
/// Return the total number of observations.
pub fn SurvivalData::length(self : SurvivalData) -> Int {
  self.observations.length()
}

///|
/// Return an observation by sorted position.
pub fn SurvivalData::get(self : SurvivalData, index : Int) -> Observation {
  self.observations[index]
}

///|
/// Return a detached copy of the sorted observations.
pub fn SurvivalData::to_array(self : SurvivalData) -> Array[Observation] {
  let result = []
  for observation in self.observations {
    result.push(observation)
  }
  result
}

///|
/// Return event times only.
pub fn SurvivalData::event_times(self : SurvivalData) -> Array[Double] {
  let result = []
  for observation in self.observations {
    if observation.event == 1 {
      result.push(observation.time)
    }
  }
  result
}

///|
/// Return censoring times only.
pub fn SurvivalData::censoring_times(self : SurvivalData) -> Array[Double] {
  let result = []
  for observation in self.observations {
    if observation.event == 0 {
      result.push(observation.time)
    }
  }
  result
}

///|
/// Count observed events.
pub fn SurvivalData::event_count(self : SurvivalData) -> Int {
  let mut count = 0
  for observation in self.observations {
    if observation.event == 1 {
      count += 1
    }
  }
  count
}

///|
/// Count right-censored observations.
pub fn SurvivalData::censored_count(self : SurvivalData) -> Int {
  self.length() - self.event_count()
}

///|
/// Small descriptive summary useful for dashboards and data quality checks.
pub(all) struct SurvivalSummary {
  n : Int
  events : Int
  censored : Int
  min_time : Double
  max_time : Double
  event_fraction : Double
} derive(Eq, Debug)

///|
/// Compute basic data quality and sample-size information.
pub fn SurvivalData::summary(self : SurvivalData) -> SurvivalSummary {
  if self.length() == 0 {
    return {
      n: 0,
      events: 0,
      censored: 0,
      min_time: 0.0,
      max_time: 0.0,
      event_fraction: 0.0,
    }
  }
  let mut minimum = self.observations[0].time
  let mut maximum = minimum
  for observation in self.observations {
    if observation.time < minimum {
      minimum = observation.time
    }
    if observation.time > maximum {
      maximum = observation.time
    }
  }
  let events = self.event_count()
  {
    n: self.length(),
    events,
    censored: self.length() - events,
    min_time: minimum,
    max_time: maximum,
    event_fraction: events.to_double() / self.length().to_double(),
  }
}

///|
/// Return a time window, retaining observations at both boundaries.
pub fn SurvivalData::between(
  self : SurvivalData,
  start : Double,
  end : Double,
) -> SurvivalData {
  let times = []
  let events = []
  for observation in self.observations {
    if observation.time >= start && observation.time <= end {
      times.push(observation.time)
      events.push(observation.event)
    }
  }
  SurvivalData::new(times, events)
}

///|
/// Collapse a dataset into unique time groups and count events and censors.
pub(all) struct EventTableRow {
  time : Double
  n_at_risk : Int
  n_events : Int
  n_censored : Int
} derive(Eq, Debug)

///|
/// Build the event table used by all non-parametric estimators.
pub fn event_table(data : SurvivalData) -> Array[EventTableRow] {
  let rows = []
  if data.length() == 0 {
    return rows
  }
  let mut index = 0
  let mut at_risk = data.length()
  while index < data.length() {
    let time = data.get(index).time
    let mut events = 0
    let mut censored = 0
    let mut cursor = index
    while cursor < data.length() && data.get(cursor).time == time {
      if data.get(cursor).event == 1 {
        events += 1
      } else {
        censored += 1
      }
      cursor += 1
    }
    rows.push({
      time,
      n_at_risk: at_risk,
      n_events: events,
      n_censored: censored,
    })
    at_risk -= events + censored
    index = cursor
  }
  rows
}

///|
/// Compare observations by time and place events before censoring at a tie.
fn observation_compare(a : Observation, b : Observation) -> Int {
  if a.time < b.time {
    -1
  } else if a.time > b.time {
    1
  } else {
    b.event.compare(a.event)
  }
}

///|
/// A weighted observation for inverse-probability and survey-weighted
/// survival analysis.
pub(all) struct WeightedObservation {
  time : Double
  event : Int
  weight : Double
} derive(Eq, Debug)

///|
/// A validated weighted dataset.
pub(all) struct WeightedSurvivalData {
  observations : Array[WeightedObservation]
} derive(Eq, Debug)

///|
/// Create a weighted dataset; all weights must be finite and positive.
pub fn WeightedSurvivalData::try_new(
  times : Array[Double],
  events : Array[Int],
  weights : Array[Double],
) -> Result[WeightedSurvivalData, DataError] {
  if times.length() != events.length() || times.length() != weights.length() {
    return Err(LengthMismatch)
  }
  let observations = []
  for i in 0..<times.length() {
    if times[i].is_nan() || times[i].is_inf() {
      return Err(InvalidTime)
    }
    if times[i] < 0.0 {
      return Err(NegativeTime)
    }
    if events[i] != 0 && events[i] != 1 {
      return Err(InvalidEvent(events[i]))
    }
    if weights[i].is_nan() || weights[i].is_inf() || weights[i] <= 0.0 {
      return Err(InvalidTime)
    }
    observations.push({ time: times[i], event: events[i], weight: weights[i] })
  }
  observations.sort_by(fn(a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      b.event.compare(a.event)
    }
  })
  Ok({ observations, })
}

///|
/// Panicking convenience constructor for weighted data.
pub fn WeightedSurvivalData::new(
  times : Array[Double],
  events : Array[Int],
  weights : Array[Double],
) -> WeightedSurvivalData {
  match WeightedSurvivalData::try_new(times, events, weights) {
    Ok(data) => data
    Err(_) => panic()
  }
}

///|
/// Number of weighted records.
pub fn WeightedSurvivalData::length(self : WeightedSurvivalData) -> Int {
  self.observations.length()
}

///|
/// Access a weighted record.
pub fn WeightedSurvivalData::get(
  self : WeightedSurvivalData,
  index : Int,
) -> WeightedObservation {
  self.observations[index]
}
