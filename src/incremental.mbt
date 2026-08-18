///|
/// Append-only collector for applications that receive observations in
/// batches and fit a curve at checkpoints.
pub(all) struct SurvivalAccumulator {
  times : Array[Double]
  events : Array[Int]
} derive(Eq, Debug)

///|
/// Create an empty accumulator.
pub fn SurvivalAccumulator::new() -> SurvivalAccumulator {
  { times: [], events: [] }
}

///|
/// Add one validated observation.
pub fn SurvivalAccumulator::add(
  self : SurvivalAccumulator,
  time : Double,
  event : Int,
) -> Unit {
  if time.is_nan() || time.is_inf() || time < 0.0 || (event != 0 && event != 1) {
    panic()
  }
  self.times.push(time)
  self.events.push(event)
}

///|
/// Add a batch of parallel observations.
pub fn SurvivalAccumulator::add_batch(
  self : SurvivalAccumulator,
  times : Array[Double],
  events : Array[Int],
) -> Unit {
  if times.length() != events.length() {
    panic()
  }
  for i in 0..<times.length() {
    self.add(times[i], events[i])
  }
}

///|
/// Merge another accumulator without sharing its arrays.
pub fn SurvivalAccumulator::merge(
  self : SurvivalAccumulator,
  other : SurvivalAccumulator,
) -> Unit {
  self.add_batch(other.times, other.events)
}

///|
/// Number of collected observations.
pub fn SurvivalAccumulator::length(self : SurvivalAccumulator) -> Int {
  self.times.length()
}

///|
/// Number of collected events.
pub fn SurvivalAccumulator::event_count(self : SurvivalAccumulator) -> Int {
  let mut count = 0
  for event in self.events {
    if event == 1 {
      count += 1
    }
  }
  count
}

///|
/// Snapshot the current observations as validated sorted data.
pub fn SurvivalAccumulator::snapshot(
  self : SurvivalAccumulator,
) -> SurvivalData {
  SurvivalData::new(self.times, self.events)
}

///|
/// Fit a KM curve from the current checkpoint.
pub fn SurvivalAccumulator::fit(self : SurvivalAccumulator) -> KMResult {
  kaplan_meier(self.snapshot())
}

///|
/// A streaming checkpoint with a named timestamp.
pub(all) struct SurvivalCheckpoint {
  label : String
  sample_size : Int
  events : Int
  median : Double?
  restricted_mean : Double
} derive(Eq, Debug)

///|
/// Materialize a checkpoint for operations dashboards.
pub fn SurvivalAccumulator::checkpoint(
  self : SurvivalAccumulator,
  label : String,
  horizon : Double,
) -> SurvivalCheckpoint {
  let fit = self.fit()
  {
    label,
    sample_size: self.length(),
    events: self.event_count(),
    median: fit.median(),
    restricted_mean: fit.restricted_mean(horizon),
  }
}

///|
/// Return the currently collected rows as CSV.
pub fn SurvivalAccumulator::to_csv(self : SurvivalAccumulator) -> String {
  survival_data_to_csv(self.snapshot())
}
