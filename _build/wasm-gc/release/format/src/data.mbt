///|
/// A single observation in survival analysis.
pub(all) struct Observation {
  /// The time of the event or censoring.
  time : Double
  /// 1 if the event occurred, 0 if censored.
  event : Int
} derive(Eq, Debug, Compare)

///|
/// SurvivalData represents a collection of observations.
pub(all) struct SurvivalData {
  observations : Array[Observation]
} derive(Eq, Debug)

///|
/// Creates a new SurvivalData instance from times and events.
///
/// # Panics
/// Panics if the lengths of `times` and `events` do not match.
pub fn SurvivalData::new(
  times : Array[Double],
  events : Array[Int],
) -> SurvivalData {
  if times.length() != events.length() {
    panic()
  }
  let observations = []
  for i in 0..<times.length() {
    let evt = events[i]
    if evt != 0 && evt != 1 {
      panic() // Events must be 0 (censored) or 1 (event)
    }
    observations.push({ time: times[i], event: evt })
  }
  // Sort by time ascending, then by event descending (events happen before censoring at same time)
  observations.sort_by(fn(a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      b.event.compare(a.event)
    }
  })
  { observations, }
}

///|
/// Returns the total number of observations.
pub fn SurvivalData::length(self : SurvivalData) -> Int {
  self.observations.length()
}

///|
/// Returns the observation at the specified index.
pub fn SurvivalData::get(self : SurvivalData, index : Int) -> Observation {
  self.observations[index]
}
