///|
/// A step in the Nelson-Aalen cumulative hazard curve.
pub(all) struct NAStep {
  time : Double
  n_at_risk : Int
  n_events : Int
  n_censored : Int
  cum_hazard : Double
} derive(Eq, Debug)

///|
/// Result of a Nelson-Aalen estimation.
pub(all) struct NAResult {
  steps : Array[NAStep]
} derive(Eq, Debug)

///|
/// Estimates the Nelson-Aalen cumulative hazard curve from the given data.
pub fn nelson_aalen(data : SurvivalData) -> NAResult {
  let steps = []
  let mut current_time = -1.0
  let mut n_at_risk = data.length()
  let mut current_events = 0
  let mut current_censored = 0
  let mut current_hazard = 0.0

  if n_at_risk == 0 {
    return { steps: [] }
  }

  for i in 0..<data.length() {
    let obs = data.get(i)

    if current_time < 0.0 {
      current_time = obs.time
    }

    if obs.time > current_time {
      // Finalize previous time step
      if current_events > 0 || current_censored > 0 {
        let hazard = current_events.to_double() / n_at_risk.to_double()
        current_hazard = current_hazard + hazard

        steps.push({
          time: current_time,
          n_at_risk,
          n_events: current_events,
          n_censored: current_censored,
          cum_hazard: current_hazard,
        })
      }

      n_at_risk = n_at_risk - current_events - current_censored
      current_time = obs.time
      current_events = 0
      current_censored = 0
    }

    if obs.event == 1 {
      current_events += 1
    } else {
      current_censored += 1
    }
  }

  // Finalize the last time step
  if current_events > 0 || current_censored > 0 {
    let hazard = current_events.to_double() / n_at_risk.to_double()
    current_hazard = current_hazard + hazard

    steps.push({
      time: current_time,
      n_at_risk,
      n_events: current_events,
      n_censored: current_censored,
      cum_hazard: current_hazard,
    })
  }

  { steps, }
}

///|
/// Predict the cumulative hazard at a given time `t`.
pub fn NAResult::predict(self : NAResult, t : Double) -> Double {
  let mut h = 0.0
  for i in 0..<self.steps.length() {
    if self.steps[i].time <= t {
      h = self.steps[i].cum_hazard
    } else {
      break
    }
  }
  h
}
