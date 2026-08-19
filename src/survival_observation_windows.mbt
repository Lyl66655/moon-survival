///|
/// A subject observation window.
pub(all) struct SurvivalObservationWindow {
  start : Double
  end : Double
  at_risk : Int
  events : Int
  censored : Int
  exposure : Double
  window_row : Bool
} derive(Eq, Debug)

///|
/// Build fixed-width observation windows.
pub fn survival_observation_windows(
  data : SurvivalData,
  width : Double,
  horizon : Double,
) -> Array[SurvivalObservationWindow] {
  let result : Array[SurvivalObservationWindow] = []
  if width <= 0.0 || horizon <= 0.0 {
    return result
  }
  let mut start = 0.0
  while start < horizon {
    let end = if start + width > horizon { horizon } else { start + width }
    let mut at_risk = 0
    let mut events = 0
    let mut censored = 0
    let mut exposure = 0.0
    for observation in data.to_array() {
      if observation.time >= start {
        at_risk += 1
      }
      if observation.time >= start && observation.time < end {
        if observation.event == 1 {
          events += 1
        } else {
          censored += 1
        }
      }
      if observation.time > start {
        exposure += if observation.time < end {
          observation.time - start
        } else {
          end - start
        }
      }
    }
    result.push({
      start,
      end,
      at_risk,
      events,
      censored,
      exposure,
      window_row: true,
    })
    start = end
  }
  result
}

///|
/// Convert windows to event rates.
pub fn survival_window_event_rates(
  windows : Array[SurvivalObservationWindow],
) -> Array[Double] {
  let result : Array[Double] = []
  for window in windows {
    result.push(
      survival_ext_safe_ratio(window.events.to_double(), window.exposure),
    )
  }
  result
}

///|
/// Compute the cumulative exposure.
pub fn survival_window_total_exposure(
  windows : Array[SurvivalObservationWindow],
) -> Double {
  let mut total = 0.0
  for window in windows {
    total += window.exposure
  }
  total
}

///|
/// Compute the total observed events.
pub fn survival_window_total_events(
  windows : Array[SurvivalObservationWindow],
) -> Int {
  let mut total = 0
  for window in windows {
    total += window.events
  }
  total
}

///|
/// Return the window containing a time.
pub fn survival_window_at(
  windows : Array[SurvivalObservationWindow],
  time : Double,
) -> SurvivalObservationWindow? {
  for window in windows {
    if time >= window.start && time < window.end {
      return Some(window)
    }
  }
  None
}

///|
/// Compute window-level Poisson confidence intervals.
pub fn survival_window_rate_intervals(
  windows : Array[SurvivalObservationWindow],
) -> Array[SurvivalEstimateInterval] {
  let result : Array[SurvivalEstimateInterval] = []
  for window in windows {
    let estimate = survival_ext_safe_ratio(
      window.events.to_double(),
      window.exposure,
    )
    let standard_error = if window.events == 0 {
      0.0
    } else {
      double_sqrt(window.events.to_double()) / clamp_positive(window.exposure)
    }
    result.push(survival_estimate_interval(estimate, standard_error, 0.95))
  }
  result
}

///|
/// Export observation windows.
pub fn survival_observation_windows_to_csv(
  windows : Array[SurvivalObservationWindow],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("start,end,at_risk,events,censored,exposure\n")
  for window in windows {
    builder.write_string(window.start.to_string())
    builder.write_char(',')
    builder.write_string(window.end.to_string())
    builder.write_char(',')
    builder.write_string(window.at_risk.to_string())
    builder.write_char(',')
    builder.write_string(window.events.to_string())
    builder.write_char(',')
    builder.write_string(window.censored.to_string())
    builder.write_char(',')
    builder.write_string(window.exposure.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate non-overlapping windows.
pub fn survival_observation_windows_is_valid(
  windows : Array[SurvivalObservationWindow],
) -> Bool {
  for i in 0..<windows.length() {
    if windows[i].start < 0.0 ||
      windows[i].end <= windows[i].start ||
      windows[i].events < 0 ||
      windows[i].censored < 0 ||
      (i > 0 && windows[i].start < windows[i - 1].end) {
      return false
    }
  }
  true
}
