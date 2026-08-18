///|
/// A competing-risk observation. Cause zero is right censoring; positive
/// causes identify mutually exclusive event types.
pub(all) struct CompetingObservation {
  time : Double
  cause : Int
} derive(Eq, Debug)

///|
/// Validated competing-risk input.
pub(all) struct CompetingRiskData {
  observations : Array[CompetingObservation]
  n_causes : Int
} derive(Eq, Debug)

///|
/// Construct competing-risk data from parallel time and cause arrays.
pub fn CompetingRiskData::new(
  times : Array[Double],
  causes : Array[Int],
) -> CompetingRiskData {
  if times.length() != causes.length() {
    panic()
  }
  let observations = []
  let mut n_causes = 0
  for i in 0..<times.length() {
    if times[i].is_nan() || times[i].is_inf() || times[i] < 0.0 {
      panic()
    }
    if causes[i] < 0 {
      panic()
    }
    if causes[i] > n_causes {
      n_causes = causes[i]
    }
    observations.push({ time: times[i], cause: causes[i] })
  }
  observations.sort_by(fn(a, b) {
    if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      b.cause.compare(a.cause)
    }
  })
  { observations, n_causes }
}

///|
/// Number of observations.
pub fn CompetingRiskData::length(self : CompetingRiskData) -> Int {
  self.observations.length()
}

///|
/// Return a sorted competing-risk observation.
pub fn CompetingRiskData::get(
  self : CompetingRiskData,
  index : Int,
) -> CompetingObservation {
  self.observations[index]
}

///|
/// Aalen-Johansen step for one or more competing causes.
pub(all) struct CIFStep {
  time : Double
  n_at_risk : Int
  n_censored : Int
  n_events : Int
  cause_events : Array[Int]
  cause_hazards : Array[Double]
  survival : Double
  cumulative_incidence : Array[Double]
} derive(Eq, Debug)

///|
/// Cumulative incidence function result.
pub(all) struct CIFResult {
  steps : Array[CIFStep]
  n_causes : Int
  sample_size : Int
} derive(Eq, Debug)

///|
/// Estimate cumulative incidence for every cause with the Aalen-Johansen
/// recursion.  The implementation handles tied causes at the same time as a
/// single risk-set transition.
pub fn cumulative_incidence(data : CompetingRiskData) -> CIFResult {
  let steps = []
  let mut index = 0
  let mut at_risk = data.length()
  let mut survival = 1.0
  let incidence = zeros(data.n_causes)
  while index < data.length() {
    let time = data.get(index).time
    let cause_events = []
    for _ in 0..<data.n_causes {
      cause_events.push(0)
    }
    let mut censored = 0
    while index < data.length() && data.get(index).time == time {
      let cause = data.get(index).cause
      if cause == 0 {
        censored += 1
      } else {
        cause_events[cause - 1] += 1
      }
      index += 1
    }
    let hazards = []
    let mut total_hazard = 0.0
    let mut event_count = 0
    for cause_event in cause_events {
      let hazard = cause_event.to_double() / at_risk.to_double()
      hazards.push(hazard)
      total_hazard += hazard
      event_count += cause_event
    }
    for i in 0..<data.n_causes {
      incidence[i] += survival * hazards[i]
    }
    survival *= 1.0 - total_hazard
    let incidence_snapshot = []
    for value in incidence {
      incidence_snapshot.push(value)
    }
    steps.push({
      time,
      n_at_risk: at_risk,
      n_censored: censored,
      n_events: event_count,
      cause_events,
      cause_hazards: hazards,
      survival,
      cumulative_incidence: incidence_snapshot,
    })
    at_risk -= censored + event_count
  }
  { steps, n_causes: data.n_causes, sample_size: data.length() }
}

///|
/// Predict cumulative incidence for one cause at one time.
pub fn CIFResult::predict(
  self : CIFResult,
  cause : Int,
  time : Double,
) -> Double {
  if cause <= 0 || cause > self.n_causes {
    return 0.0
  }
  let mut value = 0.0
  for step in self.steps {
    if step.time <= time {
      value = step.cumulative_incidence[cause - 1]
    } else {
      break
    }
  }
  value
}

///|
/// Predict the overall event-free survival under competing risks.
pub fn CIFResult::survival_at(self : CIFResult, time : Double) -> Double {
  let mut value = 1.0
  for step in self.steps {
    if step.time <= time {
      value = step.survival
    } else {
      break
    }
  }
  value
}

///|
/// Return every cause's incidence at a time.
pub fn CIFResult::predict_all(self : CIFResult, time : Double) -> Array[Double] {
  let result = zeros(self.n_causes)
  for step in self.steps {
    if step.time <= time {
      for i in 0..<self.n_causes {
        result[i] = step.cumulative_incidence[i]
      }
    } else {
      break
    }
  }
  result
}

///|
/// Cause-specific hazard increment at a time.
pub fn CIFResult::cause_hazard_at(
  self : CIFResult,
  cause : Int,
  time : Double,
) -> Double {
  if cause <= 0 || cause > self.n_causes {
    return 0.0
  }
  let mut value = 0.0
  for step in self.steps {
    if step.time <= time {
      value = step.cause_hazards[cause - 1]
    } else {
      break
    }
  }
  value
}

///|
/// Export cumulative incidence curves in long format.
pub fn CIFResult::to_csv(self : CIFResult) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,survival,cause,cause_events,cause_hazard,cumulative_incidence\n",
  )
  for step in self.steps {
    for cause in 0..<self.n_causes {
      builder.write_string(
        "\{step.time},\{step.survival},\{cause + 1},\{step.cause_events[cause]},\{step.cause_hazards[cause]},\{step.cumulative_incidence[cause]}\n",
      )
    }
  }
  builder.to_string()
}

///|
/// Simple cause-specific event counts at the end of follow-up.
pub fn competing_event_counts(data : CompetingRiskData) -> Array[Int] {
  let counts = []
  for _ in 0..<data.n_causes {
    counts.push(0)
  }
  for observation in data.observations {
    if observation.cause > 0 {
      counts[observation.cause - 1] += 1
    }
  }
  counts
}
