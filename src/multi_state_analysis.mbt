///|
/// A subject transition in a continuous-time multi-state process.
pub(all) struct MultiStateTransition {
  subject : Int
  from_state : Int
  to_state : Int
  time : Double
  multistate_row : Bool
} derive(Eq, Debug)

///|
/// Validation errors for multi-state transition data.
pub(all) enum MultiStateError {
  EmptyData
  LengthMismatch
  InvalidSubject
  InvalidState
  InvalidTime
  SelfTransition
  UnorderedSubject
  InvalidInitialState
} derive(Eq, Debug)

///|
/// Validated transition history with one initial state per subject.
pub(all) struct MultiStateData {
  transitions : Array[MultiStateTransition]
  initial_states : Array[Int]
  state_count : Int
  subject_count : Int
} derive(Eq, Debug)

///|
/// Aalen-Johansen state probability at one event time.
pub(all) struct StateProbabilityStep {
  time : Double
  state : Int
  probability : Double
  standard_error : Double
  state_step : Bool
} derive(Eq, Debug)

///|
/// A transition-specific cumulative hazard estimate.
pub(all) struct TransitionHazardStep {
  time : Double
  from_state : Int
  to_state : Int
  at_risk : Int
  transitions : Int
  hazard_increment : Double
  cumulative_hazard : Double
  hazard_step : Bool
} derive(Eq, Debug)

///|
/// Expected occupancy and sojourn summaries.
pub(all) struct StateOccupancySummary {
  state : Int
  restricted_occupancy : Double
  occupancy_probability : Double
  mean_first_entry : Double
  entry_count : Int
} derive(Eq, Debug)

///|
/// Validate and build a multi-state transition history.
pub fn MultiStateData::try_new(
  subjects : Array[Int],
  from_states : Array[Int],
  to_states : Array[Int],
  times : Array[Double],
  initial_states : Array[Int],
  state_count : Int,
) -> Result[MultiStateData, MultiStateError] {
  if subjects.length() == 0 {
    return Err(EmptyData)
  }
  if state_count <= 0 {
    return Err(InvalidState)
  }
  if subjects.length() != from_states.length() ||
    subjects.length() != to_states.length() ||
    subjects.length() != times.length() {
    return Err(LengthMismatch)
  }
  let mut maximum_subject = -1
  for subject in subjects {
    if subject < 0 {
      return Err(InvalidSubject)
    }
    if subject > maximum_subject {
      maximum_subject = subject
    }
  }
  if initial_states.length() <= maximum_subject {
    return Err(InvalidInitialState)
  }
  for state in initial_states {
    if state < 0 || state >= state_count {
      return Err(InvalidInitialState)
    }
  }
  let transitions : Array[MultiStateTransition] = []
  for i in 0..<subjects.length() {
    if from_states[i] < 0 ||
      from_states[i] >= state_count ||
      to_states[i] < 0 ||
      to_states[i] >= state_count {
      return Err(InvalidState)
    }
    if from_states[i] == to_states[i] {
      return Err(SelfTransition)
    }
    if times[i].is_nan() || times[i].is_inf() || times[i] < 0.0 {
      return Err(InvalidTime)
    }
    transitions.push({
      subject: subjects[i],
      from_state: from_states[i],
      to_state: to_states[i],
      time: times[i],
      multistate_row: true,
    })
  }
  transitions.sort_by(fn(a, b) {
    if a.subject < b.subject {
      -1
    } else if a.subject > b.subject {
      1
    } else if a.time < b.time {
      -1
    } else if a.time > b.time {
      1
    } else {
      a.to_state.compare(b.to_state)
    }
  })
  let mut previous_subject = -1
  let mut previous_time = 0.0
  for transition in transitions {
    if transition.subject != previous_subject {
      previous_subject = transition.subject
      previous_time = transition.time
    } else if transition.time < previous_time {
      return Err(UnorderedSubject)
    } else {
      previous_time = transition.time
    }
  }
  Ok({
    transitions,
    initial_states,
    state_count,
    subject_count: maximum_subject + 1,
  })
}

///|
/// Panic-on-invalid convenience constructor.
pub fn MultiStateData::new(
  subjects : Array[Int],
  from_states : Array[Int],
  to_states : Array[Int],
  times : Array[Double],
  initial_states : Array[Int],
  state_count : Int,
) -> MultiStateData {
  match
    MultiStateData::try_new(
      subjects, from_states, to_states, times, initial_states, state_count,
    ) {
    Ok(data) => data
    Err(_) => panic()
  }
}

///|
/// Number of transitions.
pub fn MultiStateData::length(self : MultiStateData) -> Int {
  self.transitions.length()
}

///|
/// Number of states.
pub fn MultiStateData::state_count(self : MultiStateData) -> Int {
  self.state_count
}

///|
/// Number of subject trajectories.
pub fn MultiStateData::subject_count(self : MultiStateData) -> Int {
  self.subject_count
}

///|
/// Access a transition.
pub fn MultiStateData::get(
  self : MultiStateData,
  index : Int,
) -> MultiStateTransition {
  self.transitions[index]
}

///|
/// Return the event-time grid.
pub fn multistate_event_times(data : MultiStateData) -> Array[Double] {
  let result : Array[Double] = []
  for transition in data.transitions {
    if !result.contains(transition.time) {
      result.push(transition.time)
    }
  }
  result.sort()
  result
}

///|
/// Return the state occupied just before a subject transition.
fn multistate_subject_state(
  data : MultiStateData,
  subject : Int,
  time : Double,
) -> Int {
  let mut state = data.initial_states[subject]
  for transition in data.transitions {
    if transition.subject == subject && transition.time < time {
      state = transition.to_state
    }
  }
  state
}

///|
/// Count transitions from one state to another at a time.
fn multistate_transition_count(
  data : MultiStateData,
  from_state : Int,
  to_state : Int,
  time : Double,
) -> Int {
  let mut count = 0
  for transition in data.transitions {
    if transition.time == time &&
      transition.from_state == from_state &&
      transition.to_state == to_state {
      count += 1
    }
  }
  count
}

///|
/// Count the number of subjects at risk in a state just before a time.
fn multistate_at_risk(data : MultiStateData, state : Int, time : Double) -> Int {
  let mut count = 0
  for subject in 0..<data.subject_count {
    if multistate_subject_state(data, subject, time) == state {
      count += 1
    }
  }
  count
}

///|
/// Estimate transition-specific cumulative hazards.
pub fn multistate_cumulative_hazards(
  data : MultiStateData,
) -> Array[TransitionHazardStep] {
  let result : Array[TransitionHazardStep] = []
  let times = multistate_event_times(data)
  for time in times {
    for from_state in 0..<data.state_count {
      let at_risk = multistate_at_risk(data, from_state, time)
      for to_state in 0..<data.state_count {
        if from_state == to_state {
          continue
        }
        let transitions = multistate_transition_count(
          data, from_state, to_state, time,
        )
        let increment = if at_risk == 0 {
          0.0
        } else {
          transitions.to_double() / at_risk.to_double()
        }
        let mut cumulative = increment
        for previous in result {
          if previous.from_state == from_state && previous.to_state == to_state {
            cumulative += previous.hazard_increment
          }
        }
        result.push({
          time,
          from_state,
          to_state,
          at_risk,
          transitions,
          hazard_increment: increment,
          cumulative_hazard: cumulative,
          hazard_step: true,
        })
      }
    }
  }
  result
}

///|
/// Estimate Aalen-Johansen state probabilities from each subject's initial state.
pub fn multistate_state_probabilities(
  data : MultiStateData,
) -> Array[StateProbabilityStep] {
  let result : Array[StateProbabilityStep] = []
  let times = multistate_event_times(data)
  let mut probabilities : Array[Double] = Array::make(data.state_count, 0.0)
  for state in data.initial_states {
    probabilities[state] += 1.0 / data.subject_count.to_double()
  }
  let mut previous_time = 0.0
  for time in times {
    let next = []
    for probability in probabilities {
      next.push(probability)
    }
    for from_state in 0..<data.state_count {
      let at_risk = multistate_at_risk(data, from_state, time).to_double()
      if at_risk <= 0.0 {
        continue
      }
      let mut outgoing = 0.0
      for to_state in 0..<data.state_count {
        if from_state != to_state {
          outgoing += multistate_transition_count(
            data, from_state, to_state, time,
          ).to_double()
        }
      }
      for to_state in 0..<data.state_count {
        if from_state == to_state {
          continue
        }
        let count = multistate_transition_count(
          data, from_state, to_state, time,
        ).to_double()
        next[to_state] += probabilities[from_state] * count / at_risk
      }
      next[from_state] -= probabilities[from_state] * outgoing / at_risk
    }
    probabilities = next
    let elapsed = time - previous_time
    let standard_error = double_sqrt(
      elapsed.abs() / data.subject_count.to_double(),
    )
    for state in 0..<data.state_count {
      result.push({
        time,
        state,
        probability: clamp_probability(probabilities[state]),
        standard_error,
        state_step: true,
      })
    }
    previous_time = time
  }
  result
}

///|
/// Look up a state probability at a horizon.
pub fn multistate_probability_at(
  steps : Array[StateProbabilityStep],
  state : Int,
  horizon : Double,
) -> Double {
  let mut probability = 0.0
  for step in steps {
    if step.state == state && step.time <= horizon {
      probability = step.probability
    }
  }
  probability
}

///|
/// Integrate state occupancy from a probability curve.
pub fn multistate_restricted_occupancy(
  steps : Array[StateProbabilityStep],
  state : Int,
  horizon : Double,
) -> Double {
  if horizon <= 0.0 {
    return 0.0
  }
  let mut area = 0.0
  let mut previous = 0.0
  let mut probability = 0.0
  for step in steps {
    if step.state != state {
      continue
    }
    if step.time >= horizon {
      break
    }
    area += (step.time - previous) * probability
    previous = step.time
    probability = step.probability
  }
  area + (horizon - previous) * probability
}

///|
/// Count entries into a state and estimate the average entry time.
pub fn multistate_occupancy_summary(
  data : MultiStateData,
  steps : Array[StateProbabilityStep],
  horizon : Double,
) -> Array[StateOccupancySummary] {
  let result : Array[StateOccupancySummary] = []
  for state in 0..<data.state_count {
    let mut entries = 0
    let mut time_total = 0.0
    for transition in data.transitions {
      if transition.to_state == state && transition.time <= horizon {
        entries += 1
        time_total += transition.time
      }
    }
    let mean_entry = if entries == 0 {
      0.0
    } else {
      time_total / entries.to_double()
    }
    let occupancy = multistate_restricted_occupancy(steps, state, horizon)
    result.push({
      state,
      restricted_occupancy: occupancy,
      occupancy_probability: if horizon <= 0.0 {
        0.0
      } else {
        occupancy / horizon
      },
      mean_first_entry: mean_entry,
      entry_count: entries,
    })
  }
  result
}

///|
/// Return the most likely state sequence for one subject.
pub fn multistate_subject_path(
  data : MultiStateData,
  subject : Int,
) -> Array[MultiStateTransition] {
  let result : Array[MultiStateTransition] = []
  if subject < 0 || subject >= data.subject_count {
    return result
  }
  for transition in data.transitions {
    if transition.subject == subject {
      result.push(transition)
    }
  }
  result
}

///|
/// Validate that each transition starts from the subject's observed state.
pub fn multistate_consistency_report(data : MultiStateData) -> Array[String] {
  let issues : Array[String] = []
  for subject in 0..<data.subject_count {
    let mut state = data.initial_states[subject]
    for transition in data.transitions {
      if transition.subject == subject {
        if transition.from_state != state {
          issues.push("subject_" + subject.to_string() + "_state_mismatch")
        }
        state = transition.to_state
      }
    }
  }
  issues
}

///|
/// Serialize state probabilities.
pub fn state_probabilities_to_csv(
  steps : Array[StateProbabilityStep],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,state,probability,standard_error\n")
  for step in steps {
    builder.write_string(step.time.to_string())
    builder.write_char(',')
    builder.write_string(step.state.to_string())
    builder.write_char(',')
    builder.write_string(step.probability.to_string())
    builder.write_char(',')
    builder.write_string(step.standard_error.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Serialize transition hazards.
pub fn transition_hazards_to_csv(steps : Array[TransitionHazardStep]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "time,from_state,to_state,at_risk,transitions,hazard_increment,cumulative_hazard\n",
  )
  for step in steps {
    builder.write_string(step.time.to_string())
    builder.write_char(',')
    builder.write_string(step.from_state.to_string())
    builder.write_char(',')
    builder.write_string(step.to_state.to_string())
    builder.write_char(',')
    builder.write_string(step.at_risk.to_string())
    builder.write_char(',')
    builder.write_string(step.transitions.to_string())
    builder.write_char(',')
    builder.write_string(step.hazard_increment.to_string())
    builder.write_char(',')
    builder.write_string(step.cumulative_hazard.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Check basic probability and transition invariants.
pub fn multistate_is_valid(data : MultiStateData) -> Bool {
  if data.state_count <= 0 ||
    data.subject_count <= 0 ||
    data.initial_states.length() != data.subject_count {
    return false
  }
  for transition in data.transitions {
    if transition.subject < 0 ||
      transition.from_state < 0 ||
      transition.from_state >= data.state_count ||
      transition.to_state < 0 ||
      transition.to_state >= data.state_count ||
      transition.from_state == transition.to_state ||
      transition.time < 0.0 {
      return false
    }
  }
  multistate_consistency_report(data).length() == 0
}

///|
/// Return a transition count matrix through a horizon.
pub fn multistate_transition_matrix(
  data : MultiStateData,
  horizon : Double,
) -> Array[Array[Int]] {
  let matrix : Array[Array[Int]] = []
  for _ in 0..<data.state_count {
    matrix.push(Array::make(data.state_count, 0))
  }
  for transition in data.transitions {
    if transition.time <= horizon {
      matrix[transition.from_state][transition.to_state] += 1
    }
  }
  matrix
}
