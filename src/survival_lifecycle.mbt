///|
/// A lifecycle state for operational survival monitoring.
pub(all) enum SurvivalLifecycleState {
  Enrolled
  AtRisk
  EventObserved
  Censored
  Withdrawn
} derive(Eq, Debug)

///|
/// A lifecycle state transition.
pub(all) struct SurvivalLifecycleTransition {
  subject : Int
  time : Double
  from_state : SurvivalLifecycleState
  to_state : SurvivalLifecycleState
  reason : String
  lifecycle_row : Bool
} derive(Eq, Debug)

///|
/// A lifecycle trace.
pub(all) struct SurvivalLifecycleTrace {
  subject : Int
  transitions : Array[SurvivalLifecycleTransition]
  current_state : SurvivalLifecycleState
  trace_row : Bool
} derive(Eq, Debug)

///|
/// Convert core data to lifecycle traces.
pub fn survival_lifecycle_traces(
  data : SurvivalData,
) -> Array[SurvivalLifecycleTrace] {
  let result : Array[SurvivalLifecycleTrace] = []
  for i in 0..<data.length() {
    let observation = data.get(i)
    let target = if observation.event == 1 { EventObserved } else { Censored }
    let transitions : Array[SurvivalLifecycleTransition] = [
      {
        subject: i,
        time: 0.0,
        from_state: Enrolled,
        to_state: AtRisk,
        reason: "enrollment",
        lifecycle_row: true,
      },
      {
        subject: i,
        time: observation.time,
        from_state: AtRisk,
        to_state: target,
        reason: if observation.event == 1 {
          "event"
        } else {
          "censoring"
        },
        lifecycle_row: true,
      },
    ]
    result.push({
      subject: i,
      transitions,
      current_state: target,
      trace_row: true,
    })
  }
  result
}

///|
/// Return a state label.
pub fn survival_lifecycle_state_name(state : SurvivalLifecycleState) -> String {
  match state {
    Enrolled => "enrolled"
    AtRisk => "at_risk"
    EventObserved => "event_observed"
    Censored => "censored"
    Withdrawn => "withdrawn"
  }
}

///|
/// Count traces by current state.
pub fn survival_lifecycle_state_counts(
  traces : Array[SurvivalLifecycleTrace],
) -> Array[Int] {
  let result = Array::make(5, 0)
  for trace in traces {
    match trace.current_state {
      Enrolled => result[0] += 1
      AtRisk => result[1] += 1
      EventObserved => result[2] += 1
      Censored => result[3] += 1
      Withdrawn => result[4] += 1
    }
  }
  result
}

///|
/// Return traces that remain under follow-up.
pub fn survival_lifecycle_at_risk(
  traces : Array[SurvivalLifecycleTrace],
) -> Array[SurvivalLifecycleTrace] {
  let result : Array[SurvivalLifecycleTrace] = []
  for trace in traces {
    if trace.current_state == AtRisk {
      result.push(trace)
    }
  }
  result
}

///|
/// Compute median event time from lifecycle traces.
pub fn survival_lifecycle_median_time(
  traces : Array[SurvivalLifecycleTrace],
) -> Double {
  let times : Array[Double] = []
  for trace in traces {
    for transition in trace.transitions {
      if transition.to_state == EventObserved {
        times.push(transition.time)
      }
    }
  }
  survival_ext_quantile(times, 0.5)
}

///|
/// Export lifecycle transitions.
pub fn survival_lifecycle_to_csv(
  traces : Array[SurvivalLifecycleTrace],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("subject,time,from_state,to_state,reason\n")
  for trace in traces {
    for transition in trace.transitions {
      builder.write_string(transition.subject.to_string())
      builder.write_char(',')
      builder.write_string(transition.time.to_string())
      builder.write_char(',')
      builder.write_string(survival_lifecycle_state_name(transition.from_state))
      builder.write_char(',')
      builder.write_string(survival_lifecycle_state_name(transition.to_state))
      builder.write_char(',')
      builder.write_string(transition.reason)
      builder.write_char('\n')
    }
  }
  builder.to_string()
}

///|
/// Validate lifecycle traces.
pub fn survival_lifecycle_is_valid(
  traces : Array[SurvivalLifecycleTrace],
) -> Bool {
  for trace in traces {
    if trace.subject < 0 || trace.transitions.length() == 0 {
      return false
    }
    for transition in trace.transitions {
      if transition.time < 0.0 {
        return false
      }
    }
  }
  true
}
