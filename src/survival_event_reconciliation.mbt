///|
/// A discrepancy between two event sources.
pub(all) struct SurvivalEventDiscrepancy {
  subject : Int
  first_time : Double
  second_time : Double
  first_event : Int
  second_event : Int
  reason : String
  discrepancy_row : Bool
} derive(Eq, Debug)

///|
/// Reconcile two event-coded sources by row identifier.
pub fn survival_reconcile_events(
  first_times : Array[Double],
  first_events : Array[Int],
  second_times : Array[Double],
  second_events : Array[Int],
  tolerance : Double,
) -> Array[SurvivalEventDiscrepancy] {
  let result : Array[SurvivalEventDiscrepancy] = []
  let count = if first_times.length() < second_times.length() {
    first_times.length()
  } else {
    second_times.length()
  }
  for i in 0..<count {
    let first_event = if i < first_events.length() {
      first_events[i]
    } else {
      -1
    }
    let second_event = if i < second_events.length() {
      second_events[i]
    } else {
      -1
    }
    let time_difference = (first_times[i] - second_times[i]).abs()
    if time_difference > tolerance || first_event != second_event {
      result.push({
        subject: i,
        first_time: first_times[i],
        second_time: second_times[i],
        first_event,
        second_event,
        reason: if first_event != second_event {
          "event_code"
        } else {
          "time_difference"
        },
        discrepancy_row: true,
      })
    }
  }
  result
}

///|
/// Count disagreements in event codes.
pub fn survival_event_disagreement_count(
  discrepancies : Array[SurvivalEventDiscrepancy],
) -> Int {
  let mut count = 0
  for discrepancy in discrepancies {
    if discrepancy.reason == "event_code" {
      count += 1
    }
  }
  count
}

///|
/// Count time disagreements.
pub fn survival_time_disagreement_count(
  discrepancies : Array[SurvivalEventDiscrepancy],
) -> Int {
  let mut count = 0
  for discrepancy in discrepancies {
    if discrepancy.reason == "time_difference" {
      count += 1
    }
  }
  count
}

///|
/// Compute the maximum time discrepancy.
pub fn survival_max_event_time_difference(
  discrepancies : Array[SurvivalEventDiscrepancy],
) -> Double {
  let mut maximum = 0.0
  for discrepancy in discrepancies {
    let difference = (discrepancy.first_time - discrepancy.second_time).abs()
    if difference > maximum {
      maximum = difference
    }
  }
  maximum
}

///|
/// Create a reconciled dataset preferring the first source.
pub fn survival_reconciled_data(
  times : Array[Double],
  events : Array[Int],
) -> SurvivalData {
  SurvivalData::new(times, events)
}

///|
/// Export discrepancies.
pub fn survival_event_discrepancies_to_csv(
  discrepancies : Array[SurvivalEventDiscrepancy],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "subject,first_time,second_time,first_event,second_event,reason\n",
  )
  for discrepancy in discrepancies {
    builder.write_string(discrepancy.subject.to_string())
    builder.write_char(',')
    builder.write_string(discrepancy.first_time.to_string())
    builder.write_char(',')
    builder.write_string(discrepancy.second_time.to_string())
    builder.write_char(',')
    builder.write_string(discrepancy.first_event.to_string())
    builder.write_char(',')
    builder.write_string(discrepancy.second_event.to_string())
    builder.write_char(',')
    builder.write_string(discrepancy.reason)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate discrepancy fields.
pub fn survival_event_discrepancies_are_valid(
  discrepancies : Array[SurvivalEventDiscrepancy],
) -> Bool {
  for discrepancy in discrepancies {
    if discrepancy.subject < 0 ||
      discrepancy.first_time < 0.0 ||
      discrepancy.second_time < 0.0 {
      return false
    }
  }
  true
}
