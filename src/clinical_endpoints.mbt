///|
/// A prespecified clinical endpoint definition.
pub(all) struct ClinicalEndpointDefinition {
  name : String
  horizon : Double
  event_codes : Array[Int]
  competing_codes : Array[Int]
  endpoint_row : Bool
} derive(Eq, Debug)

///|
/// One adjudicated endpoint record.
pub(all) struct ClinicalEndpointRecord {
  subject : Int
  time : Double
  event : Int
  endpoint : Int
  reason : String
  endpoint_record : Bool
} derive(Eq, Debug)

///|
/// Clinical endpoint summary.
pub(all) struct ClinicalEndpointSummary {
  name : String
  sample_size : Int
  events : Int
  competing_events : Int
  censored : Int
  incidence : Double
  median_event_time : Double
} derive(Eq, Debug)

///|
/// Define a composite endpoint from allowed event codes.
pub fn clinical_endpoint_definition(
  name : String,
  horizon : Double,
  event_codes : Array[Int],
  competing_codes : Array[Int],
) -> ClinicalEndpointDefinition {
  {
    name,
    horizon: if horizon < 0.0 {
      0.0
    } else {
      horizon
    },
    event_codes,
    competing_codes,
    endpoint_row: true,
  }
}

///|
/// Test membership in a code set.
fn clinical_code_contains(codes : Array[Int], value : Int) -> Bool {
  for code in codes {
    if code == value {
      return true
    }
  }
  false
}

///|
/// Adjudicate a raw cause code against an endpoint definition.
pub fn clinical_adjudicate(
  definition : ClinicalEndpointDefinition,
  subject : Int,
  time : Double,
  cause : Int,
) -> ClinicalEndpointRecord {
  let endpoint = if clinical_code_contains(definition.event_codes, cause) {
    1
  } else if clinical_code_contains(definition.competing_codes, cause) {
    2
  } else {
    0
  }
  let reason = if endpoint == 1 {
    "endpoint"
  } else if endpoint == 2 {
    "competing"
  } else {
    "censored"
  }
  { subject, time, event: cause, endpoint, reason, endpoint_record: true }
}

///|
/// Adjudicate a parallel cause-coded dataset.
pub fn clinical_adjudicate_all(
  definition : ClinicalEndpointDefinition,
  subjects : Array[Int],
  times : Array[Double],
  causes : Array[Int],
) -> Array[ClinicalEndpointRecord] {
  let result : Array[ClinicalEndpointRecord] = []
  let count = if subjects.length() < times.length() {
    subjects.length()
  } else {
    times.length()
  }
  let bounded = if count < causes.length() { count } else { causes.length() }
  for i in 0..<bounded {
    result.push(
      clinical_adjudicate(definition, subjects[i], times[i], causes[i]),
    )
  }
  result
}

///|
/// Summarize adjudicated records.
pub fn clinical_endpoint_summary(
  definition : ClinicalEndpointDefinition,
  records : Array[ClinicalEndpointRecord],
) -> ClinicalEndpointSummary {
  let mut events = 0
  let mut competing = 0
  let mut censored = 0
  let event_times : Array[Double] = []
  for record in records {
    if record.endpoint == 1 {
      events += 1
      event_times.push(record.time)
    } else if record.endpoint == 2 {
      competing += 1
    } else {
      censored += 1
    }
  }
  event_times.sort()
  {
    name: definition.name,
    sample_size: records.length(),
    events,
    competing_events: competing,
    censored,
    incidence: if records.length() == 0 {
      0.0
    } else {
      events.to_double() / records.length().to_double()
    },
    median_event_time: if event_times.length() == 0 {
      0.0
    } else {
      event_times[event_times.length() / 2]
    },
  }
}

///|
/// Construct SurvivalData for a composite endpoint.
pub fn clinical_to_survival_data(
  definition : ClinicalEndpointDefinition,
  records : Array[ClinicalEndpointRecord],
) -> SurvivalData {
  let times : Array[Double] = []
  let events : Array[Int] = []
  for record in records {
    let observed_time = if record.time > definition.horizon &&
      definition.horizon > 0.0 {
      definition.horizon
    } else {
      record.time
    }
    times.push(observed_time)
    events.push(
      if record.endpoint == 1 && record.time <= definition.horizon {
        1
      } else {
        0
      },
    )
  }
  SurvivalData::new(times, events)
}

///|
/// Calculate a restricted-mean endpoint contrast.
pub fn clinical_rmst_contrast(
  first : SurvivalData,
  second : SurvivalData,
  horizon : Double,
) -> Double {
  kaplan_meier(first).restricted_mean(horizon) -
  kaplan_meier(second).restricted_mean(horizon)
}

///|
/// Calculate an endpoint incidence contrast.
pub fn clinical_incidence_contrast(
  first : SurvivalData,
  second : SurvivalData,
  horizon : Double,
) -> Double {
  let first_risk = 1.0 - kaplan_meier(first).survival_at(horizon)
  let second_risk = 1.0 - kaplan_meier(second).survival_at(horizon)
  first_risk - second_risk
}

///|
/// Identify responders under an event-free horizon.
pub fn clinical_responder_flags(
  data : SurvivalData,
  horizon : Double,
) -> Array[Int] {
  let result : Array[Int] = []
  for observation in data.to_array() {
    result.push(
      if observation.event == 1 && observation.time <= horizon {
        0
      } else {
        1
      },
    )
  }
  result
}

///|
/// Apply a landmark window to endpoint records.
pub fn clinical_landmark_records(
  records : Array[ClinicalEndpointRecord],
  landmark : Double,
  horizon : Double,
) -> Array[ClinicalEndpointRecord] {
  let result : Array[ClinicalEndpointRecord] = []
  for record in records {
    if record.time >= landmark && record.time <= horizon {
      result.push(record)
    }
  }
  result
}

///|
/// Count endpoint records by reason.
pub fn clinical_reason_counts(
  records : Array[ClinicalEndpointRecord],
) -> Array[Int] {
  let result = Array::make(3, 0)
  for record in records {
    if record.endpoint >= 0 && record.endpoint < 3 {
      result[record.endpoint] += 1
    }
  }
  result
}

///|
/// Export adjudicated endpoint rows.
pub fn clinical_endpoint_records_to_csv(
  records : Array[ClinicalEndpointRecord],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("subject,time,event,endpoint,reason\n")
  for record in records {
    builder.write_string(record.subject.to_string())
    builder.write_char(',')
    builder.write_string(record.time.to_string())
    builder.write_char(',')
    builder.write_string(record.event.to_string())
    builder.write_char(',')
    builder.write_string(record.endpoint.to_string())
    builder.write_char(',')
    builder.write_string(record.reason)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export endpoint summary.
pub fn clinical_endpoint_summary_to_csv(
  summary : ClinicalEndpointSummary,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("name,")
  builder.write_string(summary.name)
  builder.write_char('\n')
  builder.write_string("sample_size,")
  builder.write_string(summary.sample_size.to_string())
  builder.write_char('\n')
  builder.write_string("events,")
  builder.write_string(summary.events.to_string())
  builder.write_char('\n')
  builder.write_string("competing_events,")
  builder.write_string(summary.competing_events.to_string())
  builder.write_char('\n')
  builder.write_string("censored,")
  builder.write_string(summary.censored.to_string())
  builder.write_char('\n')
  builder.write_string("incidence,")
  builder.write_string(summary.incidence.to_string())
  builder.write_char('\n')
  builder.write_string("median_event_time,")
  builder.write_string(summary.median_event_time.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Check endpoint definition and records.
pub fn clinical_endpoint_is_valid(
  definition : ClinicalEndpointDefinition,
  records : Array[ClinicalEndpointRecord],
) -> Bool {
  if definition.name == "" || definition.horizon < 0.0 {
    return false
  }
  for record in records {
    if record.subject < 0 ||
      record.time < 0.0 ||
      record.endpoint < 0 ||
      record.endpoint > 2 {
      return false
    }
  }
  true
}
