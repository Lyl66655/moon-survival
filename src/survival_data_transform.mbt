///|
/// A transform operation with a named audit label.
pub(all) struct SurvivalTransformStep {
  name : String
  input_rows : Int
  output_rows : Int
  removed_rows : Int
  transformed : Bool
} derive(Eq, Debug)

///|
/// A transformed dataset and its operation history.
pub(all) struct SurvivalTransformedData {
  data : SurvivalData
  steps : Array[SurvivalTransformStep]
  transform_version : String
} derive(Eq, Debug)

///|
/// Start a transformation chain.
pub fn survival_transform_start(data : SurvivalData) -> SurvivalTransformedData {
  { data, steps: [], transform_version: "1" }
}

///|
/// Add an operation step.
fn survival_transform_add_step(
  chain : SurvivalTransformedData,
  name : String,
  input_rows : Int,
  output_rows : Int,
) -> SurvivalTransformedData {
  let steps : Array[SurvivalTransformStep] = []
  for step in chain.steps {
    steps.push(step)
  }
  steps.push({
    name,
    input_rows,
    output_rows,
    removed_rows: input_rows - output_rows,
    transformed: true,
  })
  { data: chain.data, steps, transform_version: chain.transform_version }
}

///|
/// Remove rows beyond a follow-up horizon.
pub fn survival_transform_truncate(
  chain : SurvivalTransformedData,
  horizon : Double,
) -> SurvivalTransformedData {
  let source = chain.data
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in source.to_array() {
    if observation.time <= horizon {
      times.push(observation.time)
      events.push(observation.event)
    } else {
      times.push(horizon)
      events.push(0)
    }
  }
  let data = SurvivalData::new(times, events)
  let next = {
    data,
    steps: chain.steps,
    transform_version: chain.transform_version,
  }
  survival_transform_add_step(next, "truncate", source.length(), data.length())
}

///|
/// Recode an event to a binary indicator.
pub fn survival_transform_recode_event(
  chain : SurvivalTransformedData,
  event_code : Int,
) -> SurvivalTransformedData {
  let source = chain.data
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in source.to_array() {
    times.push(observation.time)
    events.push(if observation.event == event_code { 1 } else { 0 })
  }
  let data = SurvivalData::new(times, events)
  survival_transform_add_step(
    { data, steps: chain.steps, transform_version: chain.transform_version },
    "recode_event",
    source.length(),
    data.length(),
  )
}

///|
/// Sort and deduplicate exact duplicate rows.
pub fn survival_transform_deduplicate(
  chain : SurvivalTransformedData,
) -> SurvivalTransformedData {
  let source = chain.data
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in source.to_array() {
    if times.length() == 0 ||
      times[times.length() - 1] != observation.time ||
      events[events.length() - 1] != observation.event {
      times.push(observation.time)
      events.push(observation.event)
    }
  }
  let data = SurvivalData::new(times, events)
  survival_transform_add_step(
    { data, steps: chain.steps, transform_version: chain.transform_version },
    "deduplicate",
    source.length(),
    data.length(),
  )
}

///|
/// Apply an administrative entry delay.
pub fn survival_transform_delay_entry(
  chain : SurvivalTransformedData,
  delay : Double,
) -> SurvivalTransformedData {
  let source = chain.data
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in source.to_array() {
    let time = observation.time - delay
    if time > 0.0 {
      times.push(time)
      events.push(observation.event)
    }
  }
  let data = SurvivalData::new(times, events)
  survival_transform_add_step(
    { data, steps: chain.steps, transform_version: chain.transform_version },
    "delay_entry",
    source.length(),
    data.length(),
  )
}

///|
/// Convert durations to another unit.
pub fn survival_transform_scale_time(
  chain : SurvivalTransformedData,
  multiplier : Double,
) -> SurvivalTransformedData {
  let source = chain.data
  let factor = if multiplier <= 0.0 { 1.0 } else { multiplier }
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in source.to_array() {
    times.push(observation.time * factor)
    events.push(observation.event)
  }
  let data = SurvivalData::new(times, events)
  survival_transform_add_step(
    { data, steps: chain.steps, transform_version: chain.transform_version },
    "scale_time",
    source.length(),
    data.length(),
  )
}

///|
/// Stratify an analysis by dropping a specified group.
pub fn survival_transform_drop_group(
  chain : SurvivalTransformedData,
  groups : Array[Int],
  drop_group : Int,
) -> SurvivalTransformedData {
  let source = chain.data
  let times : Array[Double] = []
  let events : Array[Int] = []
  for i in 0..<source.length() {
    if i >= groups.length() || groups[i] != drop_group {
      times.push(source.get(i).time)
      events.push(source.get(i).event)
    }
  }
  let data = SurvivalData::new(times, events)
  survival_transform_add_step(
    { data, steps: chain.steps, transform_version: chain.transform_version },
    "drop_group",
    source.length(),
    data.length(),
  )
}

///|
/// Keep only rows in a time window.
pub fn survival_transform_window(
  chain : SurvivalTransformedData,
  start : Double,
  end : Double,
) -> SurvivalTransformedData {
  let source = chain.data
  let times : Array[Double] = []
  let events : Array[Int] = []
  for observation in source.to_array() {
    if observation.time >= start && observation.time <= end {
      times.push(observation.time)
      events.push(observation.event)
    }
  }
  let data = SurvivalData::new(times, events)
  survival_transform_add_step(
    { data, steps: chain.steps, transform_version: chain.transform_version },
    "window",
    source.length(),
    data.length(),
  )
}

///|
/// Return the number of transformation steps.
pub fn SurvivalTransformedData::step_count(
  self : SurvivalTransformedData,
) -> Int {
  self.steps.length()
}

///|
/// Return the last operation name.
pub fn SurvivalTransformedData::last_step(
  self : SurvivalTransformedData,
) -> String {
  if self.steps.length() == 0 {
    "none"
  } else {
    self.steps[self.steps.length() - 1].name
  }
}

///|
/// Compute a transformed-data checksum.
pub fn SurvivalTransformedData::checksum(
  self : SurvivalTransformedData,
) -> String {
  survival_export_checksum(survival_data_to_csv(self.data))
}

///|
/// Export transformation history.
pub fn survival_transform_steps_to_csv(
  steps : Array[SurvivalTransformStep],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("name,input_rows,output_rows,removed_rows\n")
  for step in steps {
    builder.write_string(step.name)
    builder.write_char(',')
    builder.write_string(step.input_rows.to_string())
    builder.write_char(',')
    builder.write_string(step.output_rows.to_string())
    builder.write_char(',')
    builder.write_string(step.removed_rows.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate a transform chain.
pub fn survival_transform_is_valid(chain : SurvivalTransformedData) -> Bool {
  if chain.transform_version == "" {
    return false
  }
  for observation in chain.data.to_array() {
    if observation.time < 0.0 ||
      observation.time.is_nan() ||
      observation.time.is_inf() ||
      (observation.event != 0 && observation.event != 1) {
      return false
    }
  }
  for step in chain.steps {
    if step.input_rows < step.output_rows || step.output_rows < 0 {
      return false
    }
  }
  true
}
