///|
test "survival accumulator supports checkpoints and merge" {
  let first = SurvivalAccumulator::new()
  first.add(3.0, 1)
  first.add_batch([1.0, 2.0], [1, 0])
  inspect(first.length(), content="3")
  inspect(first.event_count(), content="2")
  inspect(first.snapshot().get(0).time, content="1")
  let checkpoint = first.checkpoint("day-3", 3.0)
  inspect(checkpoint.label, content="day-3")
  inspect(checkpoint.events, content="2")
  let second = SurvivalAccumulator::new()
  second.add(4.0, 0)
  first.merge(second)
  inspect(first.length(), content="4")
  inspect(first.to_csv().contains("time,event"), content="true")
}

///|
test "hazard and survival transformations are inverse" {
  let hazards = [0.0, 0.2, 0.7, 1.5]
  let survival = hazard_to_survival(hazards)
  let restored = survival_to_hazard(survival)
  inspect(survival[0], content="1")
  inspect(survival[3] < survival[2], content="true")
  inspect(double_is_close(restored[1], hazards[1]), content="true")
  inspect(double_is_close(restored[3], hazards[3]), content="true")
}

///|
test "curve sampling interpolation and integration" {
  let fit = kaplan_meier(SurvivalData::new([1.0, 3.0], [1, 1]))
  let points = sample_km_curve(fit, 0.0, 4.0, 5)
  inspect(points.length(), content="5")
  inspect(points[0].value, content="1")
  inspect(interpolate_curve(points[0], points[1], 0.5) <= 1.0, content="true")
  inspect(
    integrate_curve([{ time: 0.0, value: 1.0 }, { time: 1.0, value: 0.0 }]),
    content="0.5",
  )
  let hazards = interval_hazards([
    { time: 0.0, value: 0.0 },
    { time: 1.0, value: 0.2 },
    { time: 2.0, value: 0.5 },
  ])
  inspect(hazards.length(), content="2")
  inspect(hazards[1].value, content="0.3")
  match
    inverse_curve([{ time: 0.0, value: 0.0 }, { time: 1.0, value: 0.5 }], 0.4) {
    Some(value) => inspect(value, content="1")
    None => fail("inverse curve did not find a crossing")
  }
}

///|
test "panic streaming invalid input" {
  let accumulator = SurvivalAccumulator::new()
  ignore(accumulator.add(-1.0, 1))
}
