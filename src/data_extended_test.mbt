///|
test "data validation returns structured errors" {
  match SurvivalData::try_new([1.0, 2.0], [1]) {
    Err(LengthMismatch) => ()
    _ => fail("length mismatch was not reported")
  }
  match SurvivalData::try_new([-1.0], [1]) {
    Err(NegativeTime) => ()
    _ => fail("negative time was not reported")
  }
  match SurvivalData::try_new([1.0], [2]) {
    Err(InvalidEvent(2)) => ()
    _ => fail("invalid event was not reported")
  }
}

///|
test "data summary and event table handle ties" {
  let data = SurvivalData::new([3.0, 1.0, 1.0, 4.0], [0, 1, 0, 1])
  let summary = data.summary()
  inspect(summary.n, content="4")
  inspect(summary.events, content="2")
  inspect(summary.censored, content="2")
  inspect(summary.min_time, content="1")
  inspect(summary.max_time, content="4")
  let table = event_table(data)
  inspect(table.length(), content="3")
  inspect(table[0].n_at_risk, content="4")
  inspect(table[0].n_events, content="1")
  inspect(table[0].n_censored, content="1")
  inspect(table[1].n_at_risk, content="2")
}

///|
test "data slicing and event/censor accessors" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  inspect(data.event_times().length(), content="2")
  inspect(data.event_times()[0], content="1")
  inspect(data.event_times()[1], content="3")
  inspect(data.censoring_times().length(), content="2")
  inspect(data.censoring_times()[0], content="2")
  inspect(data.censoring_times()[1], content="4")
  inspect(data.event_count(), content="2")
  inspect(data.censored_count(), content="2")
  inspect(data.between(2.0, 3.0).length(), content="2")
  inspect(data.to_array().length(), content="4")
}

///|
test "weighted observations preserve sorting and weights" {
  let weighted = WeightedSurvivalData::new([3.0, 1.0, 2.0], [1, 1, 0], [
    2.0, 1.0, 3.0,
  ])
  inspect(weighted.length(), content="3")
  inspect(weighted.get(0).time, content="1")
  inspect(weighted.get(0).weight, content="1")
  match WeightedSurvivalData::try_new([1.0], [1], [0.0]) {
    Err(InvalidTime) => ()
    _ => fail("non-positive weight was accepted")
  }
}

///|
test "panic constructors reject malformed input" {
  ignore(SurvivalData::new([1.0], [3]))
  ignore(WeightedSurvivalData::new([1.0], [1], [-1.0]))
}
