///|
test "SurvivalData::new normal" {
  let times = [5.0, 1.0, 3.0, 3.0, 2.0]
  let events = [0, 1, 0, 1, 1]
  let data = SurvivalData::new(times, events)
  
  // Sorted order should be:
  // t=1.0, e=1
  // t=2.0, e=1
  // t=3.0, e=1
  // t=3.0, e=0
  // t=5.0, e=0
  inspect!(data.length(), content="5")
  inspect!(data.get(0).time, content="1.0")
  inspect!(data.get(0).event, content="1")
  inspect!(data.get(2).time, content="3.0")
  inspect!(data.get(2).event, content="1")
  inspect!(data.get(3).time, content="3.0")
  inspect!(data.get(3).event, content="0")
}

///|
test "panic SurvivalData::new length mismatch" {
  ignore(SurvivalData::new([1.0, 2.0], [1]))
}

///|
test "panic SurvivalData::new invalid event" {
  ignore(SurvivalData::new([1.0], [2]))
}
