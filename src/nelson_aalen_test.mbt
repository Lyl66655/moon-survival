///|
test "nelson_aalen normal" {
  let times = [1.0, 2.0, 3.0, 4.0, 5.0]
  let events = [1, 1, 0, 1, 1]
  let data = SurvivalData::new(times, events)
  let na = nelson_aalen(data)
  
  inspect(na.steps.length(), content="5")
  inspect(na.predict(1.5), content="0.2")
  inspect(na.predict(2.5), content="0.45")
  // At t=3, censored, so hazard does not increase
  inspect(na.predict(3.5), content="0.45")
  // At t=4, 1 event out of 2 at risk. 0.45 + 0.5 = 0.95
  inspect(na.predict(4.5), content="0.95")
}
