///|
test "kaplan_meier normal" {
  let times = [1.0, 2.0, 3.0, 4.0, 5.0]
  let events = [1, 1, 0, 1, 1]
  let data = SurvivalData::new(times, events)
  let km = kaplan_meier(data)
  
  inspect(km.steps.length(), content="5")
  inspect(km.predict(1.5), content="0.8")
  inspect(km.predict(2.5), content="0.6")
  // At t=3, censored, so survival probability does not drop, but at risk pool drops
  inspect(km.predict(3.5), content="0.6")
  // At t=4, 1 event out of 2 at risk. 0.6 * (1/2) = 0.3
  inspect(km.predict(4.5), content="0.3")
}
