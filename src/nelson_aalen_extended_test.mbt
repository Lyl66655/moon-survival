///|
test "nelson aalen variance and survival transform" {
  let data = SurvivalData::new([1.0, 2.0, 3.0], [1, 1, 0])
  let na = nelson_aalen_with_confidence(data)
  inspect(na.steps.length(), content="3")
  inspect(na.steps[0].cum_hazard, content="0.3333333333333333")
  inspect(na.steps[1].variance > 0.0, content="true")
  inspect(na.survival_at(1.5) < 1.0, content="true")
  inspect(na.predict_many([0.0, 1.5, 3.5]).length(), content="3")
  let ci = na.confidence_interval(2.5)
  inspect(ci.0 >= 0.0, content="true")
  inspect(ci.1 >= ci.0, content="true")
}

///|
test "nelson aalen censor-only steps do not add hazard" {
  let na = nelson_aalen(SurvivalData::new([1.0, 2.0, 3.0], [1, 0, 0]))
  inspect(na.predict(1.5), content="0.3333333333333333")
  inspect(na.predict(3.0), content="0.3333333333333333")
  inspect(na.increment_at(2.0), content="0")
  inspect(na.step_at(2.0).unwrap().n_censored, content="1")
}
