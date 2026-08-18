///|
test "kaplan meier confidence and quantiles" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 1, 0, 1])
  let km = kaplan_meier_with_confidence(data, confidence=0.90)
  inspect(km.sample_size, content="4")
  inspect(km.event_count, content="3")
  inspect(km.observed_event_steps(), content="3")
  match km.median() {
    Some(value) => inspect(value, content="2")
    None => fail("median should be estimable")
  }
  match km.quantile(0.25) {
    Some(value) => inspect(value, content="1")
    None => fail("lower quantile should be estimable")
  }
  inspect(km.confidence_interval(2.0).0 >= 0.0, content="true")
  inspect(km.confidence_interval(2.0).1 <= 1.0, content="true")
}

///|
test "kaplan meier restricted mean and predictions" {
  let km = kaplan_meier(SurvivalData::new([2.0, 4.0], [1, 1]))
  inspect(km.predict(0.0), content="1")
  inspect(km.predict(-1.0), content="1")
  inspect(km.predict(2.0), content="0.5")
  inspect(km.restricted_mean(1.0), content="1")
  inspect(km.restricted_mean(3.0), content="2.5")
  inspect(km.restricted_mean_after(2.0, 4.0), content="1")
  inspect(km.hazard_at(2.0), content="0.5")
  let predictions = km.predict_many([0.0, 2.0, 5.0])
  inspect(predictions.length(), content="3")
  inspect(predictions[0], content="1")
  inspect(predictions[1], content="0.5")
  inspect(predictions[2], content="0")
}

///|
test "kaplan meier all censored and empty data" {
  let empty = kaplan_meier(SurvivalData::new([], []))
  inspect(empty.steps.length(), content="0")
  inspect(empty.predict(10.0), content="1")
  let censored = kaplan_meier(SurvivalData::new([1.0, 2.0], [0, 0]))
  inspect(censored.steps.length(), content="2")
  inspect(censored.predict(3.0), content="1")
  match censored.median() {
    None => ()
    Some(_) => fail("all-censored curve has a median")
  }
}

///|
test "weighted nonparametric estimators" {
  let data = WeightedSurvivalData::new([1.0, 2.0], [1, 1], [1.0, 3.0])
  let km = weighted_kaplan_meier(data)
  let na = weighted_nelson_aalen(data)
  inspect(km.predict(0.5), content="1")
  inspect(km.predict(1.5), content="0.75")
  inspect(na.predict(1.5), content="0.25")
}
