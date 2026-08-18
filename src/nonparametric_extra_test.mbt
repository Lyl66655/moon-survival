///|
test "nonparametric estimators agree on simple complete data" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 1, 1, 1])
  let km = kaplan_meier(data)
  let na = nelson_aalen(data)
  for time in [0.0, 1.0, 2.0, 3.0, 4.0] {
    let survival = km.predict(time)
    let transformed = na.survival_at(time)
    if survival < 0.0 ||
      survival > 1.0 ||
      transformed < 0.0 ||
      transformed > 1.0 {
      fail("nonparametric estimate left bounds")
    }
  }
  inspect(km.median().unwrap(), content="2")
  inspect(na.steps[3].n_events, content="1")
}

///|
test "ties and censoring are counted once per time step" {
  let data = SurvivalData::new([1.0, 1.0, 1.0, 2.0, 2.0, 3.0], [
    1, 0, 0, 1, 0, 0,
  ])
  let rows = event_table(data)
  inspect(rows.length(), content="3")
  inspect(rows[0].n_events, content="1")
  inspect(rows[0].n_censored, content="2")
  inspect(rows[1].n_at_risk, content="3")
  inspect(rows[1].n_events, content="1")
  inspect(rows[2].n_at_risk, content="1")
  inspect(kaplan_meier(data).steps.length(), content="3")
  inspect(nelson_aalen(data).steps[0].n_at_risk, content="6")
}

///|
test "weighted estimators respect weight scaling" {
  let unit = WeightedSurvivalData::new([1.0, 2.0, 3.0], [1, 1, 0], [
    1.0, 1.0, 1.0,
  ])
  let scaled = WeightedSurvivalData::new([1.0, 2.0, 3.0], [1, 1, 0], [
    5.0, 5.0, 5.0,
  ])
  let unit_km = weighted_kaplan_meier(unit)
  let scaled_km = weighted_kaplan_meier(scaled)
  inspect(
    double_is_close(unit_km.predict(2.5), scaled_km.predict(2.5)),
    content="true",
  )
  let unit_na = weighted_nelson_aalen(unit)
  let scaled_na = weighted_nelson_aalen(scaled)
  inspect(
    double_is_close(unit_na.predict(2.5), scaled_na.predict(2.5)),
    content="true",
  )
}

///|
test "confidence levels are normalized at boundaries" {
  let data = SurvivalData::new([1.0, 2.0], [1, 1])
  let low = kaplan_meier_with_confidence(data, confidence=-1.0)
  let high = kaplan_meier_with_confidence(data, confidence=2.0)
  inspect(
    low.confidence_level > 0.0 && low.confidence_level < 1.0,
    content="true",
  )
  inspect(
    high.confidence_level > 0.0 && high.confidence_level < 1.0,
    content="true",
  )
  inspect(low.steps[0].lower_ci >= 0.0, content="true")
  inspect(high.steps[0].upper_ci <= 1.0, content="true")
}
