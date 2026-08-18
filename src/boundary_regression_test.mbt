///|
test "all events tied at one time remain mathematically bounded" {
  let data = SurvivalData::new([2.0, 2.0, 2.0, 2.0], [1, 1, 1, 1])
  let km = kaplan_meier(data)
  let na = nelson_aalen(data)
  inspect(km.steps.length(), content="1")
  inspect(km.predict(1.0), content="1")
  inspect(km.predict(2.0), content="0")
  inspect(na.predict(2.0), content="1")
  inspect(na.survival_at(2.0) > 0.0, content="true")
  inspect(km.restricted_mean(2.0), content="2")
}

///|
test "large finite values do not overflow estimators" {
  let data = SurvivalData::new([1000000.0, 2000000.0, 3000000.0], [1, 0, 1])
  let km = kaplan_meier(data)
  let na = nelson_aalen(data)
  inspect(km.predict(2500000.0) >= 0.0, content="true")
  inspect(na.predict(2500000.0) >= 0.0, content="true")
  let fit = weibull_fit(data)
  inspect(fit.scale > 0.0, content="true")
}

///|
test "empty and zero width analysis helpers" {
  let empty = SurvivalData::new([], [])
  inspect(assess_data_quality(empty).valid, content="false")
  inspect(analysis_time_grid(empty, 10).length(), content="0")
  inspect(
    curve_integrated_difference(
      kaplan_meier(empty),
      kaplan_meier(empty),
      10.0,
      10,
    ),
    content="0",
  )
  inspect(
    sample_km_curve(kaplan_meier(empty), 0.0, 1.0, 0).length(),
    content="0",
  )
  inspect(integrate_curve([]), content="0")
  match inverse_curve([], 0.5) {
    None => ()
    Some(_) => fail("empty inverse curve returned a value")
  }
}

///|
test "quality report identifies ties and censoring" {
  let data = SurvivalData::new([0.0, 1.0, 1.0, 2.0], [1, 0, 1, 0])
  let report = assess_data_quality(data)
  inspect(report.duplicate_time_count, content="1")
  inspect(report.zero_time_count, content="1")
  inspect(report.valid, content="true")
  inspect(report.warnings.length() >= 2, content="true")
  inspect(validate_predictions([0.0, 0.5, 1.0]), content="true")
  inspect(validate_predictions([-0.1]), content="false")
}

///|
test "curve comparison detects expected separation" {
  let left = kaplan_meier(SurvivalData::new([1.0, 2.0], [1, 1]))
  let right = kaplan_meier(SurvivalData::new([2.0, 3.0], [1, 1]))
  let grid = [0.0, 1.0, 2.0, 3.0]
  inspect(curve_max_difference(left, right, grid) > 0.0, content="true")
  inspect(
    curve_integrated_difference(left, right, 3.0, 20) > 0.0,
    content="true",
  )
  inspect(curve_is_monotone(left, grid), content="true")
  inspect(quality_flags_to_string(["a", "b"]), content="a,b")
}
