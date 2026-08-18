///|
test "reliability fits provide monotone functions" {
  let data = SurvivalData::new([1.0, 2.0, 4.0, 5.0, 7.0], [1, 1, 0, 1, 0])
  let exponential = exponential_fit(data)
  let weibull = weibull_fit(data)
  inspect(exponential.rate > 0.0, content="true")
  inspect(exponential.survival_at(0.0), content="1")
  inspect(exponential.survival_at(5.0) < 1.0, content="true")
  inspect(weibull.scale > 0.0, content="true")
  inspect(weibull.survival_at(0.0), content="1")
  inspect(weibull.survival_at(5.0) < 1.0, content="true")
  inspect(weibull.hazard_at(5.0) > 0.0, content="true")
  inspect(weibull.quantile(0.5) > 0.0, content="true")
  inspect(reliability_report(data, 5.0).events, content="3")
}

///|
test "numeric summaries and transformations" {
  let values = [1.0, 2.0, 3.0, 4.0]
  let summary = numeric_summary(values)
  inspect(summary.count, content="4")
  inspect(summary.mean, content="2.5")
  inspect(summary.median, content="2.5")
  inspect(summary.minimum, content="1")
  inspect(summary.maximum, content="4")
  inspect(sample_variance(values), content="1.6666666666666667")
  inspect(correlation(values, values), content="1")
  inspect(covariance(values, values), content="1.6666666666666667")
  inspect(quantile_doubles(values, 0.25).unwrap(), content="1.75")
}

///|
test "standardization preserves shape" {
  let standardization = standardize_matrix([
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
  ])
  inspect(standardization.values.length(), content="3")
  inspect(standardization.values[0].length(), content="2")
  let transformed = standardization.transform([4.0, 40.0])
  inspect(transformed[0] > standardization.values[2][0], content="true")
  inspect(column_means([[1.0, 2.0], [3.0, 4.0]])[1], content="3")
}

///|
test "event rate bins and brier score are bounded" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  let bins = event_rate_bins(data, 0.0, 4.0, 4)
  inspect(bins.length(), content="4")
  let score = brier_score(data, [0.8, 0.7, 0.5, 0.2], 3.0)
  inspect(score >= 0.0 && score <= 1.0, content="true")
  let integrated = integrated_brier_score(
    data,
    [[0.8, 0.7, 0.5, 0.2], [0.7, 0.6, 0.4, 0.1]],
    [2.0, 3.0],
  )
  inspect(integrated >= 0.0 && integrated <= 1.0, content="true")
}
