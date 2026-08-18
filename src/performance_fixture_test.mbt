///|
test "medium fixture completes all core estimators" {
  let fixture = benchmark_fixture(250)
  let km = kaplan_meier(fixture.data)
  let na = nelson_aalen(fixture.data)
  let reliability = reliability_report(fixture.data, 20.0)
  inspect(km.steps.length(), content="250")
  inspect(na.steps.length(), content="250")
  inspect(reliability.sample_size, content="250")
  inspect(reliability.events > 0, content="true")
  let bootstrap = bootstrap_kaplan_meier(
    fixture.data,
    [5.0, 10.0, 20.0],
    10,
    77U,
  )
  inspect(bootstrap.replicates, content="10")
  inspect(bootstrap.mean.length(), content="3")
}

///|
test "cox medium fixture fits with stable coefficient dimensions" {
  let fixture = simulate_cox(80, [0.3, -0.25, 0.1], 0.05, 0.015, 17U)
  let fit = cox_fit(fixture.data, max_iterations=20)
  inspect(fit.coefficients.length(), content="3")
  inspect(fit.standard_errors.length(), content="3")
  inspect(fit.coefficients_table.length(), content="3")
  inspect(fit.baseline.length() > 0, content="true")
  let diagnostics = cox_ph_diagnostics(fixture.data, fit)
  inspect(diagnostics.length(), content="3")
}

///|
test "weighted and unweighted estimates expose comparable curves" {
  let unweighted = SurvivalData::new([1.0, 2.0, 3.0], [1, 0, 1])
  let weighted = WeightedSurvivalData::new([1.0, 2.0, 3.0], [1, 0, 1], [
    1.0, 1.0, 1.0,
  ])
  let km = kaplan_meier(unweighted)
  let weighted_km = weighted_kaplan_meier(weighted)
  inspect(
    double_is_close(km.predict(2.5), weighted_km.predict(2.5)),
    content="true",
  )
  let na = nelson_aalen(unweighted)
  let weighted_na = weighted_nelson_aalen(weighted)
  inspect(
    double_is_close(na.predict(2.5), weighted_na.predict(2.5)),
    content="true",
  )
}
