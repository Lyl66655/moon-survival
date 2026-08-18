///|
test "csv and markdown exports remain nonempty" {
  let data = SurvivalData::new([1.0, 2.0, 3.0], [1, 0, 1])
  inspect(kaplan_meier(data).to_csv().contains("lower_ci"), content="true")
  inspect(nelson_aalen(data).to_csv().contains("cum_hazard"), content="true")
  inspect(survival_data_to_csv(data).contains("time,event"), content="true")
  inspect(
    event_table_to_markdown(event_table(data)).contains("Censored"),
    content="true",
  )
  let fit = cox_fit(
    CoxData::new([1.0, 2.0, 3.0], [1, 1, 0], [[0.0], [1.0], [2.0]]),
  )
  inspect(fit.baseline_survival(2.0) >= 0.0, content="true")
}

///|
test "probability and quantile boundaries are clipped" {
  inspect(logit(0.0) < 0.0, content="true")
  inspect(logit(1.0) > 0.0, content="true")
  inspect(normal_quantile(0.0) < 0.0, content="true")
  inspect(normal_quantile(1.0) > 0.0, content="true")
  inspect(logrank_power(0, 2.0, 0.05), content="0")
  inspect(logrank_estimated_sample_size(10, 0.0) > 0, content="true")
}

///|
test "simulation and reliability agree on dimensions" {
  let fixture = benchmark_fixture(25)
  let report = reliability_report(fixture.data, 20.0)
  inspect(report.sample_size, content="25")
  inspect(fixture.latent_failure_times.length(), content="25")
  inspect(fixture.censoring_times.length(), content="25")
  let plan = make_power_plan(0.8, 0.05, 0.8, 0.5)
  inspect(plan.estimated_sample_size > 0, content="true")
}
