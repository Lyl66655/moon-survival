///|
test "model selection scores and likelihood ratio" {
  let first = model_score("exponential", 1, 100, -50.0)
  let second = model_score("weibull", 2, 100, -45.0)
  inspect(first.aic, content="102")
  inspect(second.aic < first.aic, content="true")
  match select_by_aic([first, second]) {
    Some(value) => inspect(value.name, content="weibull")
    None => fail("aic selection returned no model")
  }
  let ratio = likelihood_ratio_test(-50.0, -45.0, 1, 2)
  inspect(ratio.degrees_of_freedom, content="1")
  inspect(ratio.statistic, content="10")
  inspect(ratio.larger_model_preferred, content="true")
  inspect(objective_trace_is_monotone([0.0, 1.0, 1.0]), content="true")
  inspect(objective_trace_is_monotone([1.0, 0.0]), content="false")
}

///|
test "reports and rendering expose production metrics" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 1])
  let report = build_survival_report("demo", data, 4.0)
  inspect(report.summary.n, content="4")
  inspect(report.is_analysis_ready(), content="true")
  inspect(report.to_markdown().contains("Sample size"), content="true")
  inspect(report.to_csv().contains("restricted_mean"), content="true")
  let comparison = report.with_logrank(logrank_test(data, data))
  match comparison.logrank {
    Some(value) => inspect(value.events, content="6")
    None => fail("logrank was not attached")
  }
  let empty = build_survival_report("empty", SurvivalData::new([], []), 1.0)
  inspect(empty.is_analysis_ready(), content="false")
}

///|
test "power planning is monotone in target power" {
  let plan = make_power_plan(0.65, 0.05, 0.8, 0.6, loss_fraction=0.1)
  inspect(plan.required_events > 0, content="true")
  inspect(plan.estimated_sample_size >= plan.required_events, content="true")
  inspect(
    logrank_required_events(0.65, 0.05, 0.9) >= plan.required_events,
    content="true",
  )
  inspect(
    logrank_power(plan.required_events, 0.65, 0.05) >= 0.0,
    content="true",
  )
  inspect(detectable_hazard_ratio(100, 0.05, 0.8) > 1.0, content="true")
  inspect(plan.to_markdown().contains("Required events"), content="true")
}

///|
test "visualization helpers produce stable text" {
  let km = kaplan_meier(SurvivalData::new([1.0, 2.0, 3.0], [1, 1, 0]))
  let line = km.to_ascii(width=10)
  inspect(line.contains("survival"), content="true")
  inspect(sparkline([0.0, 0.5, 1.0], width=5).length(), content="3")
  inspect(doubles_to_csv([1.0, 2.0]), content="1,2")
  inspect(
    event_table_to_markdown(event_table(SurvivalData::new([1.0], [1]))).contains(
      "At risk",
    ),
    content="true",
  )
  let competing = cumulative_incidence(
    CompetingRiskData::new([1.0, 2.0], [1, 2]),
  )
  inspect(competing.to_markdown().contains("Incidence"), content="true")
}
