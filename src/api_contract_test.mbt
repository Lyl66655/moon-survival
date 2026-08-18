///|
test "public constructors retain original API behavior" {
  let data = SurvivalData::new([5.0, 1.0, 3.0], [0, 1, 1])
  inspect(data.length(), content="3")
  inspect(data.get(0).time, content="1")
  inspect(data.get(2).event, content="0")
  let km = kaplan_meier(data)
  let na = nelson_aalen(data)
  inspect(km.predict(0.0), content="1")
  inspect(na.predict(0.0), content="0")
  let cox = CoxData::new([1.0, 2.0], [1, 0], [[0.1], [0.2]])
  inspect(cox.n_covariates, content="1")
  inspect(cox_partial_likelihood(cox, [0.0]) > 0.0, content="true")
}

///|
test "all major estimators accept singleton data" {
  let data = SurvivalData::new([1.0], [1])
  inspect(kaplan_meier(data).predict(1.0), content="0")
  inspect(nelson_aalen(data).predict(1.0), content="1")
  inspect(exponential_fit(data).events, content="1")
  inspect(weibull_fit(data).events, content="1")
  inspect(piecewise_exponential_fit(data, [2.0]).total_events, content="1")
  inspect(build_survival_report("one", data, 1.0).summary.n, content="1")
}

///|
test "prediction arrays preserve requested order" {
  let fit = kaplan_meier(SurvivalData::new([1.0, 2.0, 3.0], [1, 1, 1]))
  let values = fit.predict_many([3.0, 0.0, 2.0, 1.0])
  inspect(values.length(), content="4")
  inspect(values[0], content="0")
  inspect(values[1], content="1")
  inspect(values[2], content="0.33333333333333337")
  inspect(values[3], content="0.6666666666666667")
}

///|
test "named CSV columns report missing data cleanly" {
  match survival_data_from_csv("time,event\n1,1\n", "missing", "event") {
    Err(MissingColumn("missing")) => ()
    _ => fail("missing time column was not reported")
  }
  match
    cox_data_from_csv("time,event,x\n1,1,2\n", "time", "event", ["missing"]) {
    Err(MissingColumn("missing")) => ()
    _ => fail("missing covariate was not reported")
  }
}

///|
test "result formatting is deterministic" {
  let data = SurvivalData::new([1.0, 2.0], [1, 0])
  let report = build_survival_report("stable", data, 2.0)
  inspect(report.to_markdown().has_prefix("# stable"), content="true")
  inspect(report.to_csv().has_prefix("metric,value"), content="true")
  inspect(
    kaplan_meier(data).to_ascii(width=8).has_prefix("survival"),
    content="true",
  )
}
