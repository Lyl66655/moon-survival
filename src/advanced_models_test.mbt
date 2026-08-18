///|
test "piecewise exponential fit accounts for exposure" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 5.0], [1, 0, 1, 1])
  let fit = piecewise_exponential_fit(data, [2.0, 4.0, 6.0])
  inspect(fit.segments.length(), content="3")
  inspect(fit.total_events, content="3")
  inspect(fit.total_exposure > 0.0, content="true")
  inspect(fit.hazard_at(1.0) >= 0.0, content="true")
  inspect(fit.cumulative_hazard_at(0.0), content="0")
  inspect(fit.survival_at(5.0) <= 1.0, content="true")
  inspect(fit.quantile(0.5) >= 0.0, content="true")
  inspect(fit.to_csv().contains("start,end"), content="true")
}

///|
test "stratified survival and actuarial life table" {
  let first = SurvivalData::new([1.0, 2.0, 4.0], [1, 0, 1])
  let second = SurvivalData::new([2.0, 3.0, 5.0], [1, 1, 0])
  let stratified = stratified_survival(
    ["control", "treated"],
    [first, second],
    4.0,
  )
  inspect(stratified.length(), content="2")
  inspect(stratified[0].stratum, content="control")
  inspect(stratified[1].events, content="2")
  let life_table = actuarial_life_table(first, [2.0, 4.0, 6.0])
  inspect(life_table.length(), content="3")
  inspect(life_table[0].n_entering, content="3")
  inspect(life_table[0].events, content="1")
  inspect(life_table[1].n_entering, content="1")
  inspect(life_table[0].cumulative_survival >= 0.0, content="true")
}

///|
test "piecewise model boundaries remain finite" {
  let data = SurvivalData::new([1.0, 1.0, 1.0], [0, 0, 0])
  let fit = piecewise_exponential_fit(data, [1.0, 2.0])
  inspect(fit.total_events, content="0")
  inspect(fit.survival_at(2.0), content="1")
  inspect(fit.hazard_at(10.0), content="0")
}

///|
test "panic piecewise model rejects unsorted cutpoints" {
  ignore(piecewise_exponential_fit(SurvivalData::new([1.0], [1]), [2.0, 1.0]))
}
