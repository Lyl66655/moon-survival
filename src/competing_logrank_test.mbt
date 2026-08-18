///|
test "competing risks cumulative incidence sums with survival" {
  let data = CompetingRiskData::new([1.0, 1.0, 2.0, 3.0, 4.0], [1, 2, 0, 1, 2])
  let result = cumulative_incidence(data)
  inspect(result.n_causes, content="2")
  inspect(result.steps.length(), content="4")
  inspect(result.predict(1, 1.0) > 0.0, content="true")
  inspect(result.predict(2, 1.0) > 0.0, content="true")
  let total = result.predict(1, 4.0) +
    result.predict(2, 4.0) +
    result.survival_at(4.0)
  inspect(double_is_close(total, 1.0), content="true")
  inspect(competing_event_counts(data)[0], content="2")
  inspect(result.cause_hazard_at(1, 3.0) > 0.0, content="true")
}

///|
test "logrank and restricted mean compare groups" {
  let control = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  let treatment = SurvivalData::new([2.0, 4.0, 5.0, 6.0], [1, 1, 0, 1])
  let result = logrank_test(control, treatment)
  inspect(result.degrees_of_freedom, content="1")
  inspect(result.events, content="5")
  inspect(result.p_value >= 0.0 && result.p_value <= 1.0, content="true")
  inspect(result.z_score() != 0.0, content="true")
  inspect(
    restricted_mean_difference(control, treatment, 5.0) < 0.0,
    content="true",
  )
}

///|
test "pairwise logrank matrix has unit diagonal" {
  let groups = [
    SurvivalData::new([1.0, 2.0], [1, 0]),
    SurvivalData::new([2.0, 3.0], [1, 1]),
    SurvivalData::new([3.0, 4.0], [0, 1]),
  ]
  let matrix = pairwise_logrank(groups)
  inspect(matrix.length(), content="3")
  inspect(matrix[0][0], content="1")
  inspect(matrix[1][1], content="1")
  inspect(matrix[2][2], content="1")
}
