///|
test "grouped survival summaries preserve cohort boundaries" {
  let groups = ["A", "B", "A", "B", "A", "B"]
  let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let events = [1, 1, 0, 1, 1, 0]
  let summaries = grouped_survival_summary(groups, times, events, 5.0)
  inspect(summaries.length(), content="2")
  inspect(summaries[0].group, content="A")
  inspect(summaries[0].sample_size, content="3")
  inspect(summaries[0].events, content="2")
  inspect(summaries[1].group, content="B")
  inspect(summaries[1].censoring_fraction >= 0.0, content="true")
  inspect(
    grouped_summary_to_csv(summaries).contains("group,sample_size"),
    content="true",
  )
}

///|
test "fixed effect meta analysis produces confidence limits" {
  let result = fixed_effect_meta_analysis([0.2, 0.4, 0.3], [0.1, 0.2, 0.15])
  inspect(result.studies, content="3")
  inspect(result.estimate > 0.2 && result.estimate < 0.4, content="true")
  inspect(result.standard_error > 0.0, content="true")
  inspect(result.lower < result.estimate, content="true")
  inspect(result.upper > result.estimate, content="true")
  inspect(result.p_value >= 0.0 && result.p_value <= 1.0, content="true")
  let hazard_ratio = meta_hazard_ratio(result)
  inspect(hazard_ratio.0 > 1.0, content="true")
  inspect(hazard_ratio.1 < hazard_ratio.0, content="true")
  inspect(hazard_ratio.2 > hazard_ratio.0, content="true")
}

///|
test "repeated groups are not duplicated" {
  let summaries = grouped_survival_summary(
    ["same", "same", "same"],
    [1.0, 2.0, 3.0],
    [1, 0, 1],
    3.0,
  )
  inspect(summaries.length(), content="1")
  inspect(summaries[0].events, content="2")
}
