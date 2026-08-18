///|
test "clinical arm summaries expose censoring and horizon metrics" {
  let control = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  let treatment = SurvivalData::new([2.0, 3.0, 4.0, 5.0], [1, 0, 0, 1])
  let summary = summarize_trial_arm("control", control, 4.0)
  inspect(summary.name, content="control")
  inspect(summary.sample_size, content="4")
  inspect(summary.events, content="2")
  inspect(summary.censored, content="2")
  inspect(summary.event_fraction, content="0.5")
  inspect(summary.person_time, content="10")
  inspect(summary.survival_at_horizon < 1.0, content="true")
  let comparison = compare_survival_data(
    "control", control, "treatment", treatment, 4.0,
  )
  inspect(comparison.control.sample_size, content="4")
  inspect(comparison.treatment.sample_size, content="4")
  inspect(comparison.to_markdown().contains("Hazard ratio"), content="true")
}

///|
test "clinical effect transformations are bounded and interpretable" {
  inspect(absolute_risk_reduction(0.5, 0.3), content="0.2")
  inspect(number_needed_to_treat_from_arr(0.2), content="5")
  inspect(number_needed_to_treat_from_arr(0.0), content="0")
  inspect(absolute_risk_reduction_from_nnt(5.0), content="0.2")
  inspect(survival_risk_ratio(0.8, 0.9), content="0.5")
  inspect(survival_risk_ratio(1.0, 0.9), content="0")
  inspect(hazard_ratio_label(0.5), content="strongly favorable")
  inspect(hazard_ratio_label(0.9), content="favorable")
  inspect(hazard_ratio_label(1.0), content="near null")
  inspect(hazard_ratio_label(1.5), content="unfavorable")
  inspect(
    clinical_interpretation(0.7, 0.2, 0.01).contains("below 0.05"),
    content="true",
  )
}

///|
test "hazard ratio intervals and event rates remain finite" {
  let interval = hazard_ratio_interval(0.0, 0.1)
  inspect(interval.hazard_ratio, content="1")
  inspect(interval.lower < 1.0 && interval.upper > 1.0, content="true")
  let data = SurvivalData::new([1.0, 2.0, 3.0], [1, 0, 1])
  inspect(trial_event_rate(data) > 0.3, content="true")
  match
    cox_hazard_ratio_interval(
      {
        coefficients: [0.1],
        standard_errors: [0.2],
        covariance: [[0.04]],
        log_likelihood: -1.0,
        null_log_likelihood: -1.2,
        iterations: 1,
        converged: true,
        coefficients_table: [],
        baseline: [],
      },
      0,
    ) {
    Some(value) => inspect(value.hazard_ratio > 1.0, content="true")
    None => fail("valid Cox coefficient was not reported")
  }
  inspect(
    cox_hazard_ratio_interval(
      {
        coefficients: [0.1],
        standard_errors: [0.2],
        covariance: [[0.04]],
        log_likelihood: -1.0,
        null_log_likelihood: -1.2,
        iterations: 1,
        converged: true,
        coefficients_table: [],
        baseline: [],
      },
      2,
    )
    is None,
    content="true",
  )
}

///|
test "clinical reporting handles all-censored arms and invalid horizons" {
  let data = SurvivalData::new([2.0, 3.0], [0, 0])
  let summary = summarize_trial_arm("censored", data, -1.0)
  inspect(summary.events, content="0")
  inspect(summary.censored, content="2")
  inspect(summary.event_fraction, content="0")
  inspect(summary.person_time, content="5")
  inspect(summary.median_survival is None, content="true")
  inspect(trial_event_rate(data), content="0")
  inspect(number_needed_to_treat_from_arr(-0.1), content="0")
  inspect(absolute_risk_reduction_from_nnt(0.0), content="0")
  inspect(
    clinically_favorable({
      control: summary,
      treatment: summary,
      hazard_ratio: 0.5,
      hazard_ratio_lower: 0.2,
      hazard_ratio_upper: 0.8,
      logrank_p_value: 0.01,
      restricted_mean_difference: 0.0,
      absolute_risk_reduction: 0.1,
      number_needed_to_treat: 10.0,
      interpretation: "benefit",
    }),
    content="true",
  )
}
