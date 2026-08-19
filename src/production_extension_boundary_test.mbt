///|
/// Boundary coverage for the production extension modules.
test "parametric and AFT invalid inputs are rejected" {
  let empty = SurvivalData::new([], [])
  inspect(parametric_data_is_valid(empty), content="true")
  let aft_data = AftData::new([1.0, 2.0, 4.0], [1, 0, 1], [[0.0], [1.0], [2.0]])
  inspect(aft_data.length(), content="3")
  match fit_aft(aft_data, AftWeibull, max_iterations=0) {
    Err(InvalidIteration) => ()
    _ => fail("invalid AFT iteration count was accepted")
  }
  match AftData::try_new([-1.0], [1], [[]]) {
    Err(NegativeTime) => ()
    _ => fail("negative AFT time was accepted")
  }
}

///|
test "interval censoring handles empty and reversed intervals" {
  match IntervalSurvivalData::try_new([]) {
    Err(EmptyData) => ()
    _ => fail("empty interval data was accepted")
  }
  match interval_observation_is_valid(interval_event(4.0, 2.0)) {
    Err(ReversedInterval) => ()
    _ => fail("reversed interval was accepted")
  }
  let interval_data = IntervalSurvivalData::new([
    interval_exact(1.0),
    interval_event(2.0, 3.0),
    interval_right_censored(5.0),
  ])
  let fit = turnbull_fit(interval_data, max_iterations=3)
  inspect(turnbull_is_valid(fit), content="true")
  inspect(interval_kind_counts(interval_data).length(), content="4")
}

///|
test "IPCW remains bounded for all-events and all-censored data" {
  let all_events = SurvivalData::new([1.0, 2.0, 3.0], [1, 1, 1])
  let all_censored = SurvivalData::new([1.0, 2.0, 3.0], [0, 0, 0])
  let weights_a = ipcw_weights(all_events)
  let weights_b = ipcw_weights(all_censored)
  inspect(ipcw_is_valid(weights_a), content="true")
  inspect(ipcw_is_valid(weights_b), content="true")
  inspect(ipcw_risk_at(weights_a, 2.0) >= 0.0, content="true")
  inspect(ipcw_risk_at(weights_b, 2.0) <= 1.0, content="true")
}

///|
test "recurrent and multistate constructors reject structural boundaries" {
  match RecurrentData::try_new([], [], [], []) {
    Err(EmptyData) => ()
    _ => fail("empty recurrent data was accepted")
  }
  let recurrent = RecurrentData::try_new([0, 0], [0.0, 1.0], [1.0, 2.0], [1, 0])
  match recurrent {
    Ok(data) => inspect(recurrent_is_valid(data), content="true")
    Err(_) => fail("valid recurrent data was rejected")
  }
  match MultiStateData::try_new([0], [0], [0], [1.0], [0], 2) {
    Err(SelfTransition) => ()
    _ => fail("self transition was accepted")
  }
  let multistate = MultiStateData::try_new(
    [0, 0],
    [0, 1],
    [1, 2],
    [1.0, 2.0],
    [0],
    3,
  )
  match multistate {
    Ok(data) => inspect(multistate_is_valid(data), content="true")
    Err(_) => fail("valid multistate data was rejected")
  }
}

///|
test "validation, missing data and risk strata cover zero-size boundaries" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  inspect(validation_folds(data, 1).length(), content="0")
  inspect(validation_folds(data, 2).length(), content="2")
  let missing = missing_mean_impute([[1.0, 2.0], [3.0, 4.0]])
  inspect(imputation_is_valid(missing), content="true")
  let strata = survival_risk_stratify(
    [0.1, 0.3, 0.8, 0.9],
    [
      SurvivalData::new([1.0], [1]),
      SurvivalData::new([2.0], [0]),
      SurvivalData::new([3.0], [1]),
      SurvivalData::new([4.0], [0]),
    ],
    [0.5],
    3.0,
  )
  inspect(survival_risk_strata_is_valid(strata), content="true")
  inspect(survival_risk_worst_stratum(strata) is Some(_), content="true")
}

///|
test "time windows and administrative follow-up preserve invariants" {
  let data = SurvivalData::new([0.5, 1.5, 2.5, 4.0], [1, 0, 1, 1])
  let windows = survival_time_windows(data, 0.0, 4.0, 1.0)
  inspect(survival_time_windows_is_valid(windows), content="true")
  inspect(windows.length(), content="4")
  let followup = survival_apply_followup_horizon(data, 2.0)
  inspect(survival_followup_is_valid(followup), content="true")
  inspect(survival_followup_summary(followup, 2.0).sample_size, content="4")
  inspect(survival_followup_summary(followup, 2.0).truncated_rows, content="2")
}

///|
test "curve metrics and numerical primitives are stable at endpoints" {
  let points = survival_curve_metric_points([0.0, 1.0, 2.0], [1.0, 0.5, 0.0], [
    1.0, 0.4, 0.1,
  ])
  let summary = survival_curve_metric_summary(points)
  inspect(survival_curve_metric_points_are_valid(points), content="true")
  inspect(survival_curve_metric_summary_is_valid(summary), content="true")
  inspect(
    survival_numeric_expit(survival_numeric_logit(0.5)) > 0.49,
    content="true",
  )
  inspect(
    survival_numeric_expit(survival_numeric_logit(0.5)) < 0.51,
    content="true",
  )
  let fit = survival_numeric_linear_fit([0.0, 1.0, 2.0], [1.0, 3.0, 5.0])
  inspect(survival_numeric_linear_predict(fit, 3.0), content="7")
  let (_, _, solvable) = survival_numeric_solve_2x2(
    1.0, 0.0, 0.0, 1.0, 2.0, 3.0,
  )
  inspect(solvable, content="true")
}

///|
test "feature engineering and trial operations produce bounded tables" {
  let bins = survival_feature_bins([0.0, 1.0, 2.0, 3.0], [0, 1, 0, 1], 2)
  inspect(survival_feature_bins_are_valid(bins), content="true")
  inspect(survival_feature_one_hot([0, 2], [0, 1, 2]).length(), content="2")
  let targets = survival_trial_event_targets(100, 200, [0.5, 1.0], 0.8)
  inspect(survival_trial_targets_are_valid(targets), content="true")
  let updated = survival_trial_update_targets(targets, 100, 200)
  inspect(updated[1].achieved, content="true")
  let enrollment = survival_trial_enrollment([1.0, 1.0], [10, 10], [8, 12])
  inspect(survival_trial_enrollment_is_valid(enrollment), content="true")
}

///|
test "resampling, censoring and export remain valid for small cohorts" {
  let data = SurvivalData::new([1.0, 2.0, 3.0], [1, 0, 1])
  let resampled = survival_resample_with_replacement(data, 3, 42)
  inspect(resampled.length(), content="3")
  let summary = survival_resample_rmst_difference(data, data, 3.0, 4, 42)
  inspect(survival_resample_summary_is_valid(summary), content="true")
  let censoring = survival_censoring_summary(data)
  inspect(survival_censoring_summary_is_valid(censoring), content="true")
  let bundle = survival_export_bundle("local", "boundary", "0.2.1", data)
  inspect(survival_export_is_valid(bundle), content="true")
  inspect(bundle.artifact_names().length(), content="1")
}

///|
test "group KM, quality gates and patient timelines handle empty groups" {
  let data = SurvivalData::new([1.0, 2.0], [1, 0])
  let curves = survival_group_km(
    [data, SurvivalData::new([], [])],
    ["a", "b"],
    2.0,
  )
  inspect(survival_group_curves_are_valid(curves), content="true")
  inspect(survival_group_km_contrasts(curves).length(), content="1")
  let gate = survival_quality_gate_default()
  let rows = survival_quality_batch(["demo"], [data], gate)
  inspect(survival_quality_batch_is_valid(rows), content="true")
  let timeline = survival_patient_timeline(
    0,
    [0.0, 1.0],
    [1.0, 2.0],
    ["on", "on"],
    [0, 1],
  )
  inspect(survival_patient_timeline_is_valid(timeline), content="true")
  inspect(survival_timeline_merge(timeline.entries).length(), content="2")
}
