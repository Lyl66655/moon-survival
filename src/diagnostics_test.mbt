///|
test "cox residual diagnostics have one row per event" {
  let data = CoxData::new([1.0, 2.0, 3.0, 4.0, 5.0], [1, 1, 0, 1, 0], [
    [0.0],
    [0.5],
    [1.0],
    [1.5],
    [2.0],
  ])
  let fit = cox_fit(data)
  let schoenfeld = cox_schoenfeld_residuals(data, fit.coefficients)
  inspect(schoenfeld.length(), content="3")
  let martingale = cox_martingale_residuals(data, fit)
  inspect(martingale.length(), content="5")
  inspect(martingale[0].risk > 0.0, content="true")
  let diagnostics = cox_ph_diagnostics(data, fit)
  inspect(diagnostics.length(), content="1")
  inspect(diagnostics[0].covariate, content="0")
}

///|
test "calibration bins and discrimination metrics are bounded" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0, 5.0], [1, 0, 1, 0, 1])
  let survival = [0.1, 0.3, 0.5, 0.7, 0.9]
  let calibration = calibration_curve(data, survival, 3.0, 5)
  inspect(calibration.length(), content="5")
  let scores = [0.9, 0.7, 0.6, 0.3, 0.2]
  let concordance = concordance_index(data, scores)
  inspect(concordance.comparable > 0, content="true")
  inspect(
    concordance.c_index >= 0.0 && concordance.c_index <= 1.0,
    content="true",
  )
  let auc = time_dependent_auc(data, scores, 3.5)
  inspect(auc >= 0.0 && auc <= 1.0, content="true")
}

///|
test "diagnostics handle no comparable observations" {
  let data = SurvivalData::new([1.0, 2.0], [0, 0])
  let result = concordance_index(data, [0.5, 0.2])
  inspect(result.comparable, content="0")
  inspect(result.c_index, content="0.5")
  inspect(time_dependent_auc(data, [0.5, 0.2], 1.5), content="0.5")
}

///|
test "calibration output contains all requested bins" {
  let data = SurvivalData::new([1.0, 2.0, 4.0], [1, 1, 0])
  let curve = calibration_curve(data, [1.0, 0.5, 0.0], 3.0, 3)
  let mut count = 0
  for bin in curve {
    count += bin.count
  }
  inspect(count, content="3")
}
