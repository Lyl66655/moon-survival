///|
test "cox data validation and stable likelihood" {
  let data = CoxData::new([3.0, 1.0, 2.0, 4.0], [1, 1, 0, 1], [
    [0.5],
    [1.0],
    [0.0],
    [-0.5],
  ])
  inspect(data.length(), content="4")
  inspect(data.get(0).time, content="1")
  let likelihood = cox_partial_likelihood(data, [0.0])
  inspect(likelihood > 0.0 && likelihood <= 1.0, content="true")
  let log_likelihood = cox_log_partial_likelihood(data, [0.0])
  inspect(log_likelihood < 0.0, content="true")
  let score = cox_score_vector(data, [0.0])
  inspect(score.length(), content="1")
  let information = cox_information_matrix(data, [0.0])
  inspect(information.length(), content="1")
  inspect(information[0][0] >= 0.0, content="true")
}

///|
test "cox newton fit returns diagnostics and predictions" {
  let simulated = simulate_cox(120, [0.6, -0.4], 0.08, 0.02, 88U)
  let fit = cox_fit(simulated.data, max_iterations=25)
  inspect(fit.coefficients.length(), content="2")
  inspect(fit.standard_errors.length(), content="2")
  inspect(fit.baseline.length() > 0, content="true")
  inspect(fit.iterations > 0, content="true")
  inspect(fit.log_likelihood >= fit.null_log_likelihood, content="true")
  inspect(fit.risk([1.0, 0.0]) > fit.risk([-1.0, 0.0]), content="true")
  inspect(fit.survival_at([0.0, 0.0], 2.0) <= 1.0, content="true")
  inspect(fit.baseline_survival(2.0) <= 1.0, content="true")
}

///|
test "cox singular systems are handled" {
  let singular = solve_linear_system([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0])
  match singular {
    None => ()
    Some(_) => fail("singular matrix produced a solution")
  }
  let identity = matrix_inverse([[1.0, 0.0], [0.0, 1.0]])
  match identity {
    Some(value) => inspect(value[1][1], content="1")
    None => fail("identity matrix was not invertible")
  }
  inspect(matrix_determinant([[2.0, 1.0], [1.0, 2.0]]), content="3")
}

///|
test "cox coefficient table exposes hazard ratios" {
  let data = CoxData::new([1.0, 2.0, 3.0, 4.0, 5.0], [1, 1, 0, 1, 0], [
    [0.0],
    [0.5],
    [1.0],
    [1.5],
    [2.0],
  ])
  let fit = cox_fit(data)
  inspect(fit.coefficients_table.length(), content="1")
  inspect(fit.coefficients_table[0].hazard_ratio > 0.0, content="true")
  inspect(
    fit.coefficients_table[0].upper_hr >= fit.coefficients_table[0].lower_hr,
    content="true",
  )
}
