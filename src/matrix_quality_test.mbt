///|
test "covariance and correlation matrices are symmetric" {
  let matrix = [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 5.0],
    [3.0, 6.0, 7.0],
    [4.0, 8.0, 9.0],
  ]
  let covariance = covariance_matrix(matrix)
  let correlation = correlation_matrix(matrix)
  inspect(covariance.length(), content="3")
  inspect(covariance[0].length(), content="3")
  inspect(double_is_close(covariance[0][1], covariance[1][0]), content="true")
  inspect(double_is_close(correlation[0][0], 1.0), content="true")
  inspect(double_is_close(correlation[0][1], correlation[1][0]), content="true")
  inspect(variance_inflation_factors(matrix).length(), content="3")
}

///|
test "design diagnostics expose collinearity" {
  let matrix = [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0]]
  let diagnostics = design_diagnostics(matrix)
  inspect(diagnostics.rows, content="4")
  inspect(diagnostics.columns, content="2")
  inspect(diagnostics.rank_proxy, content="2")
  inspect(diagnostics.condition_proxy, content="2")
  inspect(high_correlation_pairs(matrix).length(), content="1")
  inspect(high_correlation_pairs(matrix)[0].correlation > 0.99, content="true")
}

///|
test "matrix operations compose with solve and inverse" {
  let a = [[2.0, 1.0], [1.0, 2.0]]
  let b = [[1.0, 2.0], [3.0, 4.0]]
  let product = matrix_multiply(a, b)
  inspect(product[0][0], content="5")
  inspect(matrix_trace(a), content="4")
  let solution = solve_linear_system(a, [5.0, 4.0]).unwrap()
  inspect(solution[0], content="2")
  inspect(solution[1], content="1")
  let reconstructed = matrix_multiply(a, matrix_inverse(a).unwrap())
  inspect(double_is_close(reconstructed[0][0], 1.0), content="true")
  inspect(matrix_diagonal(a)[1], content="2")
}

///|
test "matrix utilities handle empty shapes" {
  inspect(matrix_transpose([]).length(), content="0")
  inspect(matrix_multiply([], []).length(), content="0")
  inspect(covariance_matrix([[1.0, 2.0]]).length(), content="0")
  inspect(correlation_matrix([]).length(), content="0")
}
