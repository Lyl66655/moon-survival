///|
test "vector and matrix helpers reject shape mismatches" {
  inspect(vector_add([1.0, 2.0], [3.0, 4.0])[1], content="6")
  inspect(vector_scale([1.0, 2.0], 2.0)[0], content="2")
  inspect(dot_product([1.0, 2.0], [3.0, 4.0]), content="11")
  inspect(vector_norm([3.0, 4.0]), content="5")
  inspect(identity_matrix(2)[1][1], content="1")
  inspect(matrix_to_string([[1.0, 2.0]]).contains("1,2"), content="true")
}

///|
test "numeric interpolation and distribution approximations" {
  inspect(linspace(0.0, 1.0, 1)[0], content="0")
  inspect(linspace(0.0, 1.0, 3)[1], content="0.5")
  inspect(mean_doubles([2.0, 4.0]), content="3")
  inspect(weighted_mean([1.0, 3.0], [1.0, 3.0]), content="2.5")
  inspect(
    normal_density(0.0) > 0.39 && normal_density(0.0) < 0.41,
    content="true",
  )
  inspect(double_is_close(logistic(logit(0.7)), 0.7), content="true")
  inspect(log_sum_exp([0.0, 0.0]) > 0.69, content="true")
}

///|
test "csv quote round trips preserve punctuation" {
  let table : CsvTable = {
    headers: ["id", "note"],
    rows: [["1", "hello, \"moon\""]],
  }
  let encoded = csv_stringify(table)
  let decoded = csv_parse(encoded).unwrap()
  inspect(decoded.rows[0][1], content="hello, \"moon\"")
  inspect(csv_escape("plain"), content="plain")
}

///|
test "risk and power calculations are stable for neutral effects" {
  inspect(logrank_required_events(1.0, 0.05, 0.8) > 0, content="true")
  inspect(logrank_power(100, 1.0, 0.05) < 0.6, content="true")
  inspect(detectable_hazard_ratio(1000, 0.05, 0.8) > 1.0, content="true")
  let score = model_score("null", 0, 100, 0.0)
  inspect(score.aic, content="0")
  inspect(score.bic, content="0")
}
