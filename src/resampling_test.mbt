///|
test "bootstrap mean is reproducible and bounded" {
  let values = [1.0, 2.0, 3.0, 4.0, 5.0]
  let first = bootstrap_mean(values, 50, 123U)
  let second = bootstrap_mean(values, 50, 123U)
  inspect(first.estimate, content="3")
  inspect(first.replicates, content="50")
  inspect(first.standard_error >= 0.0, content="true")
  inspect(first.lower <= first.upper, content="true")
  inspect(first == second, content="true")
}

///|
test "jackknife and rmst intervals" {
  let values = [2.0, 4.0, 6.0, 8.0]
  let jackknife = jackknife_scalar(values, mean_doubles)
  inspect(jackknife.estimate, content="5")
  inspect(jackknife.replicates, content="4")
  inspect(jackknife.standard_error > 0.0, content="true")
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  let interval = bootstrap_rmst(data, 4.0, 20, 321U)
  inspect(interval.replicates, content="20")
  inspect(interval.lower <= interval.upper, content="true")
  let difference = bootstrap_mean_difference([1.0, 2.0], [3.0, 4.0], 20, 8U)
  inspect(difference.estimate, content="-2")
  inspect(difference.replicates, content="20")
}

///|
test "resampling empty inputs are safe" {
  let empty = bootstrap_mean([], 10, 1U)
  inspect(empty.replicates, content="0")
  let jackknife = jackknife_scalar([1.0], mean_doubles)
  inspect(jackknife.standard_error, content="0")
}
