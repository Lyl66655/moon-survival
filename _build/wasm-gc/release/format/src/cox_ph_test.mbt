///|
test "CoxData::new and partial_likelihood" {
  let times = [1.0, 2.0, 3.0]
  let events = [1, 0, 1]
  let covs = [[0.5], [1.0], [1.5]]
  let data = CoxData::new(times, events, covs)

  inspect(data.n_covariates, content="1")

  // Test partial likelihood with beta = [0.0]
  // With beta=0, exp(bx) = 1.
  // At t=1 (event): risk set size 3. partial lik = 1/3
  // At t=2 (censored): ignored
  // At t=3 (event): risk set size 1. partial lik = 1/1
  // Total = 1/3
  let pl = cox_partial_likelihood(data, [0.0])
  let diff = pl - 0.333333333
  // Using an epsilon for float comparison
  let close_enough = diff < 0.0001 && diff > -0.0001
  inspect(close_enough, content="true")
}
