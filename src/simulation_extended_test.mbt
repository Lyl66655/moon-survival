///|
test "deterministic generator is reproducible" {
  let first = DeterministicRng::new(42U)
  let second = DeterministicRng::new(42U)
  inspect(first.next_u32() == second.next_u32(), content="true")
  inspect(first.next_u32() == second.next_u32(), content="true")
  inspect(first.uniform() >= 0.0 && first.uniform() < 1.0, content="true")
}

///|
test "simulated survival fixture is deterministic and valid" {
  let config = default_simulation_config(80)
  let first = simulate_survival(config)
  let second = simulate_survival(config)
  inspect(first.data.length(), content="80")
  inspect(first.data == second.data, content="true")
  inspect(first.data.event_count() > 0, content="true")
  inspect(first.censoring_fraction() >= 0.0, content="true")
  inspect(first.diagnostics(20.0).sample_size, content="80")
}

///|
test "simulated cox data retains coefficient dimension" {
  let simulated = simulate_cox(40, [0.4, -0.2], 0.05, 0.01, 7U)
  inspect(simulated.data.length(), content="40")
  inspect(simulated.data.n_covariates, content="2")
  inspect(simulated.true_coefficients.length(), content="2")
  inspect(simulated.latent_failure_times.length(), content="40")
}

///|
test "bootstrap km returns bounded confidence bands" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0, 5.0], [1, 1, 0, 1, 0])
  let bootstrap = bootstrap_kaplan_meier(data, [2.0, 4.0], 20, 99U)
  inspect(bootstrap.replicates, content="20")
  inspect(bootstrap.mean.length(), content="2")
  inspect(bootstrap.lower[0] >= 0.0, content="true")
  inspect(bootstrap.upper[0] <= 1.0, content="true")
  inspect(bootstrap.upper[1] >= bootstrap.lower[1], content="true")
}

///|
test "math helpers cover numerical boundaries" {
  inspect(double_sqrt(-1.0), content="0")
  inspect(double_exp(0.0), content="1")
  inspect(double_exp(-1000.0), content="0")
  inspect(clamp_probability(-0.2), content="0")
  inspect(clamp_probability(1.2), content="1")
  inspect(logistic(0.0), content="0.5")
  inspect(double_is_close(1.0, 1.00000001), content="true")
  inspect(normal_cdf(0.0) > 0.49 && normal_cdf(0.0) < 0.51, content="true")
  inspect(normal_quantile(0.5).abs() < 0.0001, content="true")
}
