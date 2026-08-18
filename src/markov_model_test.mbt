///|
test "two state markov model preserves probability mass" {
  let model = two_state_survival_model(0.1, 1.0, 1.0, 0.0).unwrap()
  inspect(model.state_count(), content="2")
  inspect(model.state_index("dead") is Some(1), content="true")
  let distribution = model.distribution_after(0, 10)
  inspect(distribution.length(), content="2")
  inspect(distribution[0] + distribution[1] > 0.999999, content="true")
  inspect(distribution[0] < 0.4, content="true")
  inspect(model.absorption_probability(0, 100) > 0.9999, content="true")
}

///|
test "markov validation rejects malformed matrices" {
  match MarkovModel::try_new(["a"], [[0.5]], [0.0]) {
    Err(RowNotStochastic(0, _)) => ()
    _ => fail("non-stochastic row was accepted")
  }
  match MarkovModel::try_new(["a", "b"], [[1.0]], [0.0, 0.0]) {
    Err(TransitionRowCountMismatch) => ()
    _ => fail("short transition row was accepted")
  }
  match MarkovModel::try_new(["a"], [[-0.1]], [0.0]) {
    Err(InvalidProbability(0, 0, _)) => ()
    _ => fail("negative probability was accepted")
  }
  inspect(
    match markov_model_without_rewards([], []) {
      Err(_) => true
      Ok(_) => false
    },
    content="true",
  )
}

///|
test "markov trace supports rewards and serialization" {
  let model = two_state_survival_model(0.2, 1.0, 1.0, 0.0).unwrap()
  let trace = model.trace(0, 3, discount=0.9)
  inspect(trace.length(), content="4")
  inspect(trace[0].occupancy[0], content="1")
  inspect(trace[0].total_reward, content="1")
  inspect(trace[1].total_reward > 0.7, content="true")
  inspect(model.discounted_reward(0, 3, 0.9) > 2.0, content="true")
  inspect(markov_trace_to_csv(trace).contains("state_0"), content="true")
  inspect(markov_final_markdown(model, trace).contains("dead"), content="true")
}

///|
test "markov stationary distribution and hazard conversions are stable" {
  let model = MarkovModel::try_new(["a", "b"], [[0.8, 0.2], [0.4, 0.6]], [
    0.0, 0.0,
  ]).unwrap()
  let stationary = model.stationary_distribution()
  inspect(stationary.length(), content="2")
  inspect(stationary[0] > 0.6 && stationary[0] < 0.7, content="true")
  let probability = hazard_to_probability(0.1, 1.0)
  inspect(probability > 0.09 && probability < 0.11, content="true")
  inspect(probability_to_hazard(probability, 1.0) > 0.099, content="true")
}

///|
test "markov boundary inputs return safe empty or zero values" {
  let model = two_state_survival_model(0.1, 1.0, 1.0, 0.0).unwrap()
  inspect(model.distribution_after(-1, 1).length(), content="0")
  inspect(model.distribution_after(0, -1).length(), content="0")
  inspect(model.trace(0, 2, discount=-0.1).length(), content="0")
  inspect(model.trace_from([1.0], 2).length(), content="0")
  inspect(model.expected_time_to_absorption(-1), content="0")
  inspect(model.stationary_distribution(tolerance=0.0).length(), content="0")
  inspect(hazard_to_probability(-1.0, 1.0), content="0")
  inspect(probability_to_hazard(1.0, 1.0), content="0")
}
