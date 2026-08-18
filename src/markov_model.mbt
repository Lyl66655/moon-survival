///|
/// Errors returned while validating a finite-state Markov model.
pub(all) enum MarkovError {
  EmptyStateSpace
  StateNameMismatch
  TransitionRowCountMismatch
  TransitionColumnCountMismatch(Int)
  InvalidProbability(Int, Int, Double)
  RowNotStochastic(Int, Double)
  RewardLengthMismatch
  InvalidInitialState
  InvalidHorizon
  InvalidDiscount
} derive(Eq, Debug)

///|
/// A row-stochastic discrete-time Markov model with one reward per state.
/// Rows represent the current state and columns represent the next state.
pub(all) struct MarkovModel {
  states : Array[String]
  transition : Array[Array[Double]]
  rewards : Array[Double]
} derive(Eq, Debug)

///|
/// One row of a cohort simulation trace.
pub(all) struct MarkovTraceRow {
  period : Int
  occupancy : Array[Double]
  total_reward : Double
} derive(Eq, Debug)

///|
/// Build and validate a model. Every transition row must sum to one.
pub fn MarkovModel::try_new(
  states : Array[String],
  transition : Array[Array[Double]],
  rewards : Array[Double],
) -> Result[MarkovModel, MarkovError] {
  let state_count = states.length()
  if state_count == 0 {
    return Err(EmptyStateSpace)
  }
  if transition.length() != state_count {
    return Err(TransitionRowCountMismatch)
  }
  if rewards.length() != state_count {
    return Err(RewardLengthMismatch)
  }
  for i in 0..<state_count {
    if transition[i].length() != state_count {
      return Err(TransitionColumnCountMismatch(i))
    }
    let mut row_sum = 0.0
    for j in 0..<state_count {
      let probability = transition[i][j]
      if probability.is_nan() ||
        probability.is_inf() ||
        probability < 0.0 ||
        probability > 1.0 {
        return Err(InvalidProbability(i, j, probability))
      }
      row_sum += probability
    }
    if (row_sum - 1.0).abs() > 0.00000001 {
      return Err(RowNotStochastic(i, row_sum))
    }
  }
  for reward in rewards {
    if reward.is_nan() || reward.is_inf() {
      return Err(InvalidDiscount)
    }
  }
  Ok({ states, transition, rewards })
}

///|
/// Create a model with zero rewards for all states.
pub fn markov_model_without_rewards(
  states : Array[String],
  transition : Array[Array[Double]],
) -> Result[MarkovModel, MarkovError] {
  MarkovModel::try_new(states, transition, zeros(states.length()))
}

///|
/// Number of states in the model.
pub fn MarkovModel::state_count(self : MarkovModel) -> Int {
  self.states.length()
}

///|
/// Return a copy of the state labels.
pub fn MarkovModel::state_names(self : MarkovModel) -> Array[String] {
  let result : Array[String] = []
  for state in self.states {
    result.push(state)
  }
  result
}

///|
/// Find a state index by its stable label.
pub fn MarkovModel::state_index(self : MarkovModel, name : String) -> Int? {
  for i in 0..<self.states.length() {
    if self.states[i] == name {
      return Some(i)
    }
  }
  None
}

///|
/// Return a defensive copy of the transition matrix.
pub fn MarkovModel::transition_matrix(
  self : MarkovModel,
) -> Array[Array[Double]] {
  matrix_copy(self.transition)
}

///|
/// Validate an initial distribution for a cohort simulation.
pub fn markov_validate_distribution(distribution : Array[Double]) -> Bool {
  if distribution.length() == 0 {
    return false
  }
  let mut total = 0.0
  for probability in distribution {
    if probability.is_nan() || probability.is_inf() || probability < 0.0 {
      return false
    }
    total += probability
  }
  (total - 1.0).abs() <= 0.00000001
}

///|
/// Propagate a distribution by one period.
fn markov_step(
  transition : Array[Array[Double]],
  distribution : Array[Double],
) -> Array[Double] {
  let next = zeros(distribution.length())
  for i in 0..<distribution.length() {
    for j in 0..<distribution.length() {
      next[j] += distribution[i] * transition[i][j]
    }
  }
  next
}

///|
/// Propagate a point mass from one state for a finite number of periods.
/// Invalid inputs return an empty vector so callers can handle them without a
/// panic in production data pipelines.
pub fn MarkovModel::distribution_after(
  self : MarkovModel,
  initial_state : Int,
  periods : Int,
) -> Array[Double] {
  if initial_state < 0 || initial_state >= self.state_count() || periods < 0 {
    return []
  }
  let distribution = zeros(self.state_count())
  distribution[initial_state] = 1.0
  markov_distribution_after(self, distribution, periods)
}

///|
/// Propagate an arbitrary initial distribution.
pub fn MarkovModel::distribution_from(
  self : MarkovModel,
  initial : Array[Double],
  periods : Int,
) -> Array[Double] {
  if initial.length() != self.state_count() ||
    periods < 0 ||
    !markov_validate_distribution(initial) {
    return []
  }
  markov_distribution_after(self, initial, periods)
}

///|
fn markov_distribution_after(
  model : MarkovModel,
  initial : Array[Double],
  periods : Int,
) -> Array[Double] {
  let distribution = []
  for probability in initial {
    distribution.push(probability)
  }
  let mut current = distribution
  for _ in 0..<periods {
    current = markov_step(model.transition, current)
  }
  current
}

///|
/// Produce period zero through period `periods` for a point-mass cohort.
pub fn MarkovModel::trace(
  self : MarkovModel,
  initial_state : Int,
  periods : Int,
  discount? : Double = 1.0,
) -> Array[MarkovTraceRow] {
  if initial_state < 0 ||
    initial_state >= self.state_count() ||
    periods < 0 ||
    discount < 0.0 ||
    discount > 1.0 {
    return []
  }
  let current = zeros(self.state_count())
  current[initial_state] = 1.0
  markov_trace(self, current, periods, discount)
}

///|
/// Produce a trace from an arbitrary initial distribution.
pub fn MarkovModel::trace_from(
  self : MarkovModel,
  initial : Array[Double],
  periods : Int,
  discount? : Double = 1.0,
) -> Array[MarkovTraceRow] {
  if initial.length() != self.state_count() ||
    periods < 0 ||
    discount < 0.0 ||
    discount > 1.0 ||
    !markov_validate_distribution(initial) {
    return []
  }
  markov_trace(self, initial, periods, discount)
}

///|
fn markov_trace(
  model : MarkovModel,
  initial : Array[Double],
  periods : Int,
  discount : Double,
) -> Array[MarkovTraceRow] {
  let current = []
  for probability in initial {
    current.push(probability)
  }
  let result : Array[MarkovTraceRow] = []
  let mut distribution = current
  let mut discount_factor = 1.0
  for period in 0..<(periods + 1) {
    let total_reward = dot_product(distribution, model.rewards) *
      discount_factor
    let occupancy = []
    for probability in distribution {
      occupancy.push(probability)
    }
    result.push({ period, occupancy, total_reward })
    distribution = markov_step(model.transition, distribution)
    discount_factor *= discount
  }
  result
}

///|
/// Sum expected state rewards over the first `periods` transitions.
pub fn MarkovModel::discounted_reward(
  self : MarkovModel,
  initial_state : Int,
  periods : Int,
  discount : Double,
) -> Double {
  if initial_state < 0 ||
    initial_state >= self.state_count() ||
    periods <= 0 ||
    discount < 0.0 ||
    discount > 1.0 {
    return 0.0
  }
  let trace = self.trace(initial_state, periods - 1, discount~)
  let mut total = 0.0
  for row in trace {
    total += row.total_reward
  }
  total
}

///|
/// Identify exactly absorbing states.
pub fn MarkovModel::absorbing_states(self : MarkovModel) -> Array[Int] {
  let result : Array[Int] = []
  for i in 0..<self.state_count() {
    let absorbing = self.transition[i][i] == 1.0
    let mut no_exit = true
    for j in 0..<self.state_count() {
      if j != i && self.transition[i][j] > 0.0000000001 {
        no_exit = false
      }
    }
    if absorbing && no_exit {
      result.push(i)
    }
  }
  result
}

///|
/// Probability mass in absorbing states after a finite horizon.
pub fn MarkovModel::absorption_probability(
  self : MarkovModel,
  initial_state : Int,
  periods : Int,
) -> Double {
  let distribution = self.distribution_after(initial_state, periods)
  if distribution.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for index in self.absorbing_states() {
    total += distribution[index]
  }
  total
}

///|
/// Estimate the expected time to absorption by iterating the cohort mass.
/// The result is capped at `max_periods` for non-absorbing chains.
pub fn MarkovModel::expected_time_to_absorption(
  self : MarkovModel,
  initial_state : Int,
  max_periods? : Int = 10000,
) -> Double {
  if initial_state < 0 ||
    initial_state >= self.state_count() ||
    max_periods <= 0 {
    return 0.0
  }
  let absorbing = self.absorbing_states()
  let mut distribution = self.distribution_after(initial_state, 0)
  let mut expected = 0.0
  let mut mut_period = 0
  while mut_period < max_periods {
    let mut nonabsorbed = 1.0
    for index in absorbing {
      nonabsorbed -= distribution[index]
    }
    if nonabsorbed <= 0.0000000001 {
      return expected
    }
    expected += nonabsorbed
    distribution = markov_step(self.transition, distribution)
    mut_period += 1
  }
  expected
}

///|
/// Find the stationary distribution with power iteration.
pub fn MarkovModel::stationary_distribution(
  self : MarkovModel,
  tolerance? : Double = 0.000000001,
  max_iterations? : Int = 10000,
) -> Array[Double] {
  if tolerance <= 0.0 || max_iterations <= 0 {
    return []
  }
  let initial = zeros(self.state_count())
  initial[0] = 1.0
  let mut current = initial
  let mut iteration = 0
  while iteration < max_iterations {
    let next = markov_step(self.transition, current)
    if vector_norm(vector_subtract(next, current)) <= tolerance {
      return next
    }
    current = next
    iteration += 1
  }
  current
}

///|
/// Convert a constant hazard over one period into an event probability.
pub fn hazard_to_probability(hazard : Double, period : Double) -> Double {
  if hazard <= 0.0 || period <= 0.0 || hazard.is_nan() || period.is_nan() {
    0.0
  } else {
    clamp_probability(1.0 - double_exp(-hazard * period))
  }
}

///|
/// Convert a one-period event probability into a continuous-time hazard.
pub fn probability_to_hazard(probability : Double, period : Double) -> Double {
  if probability <= 0.0 ||
    period <= 0.0 ||
    probability >= 1.0 ||
    period.is_nan() {
    0.0
  } else {
    -double_log(1.0 - probability) / period
  }
}

///|
/// Construct a two-state alive/dead model from a constant hazard.
pub fn two_state_survival_model(
  hazard : Double,
  period : Double,
  reward_alive : Double,
  reward_dead : Double,
) -> Result[MarkovModel, MarkovError] {
  let event_probability = hazard_to_probability(hazard, period)
  MarkovModel::try_new(
    ["alive", "dead"],
    [[1.0 - event_probability, event_probability], [0.0, 1.0]],
    [reward_alive, reward_dead],
  )
}

///|
/// Render a trace as CSV for audit trails and downstream plotting.
pub fn markov_trace_to_csv(trace : Array[MarkovTraceRow]) -> String {
  if trace.length() == 0 {
    return "period,total_reward\n"
  }
  let builder = StringBuilder::new()
  builder.write_string("period,total_reward")
  for i in 0..<trace[0].occupancy.length() {
    builder.write_string(",state_\{i}")
  }
  builder.write_char('\n')
  for row in trace {
    builder.write_string("\{row.period},\{row.total_reward}")
    for probability in row.occupancy {
      builder.write_string(",\{probability}")
    }
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Render the final cohort state as a compact Markdown table.
pub fn markov_final_markdown(
  model : MarkovModel,
  trace : Array[MarkovTraceRow],
) -> String {
  if trace.length() == 0 {
    return "| State | Probability |\n|---|---:|\n"
  }
  let final_row = trace[trace.length() - 1]
  let builder = StringBuilder::new()
  builder.write_string("| State | Probability |\n|---|---:|\n")
  for i in 0..<model.state_count() {
    builder.write_string("| \{model.states[i]} | \{final_row.occupancy[i]} |\n")
  }
  builder.to_string()
}
