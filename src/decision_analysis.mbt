///|
/// Health-state rewards for a cohort cost-effectiveness model.
pub(all) struct EconomicStrategy {
  name : String
  states : Array[String]
  transition : Array[Array[Double]]
  costs : Array[Double]
  utilities : Array[Double]
  initial : Array[Double]
  strategy_row : Bool
} derive(Eq, Debug)

///|
/// Validation failures for economic decision models.
pub(all) enum EconomicError {
  EmptyStates
  LengthMismatch
  InvalidProbability
  InvalidDistribution
  InvalidDiscount
  InvalidHorizon
  InvalidReward
} derive(Eq, Debug)

///|
/// Discounted outcome from a cohort simulation.
pub(all) struct EconomicOutcome {
  strategy : String
  total_cost : Double
  total_qaly : Double
  discounted_cost : Double
  discounted_qaly : Double
  life_years : Double
  cycles : Int
  outcome_row : Bool
} derive(Eq, Debug)

///|
/// Incremental economic comparison between two strategies.
pub(all) struct EconomicComparison {
  strategy : String
  comparator : String
  incremental_cost : Double
  incremental_qaly : Double
  icer : Double
  net_monetary_benefit : Double
  dominance : String
} derive(Eq, Debug)

///|
/// A deterministic probabilistic-sensitivity-analysis summary.
pub(all) struct EconomicPsaSummary {
  strategy : String
  simulations : Int
  mean_cost : Double
  mean_qaly : Double
  probability_cost_effective : Double
  mean_net_benefit : Double
  wtp : Double
} derive(Eq, Debug)

///|
/// Validate a strategy definition.
pub fn EconomicStrategy::try_new(
  name : String,
  states : Array[String],
  transition : Array[Array[Double]],
  costs : Array[Double],
  utilities : Array[Double],
  initial : Array[Double],
) -> Result[EconomicStrategy, EconomicError] {
  let state_count = states.length()
  if state_count == 0 {
    return Err(EmptyStates)
  }
  if transition.length() != state_count ||
    costs.length() != state_count ||
    utilities.length() != state_count ||
    initial.length() != state_count {
    return Err(LengthMismatch)
  }
  for i in 0..<state_count {
    if transition[i].length() != state_count {
      return Err(LengthMismatch)
    }
    let mut row_total = 0.0
    for probability in transition[i] {
      if probability.is_nan() ||
        probability.is_inf() ||
        probability < 0.0 ||
        probability > 1.0 {
        return Err(InvalidProbability)
      }
      row_total += probability
    }
    if (row_total - 1.0).abs() > 0.000001 {
      return Err(InvalidProbability)
    }
    if costs[i].is_nan() ||
      costs[i].is_inf() ||
      utilities[i].is_nan() ||
      utilities[i].is_inf() {
      return Err(InvalidReward)
    }
  }
  let mut initial_total = 0.0
  for probability in initial {
    if probability < 0.0 || probability.is_nan() || probability.is_inf() {
      return Err(InvalidDistribution)
    }
    initial_total += probability
  }
  if (initial_total - 1.0).abs() > 0.000001 {
    return Err(InvalidDistribution)
  }
  Ok({ name, states, transition, costs, utilities, initial, strategy_row: true })
}

///|
/// Construct a validated strategy, panicking only on programmer errors.
pub fn EconomicStrategy::new(
  name : String,
  states : Array[String],
  transition : Array[Array[Double]],
  costs : Array[Double],
  utilities : Array[Double],
  initial : Array[Double],
) -> EconomicStrategy {
  match
    EconomicStrategy::try_new(
      name, states, transition, costs, utilities, initial,
    ) {
    Ok(value) => value
    Err(_) => panic()
  }
}

///|
/// Multiply a row vector by a transition matrix.
fn economic_step(
  matrix : Array[Array[Double]],
  distribution : Array[Double],
) -> Array[Double] {
  let next = Array::make(distribution.length(), 0.0)
  for from_state in 0..<distribution.length() {
    for to_state in 0..<distribution.length() {
      next[to_state] += distribution[from_state] * matrix[from_state][to_state]
    }
  }
  next
}

///|
/// Present value factor for a cycle index.
pub fn economic_discount_factor(rate : Double, cycle : Int) -> Double {
  if rate < 0.0 || cycle < 0 {
    1.0
  } else {
    1.0 / double_exp(cycle.to_double() * double_log(1.0 + rate))
  }
}

///|
/// Simulate discounted cost and QALYs for one strategy.
pub fn economic_simulate(
  strategy : EconomicStrategy,
  cycles : Int,
  discount_rate : Double,
  cycle_length : Double,
) -> EconomicOutcome {
  if cycles <= 0 || cycle_length <= 0.0 || discount_rate < 0.0 {
    return {
      strategy: strategy.name,
      total_cost: 0.0,
      total_qaly: 0.0,
      discounted_cost: 0.0,
      discounted_qaly: 0.0,
      life_years: 0.0,
      cycles: 0,
      outcome_row: true,
    }
  }
  let mut distribution = strategy.initial
  let mut total_cost = 0.0
  let mut total_qaly = 0.0
  let mut discounted_cost = 0.0
  let mut discounted_qaly = 0.0
  for cycle in 0..<cycles {
    let factor = economic_discount_factor(discount_rate, cycle)
    let cycle_cost = dot_product(distribution, strategy.costs) * cycle_length
    let cycle_qaly = dot_product(distribution, strategy.utilities) *
      cycle_length
    total_cost += cycle_cost
    total_qaly += cycle_qaly
    discounted_cost += cycle_cost * factor
    discounted_qaly += cycle_qaly * factor
    distribution = economic_step(strategy.transition, distribution)
  }
  {
    strategy: strategy.name,
    total_cost,
    total_qaly,
    discounted_cost,
    discounted_qaly,
    life_years: total_qaly,
    cycles,
    outcome_row: true,
  }
}

///|
/// Simulate a collection of strategies.
pub fn economic_simulate_all(
  strategies : Array[EconomicStrategy],
  cycles : Int,
  discount_rate : Double,
  cycle_length : Double,
) -> Array[EconomicOutcome] {
  let result : Array[EconomicOutcome] = []
  for strategy in strategies {
    result.push(
      economic_simulate(strategy, cycles, discount_rate, cycle_length),
    )
  }
  result
}

///|
/// Compare one outcome against a comparator outcome.
pub fn economic_compare(
  outcome : EconomicOutcome,
  comparator : EconomicOutcome,
  wtp : Double,
) -> EconomicComparison {
  let incremental_cost = outcome.discounted_cost - comparator.discounted_cost
  let incremental_qaly = outcome.discounted_qaly - comparator.discounted_qaly
  let icer = if incremental_qaly.abs() <= survival_epsilon {
    0.0
  } else {
    incremental_cost / incremental_qaly
  }
  let nmb = incremental_qaly * wtp - incremental_cost
  let dominance = if incremental_cost <= 0.0 && incremental_qaly >= 0.0 {
    "dominant"
  } else if incremental_cost >= 0.0 && incremental_qaly <= 0.0 {
    "dominated"
  } else {
    "tradeoff"
  }
  {
    strategy: outcome.strategy,
    comparator: comparator.strategy,
    incremental_cost,
    incremental_qaly,
    icer,
    net_monetary_benefit: nmb,
    dominance,
  }
}

///|
/// Sort outcomes by discounted cost.
pub fn economic_sort_outcomes(
  outcomes : Array[EconomicOutcome],
) -> Array[EconomicOutcome] {
  let result : Array[EconomicOutcome] = []
  for outcome in outcomes {
    result.push(outcome)
  }
  result.sort_by(fn(a, b) {
    if a.discounted_cost < b.discounted_cost {
      -1
    } else if a.discounted_cost > b.discounted_cost {
      1
    } else {
      0
    }
  })
  result
}

///|
/// Compute net monetary benefit.
pub fn economic_net_benefit(outcome : EconomicOutcome, wtp : Double) -> Double {
  outcome.discounted_qaly * wtp - outcome.discounted_cost
}

///|
/// Return the strategy with the maximum net monetary benefit.
pub fn economic_best_strategy(
  outcomes : Array[EconomicOutcome],
  wtp : Double,
) -> EconomicOutcome? {
  let mut best : EconomicOutcome? = None
  for outcome in outcomes {
    match best {
      None => best = Some(outcome)
      Some(current) =>
        if economic_net_benefit(outcome, wtp) >
          economic_net_benefit(current, wtp) {
          best = Some(outcome)
        }
    }
  }
  best
}

///|
/// Produce a simple deterministic parameter perturbation for PSA.
fn economic_perturb(value : Double, index : Int) -> Double {
  let sign = if index % 2 == 0 { 1.0 } else { -1.0 }
  value * (1.0 + sign * 0.01 * (index % 7 + 1).to_double())
}

///|
/// Run a reproducible probabilistic sensitivity analysis.
pub fn economic_psa(
  strategy : EconomicStrategy,
  cycles : Int,
  discount_rate : Double,
  cycle_length : Double,
  simulations : Int,
  wtp : Double,
) -> EconomicPsaSummary {
  if simulations <= 0 {
    return {
      strategy: strategy.name,
      simulations: 0,
      mean_cost: 0.0,
      mean_qaly: 0.0,
      probability_cost_effective: 0.0,
      mean_net_benefit: 0.0,
      wtp,
    }
  }
  let mut cost_total = 0.0
  let mut qaly_total = 0.0
  let mut nmb_total = 0.0
  let mut cost_effective = 0
  for simulation in 0..<simulations {
    let transition : Array[Array[Double]] = []
    for row_index in 0..<strategy.transition.length() {
      let row : Array[Double] = []
      for column_index in 0..<strategy.transition[row_index].length() {
        row.push(strategy.transition[row_index][column_index])
      }
      transition.push(row)
    }
    let costs : Array[Double] = []
    let utilities : Array[Double] = []
    for index in 0..<strategy.costs.length() {
      costs.push(economic_perturb(strategy.costs[index], simulation + index))
      utilities.push(
        clamp_probability(
          economic_perturb(strategy.utilities[index], simulation + index),
        ),
      )
    }
    let perturbed = EconomicStrategy::new(
      strategy.name,
      strategy.states,
      transition,
      costs,
      utilities,
      strategy.initial,
    )
    let outcome = economic_simulate(
      perturbed, cycles, discount_rate, cycle_length,
    )
    cost_total += outcome.discounted_cost
    qaly_total += outcome.discounted_qaly
    let nmb = economic_net_benefit(outcome, wtp)
    nmb_total += nmb
    if nmb >= 0.0 {
      cost_effective += 1
    }
  }
  {
    strategy: strategy.name,
    simulations,
    mean_cost: cost_total / simulations.to_double(),
    mean_qaly: qaly_total / simulations.to_double(),
    probability_cost_effective: cost_effective.to_double() /
    simulations.to_double(),
    mean_net_benefit: nmb_total / simulations.to_double(),
    wtp,
  }
}

///|
/// Check a strategy's transition and reward invariants.
pub fn economic_strategy_is_valid(strategy : EconomicStrategy) -> Bool {
  match
    EconomicStrategy::try_new(
      strategy.name,
      strategy.states,
      strategy.transition,
      strategy.costs,
      strategy.utilities,
      strategy.initial,
    ) {
    Ok(_) => true
    Err(_) => false
  }
}

///|
/// Export cohort outcomes.
pub fn economic_outcomes_to_csv(outcomes : Array[EconomicOutcome]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "strategy,total_cost,total_qaly,discounted_cost,discounted_qaly,cycles\n",
  )
  for outcome in outcomes {
    builder.write_string(outcome.strategy)
    builder.write_char(',')
    builder.write_string(outcome.total_cost.to_string())
    builder.write_char(',')
    builder.write_string(outcome.total_qaly.to_string())
    builder.write_char(',')
    builder.write_string(outcome.discounted_cost.to_string())
    builder.write_char(',')
    builder.write_string(outcome.discounted_qaly.to_string())
    builder.write_char(',')
    builder.write_string(outcome.cycles.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export incremental comparisons.
pub fn economic_comparisons_to_csv(
  comparisons : Array[EconomicComparison],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "strategy,comparator,incremental_cost,incremental_qaly,icer,net_monetary_benefit,dominance\n",
  )
  for comparison in comparisons {
    builder.write_string(comparison.strategy)
    builder.write_char(',')
    builder.write_string(comparison.comparator)
    builder.write_char(',')
    builder.write_string(comparison.incremental_cost.to_string())
    builder.write_char(',')
    builder.write_string(comparison.incremental_qaly.to_string())
    builder.write_char(',')
    builder.write_string(comparison.icer.to_string())
    builder.write_char(',')
    builder.write_string(comparison.net_monetary_benefit.to_string())
    builder.write_char(',')
    builder.write_string(comparison.dominance)
    builder.write_char('\n')
  }
  builder.to_string()
}
