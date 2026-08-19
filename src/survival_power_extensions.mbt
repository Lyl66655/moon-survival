///|
/// A sample-size scenario for a survival endpoint.
pub(all) struct SurvivalPowerScenario {
  alpha : Double
  power : Double
  hazard_ratio : Double
  allocation : Double
  expected_events : Int
  total_sample_size : Int
  event_fraction : Double
  power_row : Bool
} derive(Eq, Debug)

///|
/// An achieved-power calculation.
pub(all) struct SurvivalAchievedPower {
  events : Int
  hazard_ratio : Double
  alpha : Double
  noncentrality : Double
  power : Double
  achieved_row : Bool
} derive(Eq, Debug)

///|
/// Convert a log-rank information target to event count.
pub fn survival_power_event_target(
  alpha : Double,
  power : Double,
  hazard_ratio : Double,
  allocation : Double,
) -> Int {
  if hazard_ratio <= 0.0 || allocation <= 0.0 || allocation >= 1.0 {
    return 0
  }
  let z_alpha = normal_quantile(1.0 - clamp_probability(alpha) / 2.0)
  let z_power = normal_quantile(clamp_probability(power))
  let effect = double_log(hazard_ratio).abs()
  if effect <= survival_epsilon {
    return 0
  }
  ((z_alpha + z_power) *
  (z_alpha + z_power) /
  (effect * effect * allocation * (1.0 - allocation)))
  .ceil()
  .to_int()
}

///|
/// Convert event target to total sample size using an event fraction.
pub fn survival_power_sample_target(
  events : Int,
  event_fraction : Double,
) -> Int {
  if events <= 0 || event_fraction <= 0.0 {
    0
  } else {
    (events.to_double() / clamp_probability(event_fraction)).ceil().to_int()
  }
}

///|
/// Build a planning scenario.
pub fn survival_power_scenario(
  alpha : Double,
  power : Double,
  hazard_ratio : Double,
  allocation : Double,
  event_fraction : Double,
) -> SurvivalPowerScenario {
  let events = survival_power_event_target(
    alpha, power, hazard_ratio, allocation,
  )
  {
    alpha: clamp_probability(alpha),
    power: clamp_probability(power),
    hazard_ratio,
    allocation,
    expected_events: events,
    total_sample_size: survival_power_sample_target(events, event_fraction),
    event_fraction: clamp_probability(event_fraction),
    power_row: true,
  }
}

///|
/// Compute achieved power for an observed event count.
pub fn survival_achieved_power(
  events : Int,
  hazard_ratio : Double,
  alpha : Double,
  allocation : Double,
) -> SurvivalAchievedPower {
  let effect = double_log(hazard_ratio).abs()
  let information = if events <= 0 || allocation <= 0.0 || allocation >= 1.0 {
    0.0
  } else {
    double_sqrt(events.to_double() * allocation * (1.0 - allocation)) * effect
  }
  let critical = normal_quantile(1.0 - clamp_probability(alpha) / 2.0)
  {
    events,
    hazard_ratio,
    alpha,
    noncentrality: information,
    power: clamp_probability(
      1.0 -
      normal_cdf(critical - information) +
      normal_cdf(-critical - information),
    ),
    achieved_row: true,
  }
}

///|
/// Compute power over a hazard-ratio grid.
pub fn survival_power_grid(
  alpha : Double,
  power : Double,
  hazard_ratios : Array[Double],
  allocation : Double,
  event_fraction : Double,
) -> Array[SurvivalPowerScenario] {
  let result : Array[SurvivalPowerScenario] = []
  for hazard_ratio in hazard_ratios {
    result.push(
      survival_power_scenario(
        alpha, power, hazard_ratio, allocation, event_fraction,
      ),
    )
  }
  result
}

///|
/// Compute a detectable hazard ratio for an event target.
pub fn survival_detectable_hazard_ratio(
  alpha : Double,
  power : Double,
  events : Int,
  allocation : Double,
) -> Double {
  if events <= 0 || allocation <= 0.0 || allocation >= 1.0 {
    return 1.0
  }
  let z_alpha = normal_quantile(1.0 - clamp_probability(alpha) / 2.0)
  let z_power = normal_quantile(clamp_probability(power))
  double_exp(
    (z_alpha + z_power) /
    double_sqrt(events.to_double() * allocation * (1.0 - allocation)),
  )
}

///|
/// Compute a recruitment time from sample size and accrual rate.
pub fn survival_recruitment_time(
  sample_size : Int,
  accrual_rate : Double,
  dropout_fraction : Double,
) -> Double {
  if sample_size <= 0 || accrual_rate <= 0.0 {
    0.0
  } else {
    sample_size.to_double() *
    (1.0 + clamp_probability(dropout_fraction)) /
    accrual_rate
  }
}

///|
/// Compute an event fraction from a simulated dataset.
pub fn survival_event_fraction(data : SurvivalData) -> Double {
  if data.length() == 0 {
    0.0
  } else {
    data.event_count().to_double() / data.length().to_double()
  }
}

///|
/// Compute a design's expected event count.
pub fn survival_expected_events(
  sample_size : Int,
  event_fraction : Double,
) -> Int {
  (sample_size.to_double() * clamp_probability(event_fraction)).round().to_int()
}

///|
/// Compute power under unequal allocation.
pub fn survival_power_with_allocation(
  events : Int,
  hazard_ratio : Double,
  alpha : Double,
  allocation : Double,
) -> Double {
  survival_achieved_power(events, hazard_ratio, alpha, allocation).power
}

///|
/// Compute the minimum power over a scenario set.
pub fn survival_power_minimum(
  scenarios : Array[SurvivalPowerScenario],
) -> Double {
  let values : Array[Double] = []
  for scenario in scenarios {
    values.push(scenario.power)
  }
  survival_ext_min(values)
}

///|
/// Export power scenarios.
pub fn survival_power_scenarios_to_csv(
  scenarios : Array[SurvivalPowerScenario],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "alpha,power,hazard_ratio,allocation,expected_events,total_sample_size,event_fraction\n",
  )
  for scenario in scenarios {
    builder.write_string(scenario.alpha.to_string())
    builder.write_char(',')
    builder.write_string(scenario.power.to_string())
    builder.write_char(',')
    builder.write_string(scenario.hazard_ratio.to_string())
    builder.write_char(',')
    builder.write_string(scenario.allocation.to_string())
    builder.write_char(',')
    builder.write_string(scenario.expected_events.to_string())
    builder.write_char(',')
    builder.write_string(scenario.total_sample_size.to_string())
    builder.write_char(',')
    builder.write_string(scenario.event_fraction.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate a planning scenario.
pub fn survival_power_scenario_is_valid(
  scenario : SurvivalPowerScenario,
) -> Bool {
  scenario.alpha > 0.0 &&
  scenario.alpha < 1.0 &&
  scenario.power > 0.0 &&
  scenario.power < 1.0 &&
  scenario.hazard_ratio > 0.0 &&
  scenario.allocation > 0.0 &&
  scenario.allocation < 1.0 &&
  scenario.expected_events >= 0 &&
  scenario.total_sample_size >= 0
}
