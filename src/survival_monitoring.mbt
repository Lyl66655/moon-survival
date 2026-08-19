///|
/// Group-sequential monitoring design.
pub(all) struct MonitoringDesign {
  looks : Int
  alpha : Double
  spending : String
  information_fractions : Array[Double]
  design_row : Bool
} derive(Eq, Debug)

///|
/// One interim analysis result.
pub(all) struct InterimLook {
  look : Int
  information_fraction : Double
  events : Int
  z_value : Double
  efficacy_boundary : Double
  futility_boundary : Double
  p_value : Double
  decision : String
} derive(Eq, Debug)

///|
/// Event-accrual forecast.
pub(all) struct EventAccrualForecast {
  current_events : Int
  target_events : Int
  current_time : Double
  projected_time : Double
  event_rate : Double
  completion_fraction : Double
  forecast_row : Bool
} derive(Eq, Debug)

///|
/// Create equally-spaced information looks.
pub fn monitoring_design(
  looks : Int,
  alpha : Double,
  spending : String,
) -> MonitoringDesign {
  let count = if looks <= 0 { 1 } else { looks }
  let bounded_alpha = if alpha <= 0.0 {
    0.000001
  } else if alpha >= 1.0 {
    0.999999
  } else {
    alpha
  }
  let fractions : Array[Double] = []
  for look in 1..<(count + 1) {
    fractions.push(look.to_double() / count.to_double())
  }
  {
    looks: count,
    alpha: bounded_alpha,
    spending,
    information_fractions: fractions,
    design_row: true,
  }
}

///|
/// Lan-DeMets-like alpha spending function.
pub fn monitoring_alpha_spent(
  design : MonitoringDesign,
  information_fraction : Double,
) -> Double {
  let fraction = clamp_probability(information_fraction)
  let spent = match design.spending {
    "pocock" => design.alpha * fraction
    "obrien-fleming" =>
      design.alpha *
      (
        2.0 -
        2.0 *
        normal_cdf(
          normal_quantile(1.0 - design.alpha / 2.0) /
          double_sqrt(clamp_positive(fraction)),
        )
      )
    _ => design.alpha * fraction * fraction
  }
  if spent > design.alpha {
    design.alpha
  } else {
    spent
  }
}

///|
/// Convert an alpha spent at a look to a two-sided z boundary.
pub fn monitoring_boundary(alpha_spent : Double) -> Double {
  normal_quantile(1.0 - clamp_probability(alpha_spent) / 2.0)
}

///|
/// Evaluate an interim z statistic under a monitoring design.
pub fn monitoring_interim_look(
  design : MonitoringDesign,
  look : Int,
  events : Int,
  z_value : Double,
  futility_boundary : Double,
) -> InterimLook {
  let index = if look <= 0 {
    0
  } else if look > design.looks {
    design.looks - 1
  } else {
    look - 1
  }
  let fraction = design.information_fractions[index]
  let spent = monitoring_alpha_spent(design, fraction)
  let efficacy = monitoring_boundary(spent)
  let p_value = 2.0 * (1.0 - normal_cdf(z_value.abs()))
  let decision = if z_value >= efficacy {
    "stop_efficacy"
  } else if z_value <= futility_boundary {
    "stop_futility"
  } else {
    "continue"
  }
  {
    look: index + 1,
    information_fraction: fraction,
    events,
    z_value,
    efficacy_boundary: efficacy,
    futility_boundary,
    p_value,
    decision,
  }
}

///|
/// Forecast calendar time needed to reach an event target.
pub fn monitoring_forecast(
  times : Array[Double],
  events : Array[Int],
  target_events : Int,
) -> EventAccrualForecast {
  if times.length() != events.length() ||
    times.length() == 0 ||
    target_events <= 0 {
    return {
      current_events: 0,
      target_events,
      current_time: 0.0,
      projected_time: 0.0,
      event_rate: 0.0,
      completion_fraction: 0.0,
      forecast_row: false,
    }
  }
  let mut count = 0
  let mut last_time = 0.0
  for i in 0..<times.length() {
    count += events[i]
    if times[i] > last_time {
      last_time = times[i]
    }
  }
  let rate = if last_time <= survival_epsilon {
    0.0
  } else {
    count.to_double() / last_time
  }
  let projected = if rate <= survival_epsilon {
    0.0
  } else {
    target_events.to_double() / rate
  }
  {
    current_events: count,
    target_events,
    current_time: last_time,
    projected_time: projected,
    event_rate: rate,
    completion_fraction: clamp_probability(
      count.to_double() / target_events.to_double(),
    ),
    forecast_row: true,
  }
}

///|
/// Compute information fraction from observed and target events.
pub fn monitoring_information_fraction(
  events : Int,
  target_events : Int,
) -> Double {
  if target_events <= 0 {
    0.0
  } else {
    clamp_probability(events.to_double() / target_events.to_double())
  }
}

///|
/// Build interim looks from an event sequence and z-statistic sequence.
pub fn monitoring_sequence(
  design : MonitoringDesign,
  event_counts : Array[Int],
  z_values : Array[Double],
  futility_boundary : Double,
) -> Array[InterimLook] {
  let result : Array[InterimLook] = []
  let count = if event_counts.length() < z_values.length() {
    event_counts.length()
  } else {
    z_values.length()
  }
  for i in 0..<count {
    result.push(
      monitoring_interim_look(
        design,
        i + 1,
        event_counts[i],
        z_values[i],
        futility_boundary,
      ),
    )
  }
  result
}

///|
/// Return whether any efficacy boundary has been crossed.
pub fn monitoring_efficacy_crossed(looks : Array[InterimLook]) -> Bool {
  for look in looks {
    if look.decision == "stop_efficacy" {
      return true
    }
  }
  false
}

///|
/// Return the first stopping look, when any.
pub fn monitoring_first_stop(looks : Array[InterimLook]) -> InterimLook? {
  for look in looks {
    if look.decision != "continue" {
      return Some(look)
    }
  }
  None
}

///|
/// Compute a conditional power approximation at an interim look.
pub fn monitoring_conditional_power(
  z_value : Double,
  information_fraction : Double,
  planned_effect : Double,
  alpha : Double,
) -> Double {
  let fraction = clamp_positive(information_fraction)
  let boundary = monitoring_boundary(alpha)
  let future_mean = planned_effect * double_sqrt(1.0 - fraction)
  clamp_probability(
    1.0 -
    normal_cdf(
      (boundary - z_value * double_sqrt(fraction) - future_mean) /
      double_sqrt(1.0 - fraction),
    ),
  )
}

///|
/// Compute a simple sample-size event target for a log-rank comparison.
pub fn monitoring_event_target(
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
  let information = (z_alpha + z_power) *
    (z_alpha + z_power) /
    (effect * effect)
  (information / (allocation * (1.0 - allocation))).ceil().to_int()
}

///|
/// Serialize interim decisions.
pub fn monitoring_looks_to_csv(looks : Array[InterimLook]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "look,information_fraction,events,z_value,efficacy_boundary,futility_boundary,p_value,decision\n",
  )
  for look in looks {
    builder.write_string(look.look.to_string())
    builder.write_char(',')
    builder.write_string(look.information_fraction.to_string())
    builder.write_char(',')
    builder.write_string(look.events.to_string())
    builder.write_char(',')
    builder.write_string(look.z_value.to_string())
    builder.write_char(',')
    builder.write_string(look.efficacy_boundary.to_string())
    builder.write_char(',')
    builder.write_string(look.futility_boundary.to_string())
    builder.write_char(',')
    builder.write_string(look.p_value.to_string())
    builder.write_char(',')
    builder.write_string(look.decision)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Check a monitoring design.
pub fn monitoring_design_is_valid(design : MonitoringDesign) -> Bool {
  if design.looks <= 0 ||
    design.alpha <= 0.0 ||
    design.alpha >= 1.0 ||
    design.information_fractions.length() != design.looks {
    return false
  }
  let mut previous = 0.0
  for fraction in design.information_fractions {
    if fraction <= previous || fraction > 1.0 {
      return false
    }
    previous = fraction
  }
  true
}
