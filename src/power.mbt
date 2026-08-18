///|
/// Planning result for a two-arm time-to-event study.
pub(all) struct PowerPlan {
  hazard_ratio : Double
  alpha : Double
  target_power : Double
  allocation : Double
  required_events : Int
  estimated_sample_size : Int
} derive(Eq, Debug)

///|
/// Approximate required event count for a two-arm log-rank design.
pub fn logrank_required_events(
  hazard_ratio : Double,
  alpha : Double,
  power : Double,
  allocation? : Double = 0.5,
) -> Int {
  let hr = if hazard_ratio <= 0.0 { 1.0 } else { hazard_ratio }
  let q = if allocation <= 0.0 || allocation >= 1.0 { 0.5 } else { allocation }
  if (hr - 1.0).abs() <= survival_epsilon {
    return 2147483647
  }
  let alpha_probability = clamp_probability(1.0 - alpha / 2.0)
  let power_probability = clamp_probability(power)
  let z_alpha = normal_quantile(alpha_probability)
  let z_power = normal_quantile(power_probability)
  let information = (z_alpha + z_power) *
    (z_alpha + z_power) /
    (double_log(hr) * double_log(hr) * q * (1.0 - q))
  let events = (information + 0.999999).to_int()
  if events < 1 {
    1
  } else {
    events
  }
}

///|
/// Approximate power for an event count and hazard ratio.
pub fn logrank_power(
  events : Int,
  hazard_ratio : Double,
  alpha : Double,
  allocation? : Double = 0.5,
) -> Double {
  if events <= 0 || hazard_ratio <= 0.0 {
    0.0
  } else {
    let q = if allocation <= 0.0 || allocation >= 1.0 {
      0.5
    } else {
      allocation
    }
    let noncentrality = double_sqrt(events.to_double() * q * (1.0 - q)) *
      double_log(hazard_ratio).abs()
    let cutoff = normal_quantile(1.0 - alpha / 2.0)
    clamp_probability(1.0 - normal_cdf(cutoff - noncentrality))
  }
}

///|
/// Convert required events to a total sample size under an expected event
/// fraction and a small loss-to-follow-up allowance.
pub fn logrank_estimated_sample_size(
  required_events : Int,
  expected_event_fraction : Double,
  loss_fraction? : Double = 0.0,
) -> Int {
  let event_fraction = if expected_event_fraction <= 0.0 {
    0.01
  } else if expected_event_fraction > 1.0 {
    1.0
  } else {
    expected_event_fraction
  }
  let loss = if loss_fraction < 0.0 {
    0.0
  } else if loss_fraction >= 1.0 {
    0.95
  } else {
    loss_fraction
  }
  let size = required_events.to_double() / (event_fraction * (1.0 - loss))
  if size < 1.0 {
    1
  } else {
    (size + 0.999999).to_int()
  }
}

///|
/// Build a complete planning object.
pub fn make_power_plan(
  hazard_ratio : Double,
  alpha : Double,
  target_power : Double,
  expected_event_fraction : Double,
  allocation? : Double = 0.5,
  loss_fraction? : Double = 0.0,
) -> PowerPlan {
  let required = logrank_required_events(
    hazard_ratio,
    alpha,
    target_power,
    allocation~,
  )
  {
    hazard_ratio,
    alpha,
    target_power,
    allocation: if allocation <= 0.0 || allocation >= 1.0 {
      0.5
    } else {
      allocation
    },
    required_events: required,
    estimated_sample_size: logrank_estimated_sample_size(
      required,
      expected_event_fraction,
      loss_fraction~,
    ),
  }
}

///|
/// Find a hazard ratio that can be detected with a fixed event budget.
pub fn detectable_hazard_ratio(
  events : Int,
  alpha : Double,
  target_power : Double,
  allocation? : Double = 0.5,
) -> Double {
  if events <= 0 {
    1.0
  } else {
    let q = if allocation <= 0.0 || allocation >= 1.0 {
      0.5
    } else {
      allocation
    }
    let z = normal_quantile(1.0 - alpha / 2.0) + normal_quantile(target_power)
    double_exp(z / double_sqrt(events.to_double() * q * (1.0 - q)))
  }
}

///|
/// Render a power plan as a small Markdown table.
pub fn PowerPlan::to_markdown(self : PowerPlan) -> String {
  "| Parameter | Value |\n|---|---:|\n" +
  "| Hazard ratio | \{self.hazard_ratio} |\n" +
  "| Alpha | \{self.alpha} |\n" +
  "| Target power | \{self.target_power} |\n" +
  "| Required events | \{self.required_events} |\n" +
  "| Estimated sample size | \{self.estimated_sample_size} |\n"
}
