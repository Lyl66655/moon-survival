///|
/// A sensitivity scenario for unmeasured censoring or outcome shift.
pub(all) struct SurvivalSensitivityScenario {
  name : String
  hazard_multiplier : Double
  time_shift : Double
  censoring_fraction : Double
  adjusted_risk : Double
  adjusted_survival : Double
  scenario_row : Bool
} derive(Eq, Debug)

///|
/// A tipping-point result over a two-dimensional grid.
pub(all) struct SurvivalTippingPoint {
  treatment_multiplier : Double
  control_multiplier : Double
  treatment_risk : Double
  control_risk : Double
  risk_difference : Double
  conclusion : String
} derive(Eq, Debug)

///|
/// An E-value style robustness summary for a risk ratio.
pub(all) struct SurvivalEValue {
  risk_ratio : Double
  e_value : Double
  lower_bound : Double
  interpretation : String
} derive(Eq, Debug)

///|
/// An interval of plausible adjusted estimates.
pub(all) struct SurvivalSensitivityInterval {
  lower : Double
  upper : Double
  center : Double
  width : Double
  interval_row : Bool
} derive(Eq, Debug)

///|
/// Adjust a cumulative incidence for a multiplicative latent-risk factor.
pub fn sensitivity_adjusted_risk(
  observed_risk : Double,
  hazard_multiplier : Double,
  censoring_fraction : Double,
) -> Double {
  let risk = clamp_probability(observed_risk)
  let multiplier = if hazard_multiplier <= 0.0 {
    1.0
  } else {
    hazard_multiplier
  }
  let censoring = clamp_probability(censoring_fraction)
  clamp_probability(risk * multiplier * (1.0 - censoring) + censoring * risk)
}

///|
/// Adjust a survival curve by a time shift.
pub fn sensitivity_adjusted_survival(
  survival : Double,
  time_shift : Double,
  scale : Double,
) -> Double {
  let safe_scale = clamp_positive(scale)
  clamp_probability(
    double_exp(double_log(clamp_positive(survival)) - time_shift / safe_scale),
  )
}

///|
/// Evaluate a single sensitivity scenario.
pub fn sensitivity_scenario(
  name : String,
  observed_risk : Double,
  observed_survival : Double,
  hazard_multiplier : Double,
  time_shift : Double,
  censoring_fraction : Double,
  scale : Double,
) -> SurvivalSensitivityScenario {
  {
    name,
    hazard_multiplier,
    time_shift,
    censoring_fraction,
    adjusted_risk: sensitivity_adjusted_risk(
      observed_risk, hazard_multiplier, censoring_fraction,
    ),
    adjusted_survival: sensitivity_adjusted_survival(
      observed_survival, time_shift, scale,
    ),
    scenario_row: true,
  }
}

///|
/// Generate a rectangular sensitivity grid.
pub fn sensitivity_grid(
  observed_risk_treatment : Double,
  observed_risk_control : Double,
  treatment_multipliers : Array[Double],
  control_multipliers : Array[Double],
) -> Array[SurvivalTippingPoint] {
  let result : Array[SurvivalTippingPoint] = []
  for treatment_multiplier in treatment_multipliers {
    for control_multiplier in control_multipliers {
      let treatment_risk = sensitivity_adjusted_risk(
        observed_risk_treatment, treatment_multiplier, 0.0,
      )
      let control_risk = sensitivity_adjusted_risk(
        observed_risk_control, control_multiplier, 0.0,
      )
      let difference = treatment_risk - control_risk
      result.push({
        treatment_multiplier,
        control_multiplier,
        treatment_risk,
        control_risk,
        risk_difference: difference,
        conclusion: if difference < 0.0 {
          "treatment_lower_risk"
        } else if difference > 0.0 {
          "control_lower_risk"
        } else {
          "tie"
        },
      })
    }
  }
  result
}

///|
/// Compute an E-value for a risk ratio away from the null.
pub fn sensitivity_e_value(
  risk_ratio : Double,
  lower_bound : Double,
) -> SurvivalEValue {
  let ratio = if risk_ratio >= 1.0 {
    risk_ratio
  } else {
    1.0 / clamp_positive(risk_ratio)
  }
  let e_value = ratio + double_sqrt(ratio * (ratio - 1.0))
  {
    risk_ratio,
    e_value,
    lower_bound,
    interpretation: if lower_bound > 1.0 {
      "association remains above null"
    } else {
      "uncertainty includes null"
    },
  }
}

///|
/// Compute a tipping point for a risk-difference decision.
pub fn sensitivity_tipping_point(
  treatment_risk : Double,
  control_risk : Double,
  treatment_shift : Double,
  control_shift : Double,
) -> SurvivalTippingPoint {
  let adjusted_treatment = clamp_probability(treatment_risk + treatment_shift)
  let adjusted_control = clamp_probability(control_risk + control_shift)
  let difference = adjusted_treatment - adjusted_control
  {
    treatment_multiplier: 1.0,
    control_multiplier: 1.0,
    treatment_risk: adjusted_treatment,
    control_risk: adjusted_control,
    risk_difference: difference,
    conclusion: if difference < 0.0 {
      "treatment_lower_risk"
    } else {
      "not_lower_risk"
    },
  }
}

///|
/// Compute lower and upper bounds under an additive departure from MAR.
pub fn sensitivity_bounds(
  estimate : Double,
  standard_error : Double,
  departure : Double,
) -> SurvivalSensitivityInterval {
  let margin = 1.96 * standard_error.abs() + departure.abs()
  let lower = estimate - margin
  let upper = estimate + margin
  { lower, upper, center: estimate, width: upper - lower, interval_row: true }
}

///|
/// Apply a monotone delta adjustment to subject-level risks.
pub fn sensitivity_delta_adjust(
  risks : Array[Double],
  delta : Double,
) -> Array[Double] {
  let result : Array[Double] = []
  for risk in risks {
    result.push(clamp_probability(risk + delta))
  }
  result
}

///|
/// Summarize an adjusted risk vector.
pub fn sensitivity_risk_summary(
  risks : Array[Double],
) -> SurvivalSensitivityInterval {
  if risks.length() == 0 {
    return {
      lower: 0.0,
      upper: 0.0,
      center: 0.0,
      width: 0.0,
      interval_row: true,
    }
  }
  let sorted : Array[Double] = []
  for risk in risks {
    sorted.push(risk)
  }
  sorted.sort()
  {
    lower: sorted[0],
    upper: sorted[sorted.length() - 1],
    center: mean_doubles(sorted),
    width: sorted[sorted.length() - 1] - sorted[0],
    interval_row: true,
  }
}

///|
/// Export sensitivity scenarios.
pub fn sensitivity_scenarios_to_csv(
  scenarios : Array[SurvivalSensitivityScenario],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "name,hazard_multiplier,time_shift,censoring_fraction,adjusted_risk,adjusted_survival\n",
  )
  for scenario in scenarios {
    builder.write_string(scenario.name)
    builder.write_char(',')
    builder.write_string(scenario.hazard_multiplier.to_string())
    builder.write_char(',')
    builder.write_string(scenario.time_shift.to_string())
    builder.write_char(',')
    builder.write_string(scenario.censoring_fraction.to_string())
    builder.write_char(',')
    builder.write_string(scenario.adjusted_risk.to_string())
    builder.write_char(',')
    builder.write_string(scenario.adjusted_survival.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Check an E-value and scenario grid.
pub fn sensitivity_e_value_is_valid(value : SurvivalEValue) -> Bool {
  value.e_value >= 1.0 && !value.e_value.is_nan() && !value.e_value.is_inf()
}
