///|
/// A clinically oriented summary of one right-censored study arm.
pub(all) struct TrialArmSummary {
  name : String
  sample_size : Int
  events : Int
  censored : Int
  event_fraction : Double
  median_survival : Double?
  survival_at_horizon : Double
  restricted_mean : Double
  person_time : Double
} derive(Eq, Debug)

///|
/// Effect measures and decision aids for a two-arm survival comparison.
pub(all) struct TrialComparison {
  control : TrialArmSummary
  treatment : TrialArmSummary
  hazard_ratio : Double
  hazard_ratio_lower : Double
  hazard_ratio_upper : Double
  logrank_p_value : Double
  restricted_mean_difference : Double
  absolute_risk_reduction : Double
  number_needed_to_treat : Double
  interpretation : String
} derive(Eq, Debug)

///|
/// A confidence interval for a log hazard coefficient.
pub(all) struct HazardRatioInterval {
  hazard_ratio : Double
  lower : Double
  upper : Double
  standard_error : Double
} derive(Eq, Debug)

///|
/// Compute total observed follow-up time for an arm.
pub fn trial_person_time(data : SurvivalData) -> Double {
  let mut total = 0.0
  for observation in data.to_array() {
    total += observation.time
  }
  total
}

///|
/// Build a report-ready arm summary from the same data used by KM.
pub fn summarize_trial_arm(
  name : String,
  data : SurvivalData,
  horizon : Double,
) -> TrialArmSummary {
  let curve = kaplan_meier(data)
  let sample_size = data.length()
  let events = data.event_count()
  let median_survival = curve.median()
  let restricted_mean = curve.restricted_mean(horizon)
  {
    name,
    sample_size,
    events,
    censored: sample_size - events,
    event_fraction: if sample_size == 0 {
      0.0
    } else {
      events.to_double() / sample_size.to_double()
    },
    median_survival,
    survival_at_horizon: curve.survival_at(horizon),
    restricted_mean,
    person_time: trial_person_time(data),
  }
}

///|
/// Compute the absolute reduction in event risk by a fixed horizon.
pub fn absolute_risk_reduction(
  control_risk : Double,
  treatment_risk : Double,
) -> Double {
  clamp_probability(control_risk) - clamp_probability(treatment_risk)
}

///|
/// Convert an absolute risk reduction into an NNT.
pub fn number_needed_to_treat_from_arr(arr : Double) -> Double {
  if arr <= 0.000000000001 {
    0.0
  } else {
    1.0 / arr
  }
}

///|
/// Estimate events per unit of observed person-time.
pub fn trial_event_rate(data : SurvivalData) -> Double {
  let follow_up = trial_person_time(data)
  if follow_up <= survival_epsilon {
    0.0
  } else {
    data.event_count().to_double() / follow_up
  }
}

///|
/// Construct a confidence interval on the hazard-ratio scale.
pub fn hazard_ratio_interval(
  log_hazard_ratio : Double,
  standard_error : Double,
  confidence? : Double = 0.95,
) -> HazardRatioInterval {
  let level = if confidence <= 0.0 || confidence >= 1.0 {
    0.95
  } else {
    confidence
  }
  let z = normal_quantile(0.5 + level / 2.0)
  let center = double_exp(log_hazard_ratio)
  let lower = double_exp(log_hazard_ratio - z * standard_error)
  let upper = double_exp(log_hazard_ratio + z * standard_error)
  { hazard_ratio: center, lower, upper, standard_error }
}

///|
/// Extract one hazard-ratio interval from a fitted Cox model.
pub fn cox_hazard_ratio_interval(
  fit : CoxFitResult,
  coefficient_index : Int,
  confidence? : Double = 0.95,
) -> HazardRatioInterval? {
  if coefficient_index < 0 || coefficient_index >= fit.coefficients.length() {
    return None
  }
  Some(
    hazard_ratio_interval(
      fit.coefficients[coefficient_index],
      fit.standard_errors[coefficient_index],
      confidence~,
    ),
  )
}

///|
/// Compare two arms at a clinically meaningful horizon.
pub fn compare_trial_arms(
  control : TrialArmSummary,
  treatment : TrialArmSummary,
  hazard_ratio : Double,
  hazard_ratio_lower : Double,
  hazard_ratio_upper : Double,
  logrank_p_value : Double,
) -> TrialComparison {
  let arr = absolute_risk_reduction(
    1.0 - control.survival_at_horizon,
    1.0 - treatment.survival_at_horizon,
  )
  let nnt = number_needed_to_treat_from_arr(arr)
  let interpretation = clinical_interpretation(
    hazard_ratio, arr, logrank_p_value,
  )
  {
    control,
    treatment,
    hazard_ratio,
    hazard_ratio_lower,
    hazard_ratio_upper,
    logrank_p_value,
    restricted_mean_difference: treatment.restricted_mean -
    control.restricted_mean,
    absolute_risk_reduction: arr,
    number_needed_to_treat: nnt,
    interpretation,
  }
}

///|
/// Compare two raw arms without requiring callers to calculate KM first.
pub fn compare_survival_data(
  control_name : String,
  control : SurvivalData,
  treatment_name : String,
  treatment : SurvivalData,
  horizon : Double,
) -> TrialComparison {
  let control_summary = summarize_trial_arm(control_name, control, horizon)
  let treatment_summary = summarize_trial_arm(
    treatment_name, treatment, horizon,
  )
  let logrank = logrank_test(control, treatment)
  compare_trial_arms(
    control_summary,
    treatment_summary,
    1.0,
    1.0,
    1.0,
    logrank.p_value,
  )
}

///|
/// Generate a short, direction-aware interpretation for a report.
pub fn clinical_interpretation(
  hazard_ratio : Double,
  absolute_risk_reduction : Double,
  p_value : Double,
) -> String {
  let direction = if hazard_ratio < 1.0 {
    "lower estimated hazard"
  } else if hazard_ratio > 1.0 {
    "higher estimated hazard"
  } else {
    "similar estimated hazard"
  }
  let effect = if absolute_risk_reduction > 0.0 {
    "fewer events by the horizon"
  } else if absolute_risk_reduction < 0.0 {
    "more events by the horizon"
  } else {
    "no absolute risk separation by the horizon"
  }
  let evidence = if p_value < 0.05 {
    "the log-rank result is below 0.05"
  } else {
    "the log-rank result is not below 0.05"
  }
  "Treatment shows \{direction}, \{effect}; \{evidence}."
}

///|
/// Flag whether an estimated benefit is directionally favorable and passes alpha.
pub fn clinically_favorable(
  comparison : TrialComparison,
  alpha? : Double = 0.05,
) -> Bool {
  comparison.hazard_ratio < 1.0 &&
  comparison.absolute_risk_reduction > 0.0 &&
  comparison.logrank_p_value < alpha
}

///|
/// Convert one arm summary to a compact Markdown row.
pub fn TrialArmSummary::to_markdown_row(self : TrialArmSummary) -> String {
  let median = match self.median_survival {
    Some(value) => "\{value}"
    None => "not reached"
  }
  "| \{self.name} | \{self.sample_size} | \{self.events} | \{self.event_fraction} | \{median} | \{self.survival_at_horizon} | \{self.restricted_mean} |\n"
}

///|
/// Render a complete comparison with enough context for a review record.
pub fn TrialComparison::to_markdown(self : TrialComparison) -> String {
  let builder = StringBuilder::new()
  builder.write_string("## Survival comparison\n\n")
  builder.write_string(
    "| Arm | N | Events | Event fraction | Median | Survival at horizon | RMST |\n",
  )
  builder.write_string("|---|---:|---:|---:|---:|---:|---:|\n")
  builder.write_string(self.control.to_markdown_row())
  builder.write_string(self.treatment.to_markdown_row())
  builder.write_string("\n| Measure | Estimate |\n|---|---:|\n")
  builder.write_string("| Hazard ratio | \{self.hazard_ratio} |\n")
  builder.write_string(
    "| Hazard ratio lower CI | \{self.hazard_ratio_lower} |\n",
  )
  builder.write_string(
    "| Hazard ratio upper CI | \{self.hazard_ratio_upper} |\n",
  )
  builder.write_string("| Log-rank p-value | \{self.logrank_p_value} |\n")
  builder.write_string(
    "| RMST difference | \{self.restricted_mean_difference} |\n",
  )
  builder.write_string(
    "| Absolute risk reduction | \{self.absolute_risk_reduction} |\n",
  )
  builder.write_string(
    "| Number needed to treat | \{self.number_needed_to_treat} |\n\n",
  )
  builder.write_string(self.interpretation)
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Calculate a survival-based risk ratio at a fixed horizon.
pub fn survival_risk_ratio(
  control_survival : Double,
  treatment_survival : Double,
) -> Double {
  let control_risk = 1.0 - clamp_probability(control_survival)
  let treatment_risk = 1.0 - clamp_probability(treatment_survival)
  if control_risk <= survival_epsilon {
    0.0
  } else {
    treatment_risk / control_risk
  }
}

///|
/// Convert an NNT to an absolute risk reduction safely.
pub fn absolute_risk_reduction_from_nnt(nnt : Double) -> Double {
  if nnt <= 0.0 {
    0.0
  } else {
    1.0 / nnt
  }
}

///|
/// Return a stable textual label for an effect estimate.
pub fn hazard_ratio_label(hazard_ratio : Double) -> String {
  if hazard_ratio < 0.8 {
    "strongly favorable"
  } else if hazard_ratio < 1.0 {
    "favorable"
  } else if hazard_ratio <= 1.25 {
    "near null"
  } else {
    "unfavorable"
  }
}
