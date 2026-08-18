///|
/// Schoenfeld residual for one event time and one covariate.
pub(all) struct SchoenfeldResidual {
  time : Double
  covariate : Int
  observed : Double
  expected : Double
  residual : Double
} derive(Eq, Debug)

///|
/// Compute tie-aware Schoenfeld residuals under a fitted Cox coefficient.
pub fn cox_schoenfeld_residuals(
  data : CoxData,
  beta : Array[Double],
) -> Array[SchoenfeldResidual] {
  let result = []
  let mut index = 0
  while index < data.length() {
    let time = data.get(index).time
    let event_rows = []
    while index < data.length() && data.get(index).time == time {
      if data.get(index).event == 1 {
        event_rows.push(data.get(index).covariates)
      }
      index += 1
    }
    if event_rows.length() > 0 {
      let risk_sum = cox_risk_sum(data, beta, time)
      let expected = zeros(data.n_covariates)
      for observation in data.observations {
        if observation.time >= time {
          let weight = double_exp(dot_product(observation.covariates, beta))
          for j in 0..<data.n_covariates {
            expected[j] += weight * observation.covariates[j]
          }
        }
      }
      for event in event_rows {
        for j in 0..<data.n_covariates {
          let expected_value = expected[j] / risk_sum
          result.push({
            time,
            covariate: j,
            observed: event[j],
            expected: expected_value,
            residual: event[j] - expected_value,
          })
        }
      }
    }
  }
  result
}

///|
/// Martingale residual for one subject-level Cox observation.
pub(all) struct MartingaleResidual {
  time : Double
  event : Int
  risk : Double
  expected_events : Double
  residual : Double
} derive(Eq, Debug)

///|
/// Compute martingale residuals from a fitted baseline hazard.
pub fn cox_martingale_residuals(
  data : CoxData,
  fit : CoxFitResult,
) -> Array[MartingaleResidual] {
  let result = []
  for observation in data.observations {
    let risk = double_exp(dot_product(observation.covariates, fit.coefficients))
    let expected = cox_baseline_cumulative_hazard(
        fit.baseline,
        observation.time,
      ) *
      risk
    result.push({
      time: observation.time,
      event: observation.event,
      risk,
      expected_events: expected,
      residual: observation.event.to_double() - expected,
    })
  }
  result
}

///|
/// Proportional-hazards diagnostic from residual/time correlations.
pub(all) struct PHDiagnostic {
  covariate : Int
  residual_time_correlation : Double
  residual_mean : Double
  residual_variance : Double
  warning : Bool
} derive(Eq, Debug)

///|
/// Calculate a simple residual trend diagnostic. A large absolute correlation
/// warns that the proportional-hazards assumption may need stratification or
/// time interactions.
pub fn cox_ph_diagnostics(
  data : CoxData,
  fit : CoxFitResult,
  warning_threshold? : Double = 0.25,
) -> Array[PHDiagnostic] {
  let residuals = cox_schoenfeld_residuals(data, fit.coefficients)
  let result = []
  for covariate in 0..<data.n_covariates {
    let times = []
    let values = []
    for residual in residuals {
      if residual.covariate == covariate {
        times.push(residual.time)
        values.push(residual.residual)
      }
    }
    let association = if times.length() < 2 {
      0.0
    } else {
      correlation(times, values)
    }
    let summary = residual_summary(values)
    result.push({
      covariate,
      residual_time_correlation: association,
      residual_mean: summary.mean,
      residual_variance: summary.variance,
      warning: association.abs() > warning_threshold,
    })
  }
  result
}

///|
/// A calibration bin for predicted event risk.
pub(all) struct CalibrationBin {
  lower_risk : Double
  upper_risk : Double
  count : Int
  mean_predicted_risk : Double
  observed_event_fraction : Double
  absolute_error : Double
} derive(Eq, Debug)

///|
/// Compare predicted event probabilities with observed outcomes by risk bin.
pub fn calibration_curve(
  data : SurvivalData,
  predicted_survival : Array[Double],
  horizon : Double,
  bins : Int,
) -> Array[CalibrationBin] {
  if data.length() != predicted_survival.length() || bins <= 0 {
    panic()
  }
  let result = []
  for bin in 0..<bins {
    let lower = bin.to_double() / bins.to_double()
    let upper = (bin + 1).to_double() / bins.to_double()
    let mut count = 0
    let mut predicted = 0.0
    let mut observed = 0
    for i in 0..<data.length() {
      let risk = 1.0 - clamp_probability(predicted_survival[i])
      let in_bin = if bin == bins - 1 {
        risk >= lower && risk <= upper
      } else {
        risk >= lower && risk < upper
      }
      if in_bin {
        count += 1
        predicted += risk
        if data.get(i).event == 1 && data.get(i).time <= horizon {
          observed += 1
        }
      }
    }
    let mean_predicted = if count == 0 {
      0.0
    } else {
      predicted / count.to_double()
    }
    let observed_fraction = if count == 0 {
      0.0
    } else {
      observed.to_double() / count.to_double()
    }
    result.push({
      lower_risk: lower,
      upper_risk: upper,
      count,
      mean_predicted_risk: mean_predicted,
      observed_event_fraction: observed_fraction,
      absolute_error: (mean_predicted - observed_fraction).abs(),
    })
  }
  result
}

///|
/// Concordance-index summary for a risk score.
pub(all) struct ConcordanceResult {
  concordant : Int
  discordant : Int
  tied : Int
  comparable : Int
  c_index : Double
} derive(Eq, Debug)

///|
/// Harrell's concordance index for pairs where the earlier observation is an
/// observed event. Larger scores must indicate higher risk.
pub fn concordance_index(
  data : SurvivalData,
  risk_scores : Array[Double],
) -> ConcordanceResult {
  if data.length() != risk_scores.length() {
    panic()
  }
  let mut concordant = 0
  let mut discordant = 0
  let mut tied = 0
  for i in 0..<data.length() {
    if data.get(i).event != 1 {
      continue
    }
    for j in 0..<data.length() {
      if data.get(j).time <= data.get(i).time {
        continue
      }
      if risk_scores[i] > risk_scores[j] {
        concordant += 1
      } else if risk_scores[i] < risk_scores[j] {
        discordant += 1
      } else {
        tied += 1
      }
    }
  }
  let comparable = concordant + discordant + tied
  let c_index = if comparable == 0 {
    0.5
  } else {
    (concordant.to_double() + tied.to_double() / 2.0) / comparable.to_double()
  }
  { concordant, discordant, tied, comparable, c_index }
}

///|
/// A simple time-dependent AUC treating known outcomes by a horizon.
pub fn time_dependent_auc(
  data : SurvivalData,
  risk_scores : Array[Double],
  horizon : Double,
) -> Double {
  if data.length() != risk_scores.length() {
    panic()
  }
  let cases = []
  let controls = []
  for i in 0..<data.length() {
    let observation = data.get(i)
    if observation.event == 1 && observation.time <= horizon {
      cases.push(risk_scores[i])
    } else if observation.time > horizon {
      controls.push(risk_scores[i])
    }
  }
  if cases.length() == 0 || controls.length() == 0 {
    return 0.5
  }
  let mut wins = 0.0
  for case_score in cases {
    for control_score in controls {
      if case_score > control_score {
        wins += 1.0
      } else if case_score == control_score {
        wins += 0.5
      }
    }
  }
  wins / (cases.length() * controls.length()).to_double()
}
