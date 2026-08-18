///|
/// Information criteria for survival-model comparison.
pub(all) struct ModelScore {
  name : String
  parameters : Int
  observations : Int
  log_likelihood : Double
  aic : Double
  bic : Double
  aicc : Double
} derive(Eq, Debug)

///|
/// Calculate AIC, BIC, and finite-sample corrected AIC.
pub fn model_score(
  name : String,
  parameters : Int,
  observations : Int,
  log_likelihood : Double,
) -> ModelScore {
  let aic = 2.0 * parameters.to_double() - 2.0 * log_likelihood
  let bic = parameters.to_double() * double_log(observations.to_double()) -
    2.0 * log_likelihood
  let correction = if observations - parameters - 1 <= 0 {
    0.0
  } else {
    2.0 *
    parameters.to_double() *
    (parameters + 1).to_double() /
    (observations - parameters - 1).to_double()
  }
  {
    name,
    parameters,
    observations,
    log_likelihood,
    aic,
    bic,
    aicc: aic + correction,
  }
}

///|
/// Select the lowest information-criterion model.
pub fn select_by_aic(scores : Array[ModelScore]) -> ModelScore? {
  if scores.length() == 0 {
    None
  } else {
    let mut best = scores[0]
    for score in scores {
      if score.aic < best.aic {
        best = score
      }
    }
    Some(best)
  }
}

///|
/// Select the lowest BIC model.
pub fn select_by_bic(scores : Array[ModelScore]) -> ModelScore? {
  if scores.length() == 0 {
    None
  } else {
    let mut best = scores[0]
    for score in scores {
      if score.bic < best.bic {
        best = score
      }
    }
    Some(best)
  }
}

///|
/// Likelihood-ratio test between nested models.
pub(all) struct LikelihoodRatioResult {
  statistic : Double
  degrees_of_freedom : Int
  p_value : Double
  larger_model_preferred : Bool
} derive(Eq, Debug)

///|
/// Compare nested models using twice the log-likelihood difference.
pub fn likelihood_ratio_test(
  smaller_log_likelihood : Double,
  larger_log_likelihood : Double,
  smaller_parameters : Int,
  larger_parameters : Int,
) -> LikelihoodRatioResult {
  let statistic = 2.0 * (larger_log_likelihood - smaller_log_likelihood)
  let degrees = larger_parameters - smaller_parameters
  {
    statistic,
    degrees_of_freedom: if degrees < 1 {
      1
    } else {
      degrees
    },
    p_value: chi_square_survival(
      statistic,
      if degrees < 1 {
        1
      } else {
        degrees
      },
    ),
    larger_model_preferred: statistic > 0.0,
  }
}

///|
/// Validate that an objective did not decrease during a fitting trace.
pub fn objective_trace_is_monotone(
  trace : Array[Double],
  tolerance? : Double = 0.000001,
) -> Bool {
  for i in 1..<trace.length() {
    if trace[i] + tolerance < trace[i - 1] {
      return false
    }
  }
  true
}

///|
/// Render model scores as CSV.
pub fn model_scores_to_csv(scores : Array[ModelScore]) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "name,parameters,observations,log_likelihood,aic,bic,aicc\n",
  )
  for score in scores {
    builder.write_string(
      "\{score.name},\{score.parameters},\{score.observations},\{score.log_likelihood},\{score.aic},\{score.bic},\{score.aicc}\n",
    )
  }
  builder.to_string()
}
