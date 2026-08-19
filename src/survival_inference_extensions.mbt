///|
/// A generic estimate with an uncertainty interval.
pub(all) struct SurvivalEstimateInterval {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  confidence : Double
  inference_row : Bool
} derive(Eq, Debug)

///|
/// A fixed-effect meta-analysis result.
pub(all) struct SurvivalMetaAnalysis {
  estimate : Double
  standard_error : Double
  lower : Double
  upper : Double
  heterogeneity : Double
  study_count : Int
} derive(Eq, Debug)

///|
/// Build a normal confidence interval.
pub fn survival_estimate_interval(
  estimate : Double,
  standard_error : Double,
  confidence : Double,
) -> SurvivalEstimateInterval {
  let level = if confidence <= 0.0 {
    0.95
  } else if confidence >= 1.0 {
    0.999999
  } else {
    confidence
  }
  let z = normal_quantile(0.5 + level / 2.0)
  {
    estimate,
    standard_error: standard_error.abs(),
    lower: estimate - z * standard_error.abs(),
    upper: estimate + z * standard_error.abs(),
    confidence: level,
    inference_row: true,
  }
}

///|
/// Convert an interval to a risk ratio scale.
pub fn survival_interval_exp(
  interval : SurvivalEstimateInterval,
) -> SurvivalEstimateInterval {
  {
    estimate: double_exp(interval.estimate),
    standard_error: interval.standard_error,
    lower: double_exp(interval.lower),
    upper: double_exp(interval.upper),
    confidence: interval.confidence,
    inference_row: true,
  }
}

///|
/// Compute a two-sided normal p-value.
pub fn survival_two_sided_p_value(
  estimate : Double,
  standard_error : Double,
) -> Double {
  if standard_error.abs() <= survival_epsilon {
    0.0
  } else {
    2.0 * (1.0 - normal_cdf((estimate / standard_error).abs()))
  }
}

///|
/// Compute a z-statistic.
pub fn survival_z_value(estimate : Double, standard_error : Double) -> Double {
  survival_ext_safe_ratio(estimate, standard_error)
}

///|
/// Fixed-effect inverse-variance meta-analysis.
pub fn survival_meta_analysis(
  estimates : Array[Double],
  standard_errors : Array[Double],
) -> SurvivalMetaAnalysis {
  let count = if estimates.length() < standard_errors.length() {
    estimates.length()
  } else {
    standard_errors.length()
  }
  if count == 0 {
    return {
      estimate: 0.0,
      standard_error: 0.0,
      lower: 0.0,
      upper: 0.0,
      heterogeneity: 0.0,
      study_count: 0,
    }
  }
  let mut weighted_total = 0.0
  let mut weight_total = 0.0
  for i in 0..<count {
    let weight = survival_ext_safe_ratio(
      1.0,
      standard_errors[i] * standard_errors[i],
    )
    weighted_total += weight * estimates[i]
    weight_total += weight
  }
  let estimate = survival_ext_safe_ratio(weighted_total, weight_total)
  let standard_error = double_sqrt(survival_ext_safe_ratio(1.0, weight_total))
  let mut heterogeneity = 0.0
  for i in 0..<count {
    let delta = estimates[i] - estimate
    heterogeneity += delta *
      delta /
      clamp_positive(standard_errors[i] * standard_errors[i])
  }
  let interval = survival_estimate_interval(estimate, standard_error, 0.95)
  {
    estimate,
    standard_error,
    lower: interval.lower,
    upper: interval.upper,
    heterogeneity,
    study_count: count,
  }
}

///|
/// Compute a random-effects variance-of-estimates approximation.
pub fn survival_between_study_variance(
  estimates : Array[Double],
  standard_errors : Array[Double],
) -> Double {
  let meta = survival_meta_analysis(estimates, standard_errors)
  if meta.study_count <= 1 {
    return 0.0
  }
  let mut total = 0.0
  for estimate in estimates {
    let delta = estimate - meta.estimate
    total += delta * delta
  }
  let average_variance = survival_ext_mean(standard_errors) *
    survival_ext_mean(standard_errors)
  clamp_positive(
    total / (meta.study_count - 1).to_double() - average_variance,
    floor=0.0,
  )
}

///|
/// Compute a Wald test for a difference.
pub fn survival_difference_test(
  first : Double,
  first_se : Double,
  second : Double,
  second_se : Double,
) -> SurvivalEstimateInterval {
  let difference = first - second
  let se = double_sqrt(first_se * first_se + second_se * second_se)
  survival_estimate_interval(difference, se, 0.95)
}

///|
/// Compute a multiplicity-adjusted alpha.
pub fn survival_bonferroni_alpha(alpha : Double, comparisons : Int) -> Double {
  if comparisons <= 0 {
    clamp_probability(alpha)
  } else {
    clamp_probability(alpha) / comparisons.to_double()
  }
}

///|
/// Compute a family of corrected intervals.
pub fn survival_adjusted_intervals(
  estimates : Array[Double],
  standard_errors : Array[Double],
  alpha : Double,
) -> Array[SurvivalEstimateInterval] {
  let result : Array[SurvivalEstimateInterval] = []
  let comparisons = if estimates.length() < standard_errors.length() {
    estimates.length()
  } else {
    standard_errors.length()
  }
  let confidence = 1.0 - survival_bonferroni_alpha(1.0 - alpha, comparisons)
  for i in 0..<comparisons {
    result.push(
      survival_estimate_interval(estimates[i], standard_errors[i], confidence),
    )
  }
  result
}

///|
/// Pool probabilities on the logit scale.
pub fn survival_pool_probabilities(
  probabilities : Array[Double],
  weights : Array[Double],
) -> Double {
  let logits : Array[Double] = []
  let count = if probabilities.length() < weights.length() {
    probabilities.length()
  } else {
    weights.length()
  }
  for i in 0..<count {
    logits.push(survival_ext_logit(probabilities[i]))
  }
  survival_ext_inv_logit(survival_ext_weighted_mean(logits, weights))
}

///|
/// Compute a confidence interval for a probability by the logit delta method.
pub fn survival_probability_interval(
  probability : Double,
  standard_error : Double,
  confidence : Double,
) -> SurvivalEstimateInterval {
  let logit = survival_ext_logit(probability)
  let scale = clamp_positive(probability * (1.0 - probability))
  let interval = survival_estimate_interval(
    logit,
    standard_error / scale,
    confidence,
  )
  {
    estimate: survival_ext_inv_logit(interval.estimate),
    standard_error,
    lower: survival_ext_inv_logit(interval.lower),
    upper: survival_ext_inv_logit(interval.upper),
    confidence: interval.confidence,
    inference_row: true,
  }
}

///|
/// Export a generic interval.
pub fn survival_estimate_interval_to_csv(
  interval : SurvivalEstimateInterval,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("estimate,standard_error,lower,upper,confidence\n")
  builder.write_string(interval.estimate.to_string())
  builder.write_char(',')
  builder.write_string(interval.standard_error.to_string())
  builder.write_char(',')
  builder.write_string(interval.lower.to_string())
  builder.write_char(',')
  builder.write_string(interval.upper.to_string())
  builder.write_char(',')
  builder.write_string(interval.confidence.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Validate an interval.
pub fn survival_estimate_interval_is_valid(
  interval : SurvivalEstimateInterval,
) -> Bool {
  interval.confidence > 0.0 &&
  interval.confidence < 1.0 &&
  interval.lower <= interval.estimate &&
  interval.estimate <= interval.upper &&
  interval.standard_error >= 0.0
}
