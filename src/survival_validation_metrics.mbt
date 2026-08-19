///|
/// A validation metric with a direction and display label.
pub(all) struct SurvivalValidationMetric {
  name : String
  value : Double
  standard_error : Double
  higher_is_better : Bool
  metric_row : Bool
} derive(Eq, Debug)

///|
/// Compute accuracy for a horizon outcome.
pub fn survival_accuracy(
  data : SurvivalData,
  risks : Array[Double],
  threshold : Double,
  horizon : Double,
) -> Double {
  let outcomes = survival_ext_horizon_outcome(data, horizon)
  if outcomes.length() != risks.length() || outcomes.length() == 0 {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<outcomes.length() {
    if (risks[i] >= threshold && outcomes[i] == 1) ||
      (risks[i] < threshold && outcomes[i] == 0) {
      correct += 1
    }
  }
  correct.to_double() / outcomes.length().to_double()
}

///|
/// Compute sensitivity and specificity together.
pub fn survival_sensitivity_specificity(
  data : SurvivalData,
  risks : Array[Double],
  threshold : Double,
  horizon : Double,
) -> (Double, Double) {
  let metrics = survival_classification_metrics(risks, data, threshold, horizon)
  (metrics.sensitivity, metrics.specificity)
}

///|
/// Compute a simple F1 score.
pub fn survival_f1_score(
  data : SurvivalData,
  risks : Array[Double],
  threshold : Double,
  horizon : Double,
) -> Double {
  let metrics = survival_classification_metrics(risks, data, threshold, horizon)
  survival_ext_safe_ratio(
    2.0 * metrics.positive_predictive_value * metrics.sensitivity,
    metrics.positive_predictive_value + metrics.sensitivity,
  )
}

///|
/// Compute a Matthews correlation coefficient.
pub fn survival_matthews_correlation(
  data : SurvivalData,
  risks : Array[Double],
  threshold : Double,
  horizon : Double,
) -> Double {
  let outcomes = survival_ext_horizon_outcome(data, horizon)
  let mut tp = 0.0
  let mut tn = 0.0
  let mut fp = 0.0
  let mut fn_value = 0.0
  for i in 0..<outcomes.length() {
    let predicted = risks[i] >= threshold
    if predicted && outcomes[i] == 1 {
      tp += 1.0
    } else if !predicted && outcomes[i] == 0 {
      tn += 1.0
    } else if predicted {
      fp += 1.0
    } else {
      fn_value += 1.0
    }
  }
  survival_ext_safe_ratio(
    tp * tn - fp * fn_value,
    double_sqrt((tp + fp) * (tp + fn_value) * (tn + fp) * (tn + fn_value)),
  )
}

///|
/// Compute a precision-recall area on a threshold grid.
pub fn survival_precision_recall_area(
  data : SurvivalData,
  risks : Array[Double],
  thresholds : Array[Double],
  horizon : Double,
) -> Double {
  if thresholds.length() < 2 {
    return 0.0
  }
  let mut area = 0.0
  for i in 1..<thresholds.length() {
    let first = survival_classification_metrics(
      risks,
      data,
      thresholds[i - 1],
      horizon,
    )
    let second = survival_classification_metrics(
      risks,
      data,
      thresholds[i],
      horizon,
    )
    area += (thresholds[i] - thresholds[i - 1]).abs() *
      (first.positive_predictive_value + second.positive_predictive_value) /
      2.0
  }
  area
}

///|
/// Compute a collection of standard validation metrics.
pub fn survival_validation_metrics(
  data : SurvivalData,
  risks : Array[Double],
  horizon : Double,
) -> Array[SurvivalValidationMetric] {
  let result : Array[SurvivalValidationMetric] = []
  let c_index = concordance_index(data, risks)
  result.push({
    name: "c_index",
    value: c_index.c_index,
    standard_error: 0.0,
    higher_is_better: true,
    metric_row: true,
  })
  result.push({
    name: "auc",
    value: time_dependent_auc(data, risks, horizon),
    standard_error: 0.0,
    higher_is_better: true,
    metric_row: true,
  })
  result.push({
    name: "brier",
    value: ipcw_brier_score(data, survival_ext_scale(risks, -1.0), horizon),
    standard_error: 0.0,
    higher_is_better: false,
    metric_row: true,
  })
  result.push({
    name: "accuracy",
    value: survival_accuracy(data, risks, 0.5, horizon),
    standard_error: 0.0,
    higher_is_better: true,
    metric_row: true,
  })
  result
}

///|
/// Compute an optimism-corrected metric.
pub fn survival_metric_corrected(
  apparent : Double,
  optimism : Double,
  standard_error : Double,
) -> SurvivalValidationMetric {
  {
    name: "corrected",
    value: apparent - optimism,
    standard_error,
    higher_is_better: true,
    metric_row: true,
  }
}

///|
/// Select the best threshold by net benefit.
pub fn survival_best_threshold(
  metrics : Array[SurvivalClassificationMetrics],
) -> Double {
  if metrics.length() == 0 {
    return 0.5
  }
  let mut best = metrics[0]
  for metric in metrics {
    if metric.net_benefit > best.net_benefit {
      best = metric
    }
  }
  best.threshold
}

///|
/// Export validation metrics.
pub fn survival_validation_metrics_to_csv(
  metrics : Array[SurvivalValidationMetric],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("name,value,standard_error,higher_is_better\n")
  for metric in metrics {
    builder.write_string(metric.name)
    builder.write_char(',')
    builder.write_string(metric.value.to_string())
    builder.write_char(',')
    builder.write_string(metric.standard_error.to_string())
    builder.write_char(',')
    builder.write_string(metric.higher_is_better.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Verify validation metric values.
pub fn survival_validation_metric_is_valid(
  metric : SurvivalValidationMetric,
) -> Bool {
  !metric.value.is_nan() &&
  !metric.value.is_inf() &&
  metric.standard_error >= 0.0
}
