///|
/// A patient-level risk score with a fixed prediction horizon.
pub(all) struct SurvivalRiskScore {
  row_id : Int
  horizon : Double
  risk : Double
  survival : Double
  linear_predictor : Double
  score_row : Bool
} derive(Eq, Debug)

///|
/// A risk band summary for operational triage.
pub(all) struct SurvivalRiskBand {
  lower : Double
  upper : Double
  count : Int
  events : Int
  event_rate : Double
  mean_risk : Double
  band_row : Bool
} derive(Eq, Debug)

///|
/// Classification and reclassification metrics.
pub(all) struct SurvivalClassificationMetrics {
  threshold : Double
  sensitivity : Double
  specificity : Double
  positive_predictive_value : Double
  negative_predictive_value : Double
  net_benefit : Double
  sample_size : Int
} derive(Eq, Debug)

///|
/// Net reclassification improvement between two risk models.
pub(all) struct SurvivalReclassification {
  event_up : Int
  event_down : Int
  nonevent_up : Int
  nonevent_down : Int
  nri_events : Double
  nri_nonevents : Double
  nri : Double
} derive(Eq, Debug)

///|
/// Convert probabilities to stable risk-score rows.
pub fn survival_risk_scores(
  probabilities : Array[Double],
  horizon : Double,
) -> Array[SurvivalRiskScore] {
  let result : Array[SurvivalRiskScore] = []
  for i in 0..<probabilities.length() {
    let survival = clamp_probability(probabilities[i])
    result.push({
      row_id: i,
      horizon,
      risk: 1.0 - survival,
      survival,
      linear_predictor: -double_log(clamp_positive(survival)),
      score_row: true,
    })
  }
  result
}

///|
/// Return a numeric risk vector.
pub fn survival_risk_values(scores : Array[SurvivalRiskScore]) -> Array[Double] {
  let result : Array[Double] = []
  for score in scores {
    result.push(score.risk)
  }
  result
}

///|
/// Build equal-width risk bands.
pub fn survival_risk_bands(
  risks : Array[Double],
  data : SurvivalData,
  bins : Int,
  horizon : Double,
) -> Array[SurvivalRiskBand] {
  if risks.length() != data.length() || bins <= 0 {
    return []
  }
  let result : Array[SurvivalRiskBand] = []
  for bin in 0..<bins {
    let lower = bin.to_double() / bins.to_double()
    let upper = (bin + 1).to_double() / bins.to_double()
    let mut count = 0
    let mut events = 0
    let mut risk_total = 0.0
    for i in 0..<risks.length() {
      let risk = clamp_probability(risks[i])
      let in_band = if bin == bins - 1 {
        risk >= lower && risk <= upper
      } else {
        risk >= lower && risk < upper
      }
      if in_band {
        count += 1
        risk_total += risk
        if data.get(i).event == 1 && data.get(i).time <= horizon {
          events += 1
        }
      }
    }
    result.push({
      lower,
      upper,
      count,
      events,
      event_rate: if count == 0 {
        0.0
      } else {
        events.to_double() / count.to_double()
      },
      mean_risk: if count == 0 {
        0.0
      } else {
        risk_total / count.to_double()
      },
      band_row: true,
    })
  }
  result
}

///|
/// Compute threshold metrics for a binary horizon endpoint.
pub fn survival_classification_metrics(
  risks : Array[Double],
  data : SurvivalData,
  threshold : Double,
  horizon : Double,
) -> SurvivalClassificationMetrics {
  if risks.length() != data.length() {
    return {
      threshold,
      sensitivity: 0.0,
      specificity: 0.0,
      positive_predictive_value: 0.0,
      negative_predictive_value: 0.0,
      net_benefit: 0.0,
      sample_size: 0,
    }
  }
  let mut true_positive = 0
  let mut false_positive = 0
  let mut true_negative = 0
  let mut false_negative = 0
  for i in 0..<data.length() {
    let outcome = data.get(i).event == 1 && data.get(i).time <= horizon
    let predicted = risks[i] >= threshold
    if predicted && outcome {
      true_positive += 1
    } else if predicted {
      false_positive += 1
    } else if outcome {
      false_negative += 1
    } else {
      true_negative += 1
    }
  }
  let cases = true_positive + false_negative
  let controls = true_negative + false_positive
  let sensitivity = if cases == 0 {
    0.0
  } else {
    true_positive.to_double() / cases.to_double()
  }
  let specificity = if controls == 0 {
    0.0
  } else {
    true_negative.to_double() / controls.to_double()
  }
  let predicted_positive = true_positive + false_positive
  let predicted_negative = true_negative + false_negative
  let ppv = if predicted_positive == 0 {
    0.0
  } else {
    true_positive.to_double() / predicted_positive.to_double()
  }
  let npv = if predicted_negative == 0 {
    0.0
  } else {
    true_negative.to_double() / predicted_negative.to_double()
  }
  let net_benefit = if data.length() == 0 {
    0.0
  } else {
    true_positive.to_double() / data.length().to_double() -
    false_positive.to_double() /
    data.length().to_double() *
    threshold /
    clamp_positive(1.0 - threshold)
  }
  {
    threshold,
    sensitivity,
    specificity,
    positive_predictive_value: ppv,
    negative_predictive_value: npv,
    net_benefit,
    sample_size: data.length(),
  }
}

///|
/// Evaluate a threshold grid for decision-curve analysis.
pub fn survival_decision_curve(
  risks : Array[Double],
  data : SurvivalData,
  thresholds : Array[Double],
  horizon : Double,
) -> Array[SurvivalClassificationMetrics] {
  let result : Array[SurvivalClassificationMetrics] = []
  for threshold in thresholds {
    result.push(
      survival_classification_metrics(risks, data, threshold, horizon),
    )
  }
  result
}

///|
/// Calculate net reclassification improvement for a binary endpoint.
pub fn survival_reclassification(
  old_risks : Array[Double],
  new_risks : Array[Double],
  data : SurvivalData,
  horizon : Double,
  threshold : Double,
) -> SurvivalReclassification {
  if old_risks.length() != data.length() || new_risks.length() != data.length() {
    return {
      event_up: 0,
      event_down: 0,
      nonevent_up: 0,
      nonevent_down: 0,
      nri_events: 0.0,
      nri_nonevents: 0.0,
      nri: 0.0,
    }
  }
  let mut event_up = 0
  let mut event_down = 0
  let mut nonevent_up = 0
  let mut nonevent_down = 0
  for i in 0..<data.length() {
    let old_high = old_risks[i] >= threshold
    let new_high = new_risks[i] >= threshold
    let event = data.get(i).event == 1 && data.get(i).time <= horizon
    if event && new_high && !old_high {
      event_up += 1
    }
    if event && !new_high && old_high {
      event_down += 1
    }
    if !event && new_high && !old_high {
      nonevent_up += 1
    }
    if !event && !new_high && old_high {
      nonevent_down += 1
    }
  }
  let event_total = event_up + event_down
  let nonevent_total = nonevent_up + nonevent_down
  let event_component = if event_total == 0 {
    0.0
  } else {
    (event_up - event_down).to_double() / event_total.to_double()
  }
  let nonevent_component = if nonevent_total == 0 {
    0.0
  } else {
    (nonevent_down - nonevent_up).to_double() / nonevent_total.to_double()
  }
  {
    event_up,
    event_down,
    nonevent_up,
    nonevent_down,
    nri_events: event_component,
    nri_nonevents: nonevent_component,
    nri: event_component + nonevent_component,
  }
}

///|
/// Compute a score distribution summary.
pub fn survival_score_summary(
  scores : Array[SurvivalRiskScore],
) -> (Double, Double, Double, Double) {
  if scores.length() == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }
  let values : Array[Double] = []
  for score in scores {
    values.push(score.risk)
  }
  values.sort()
  (
    values[0],
    values[values.length() / 2],
    mean_doubles(values),
    values[values.length() - 1],
  )
}

///|
/// Rank risk scores with deterministic row ids for queue generation.
pub fn survival_rank_scores(
  scores : Array[SurvivalRiskScore],
) -> Array[SurvivalRiskScore] {
  let result : Array[SurvivalRiskScore] = []
  for score in scores {
    result.push(score)
  }
  result.sort_by(fn(a, b) {
    if a.risk > b.risk {
      -1
    } else if a.risk < b.risk {
      1
    } else {
      a.row_id.compare(b.row_id)
    }
  })
  result
}

///|
/// Select the top fraction of risk scores for an intervention queue.
pub fn survival_top_risk(
  scores : Array[SurvivalRiskScore],
  fraction : Double,
) -> Array[SurvivalRiskScore] {
  let ranked = survival_rank_scores(scores)
  let count = (ranked.length().to_double() * clamp_probability(fraction))
    .ceil()
    .to_int()
  let result : Array[SurvivalRiskScore] = []
  for i in 0..<count {
    if i < ranked.length() {
      result.push(ranked[i])
    }
  }
  result
}

///|
/// Export risk bands.
pub fn survival_risk_bands_to_csv(bands : Array[SurvivalRiskBand]) -> String {
  let builder = StringBuilder::new()
  builder.write_string("lower,upper,count,events,event_rate,mean_risk\n")
  for band in bands {
    builder.write_string(band.lower.to_string())
    builder.write_char(',')
    builder.write_string(band.upper.to_string())
    builder.write_char(',')
    builder.write_string(band.count.to_string())
    builder.write_char(',')
    builder.write_string(band.events.to_string())
    builder.write_char(',')
    builder.write_string(band.event_rate.to_string())
    builder.write_char(',')
    builder.write_string(band.mean_risk.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export classification metrics.
pub fn survival_classification_to_csv(
  metrics : Array[SurvivalClassificationMetrics],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "threshold,sensitivity,specificity,ppv,npv,net_benefit,sample_size\n",
  )
  for metric in metrics {
    builder.write_string(metric.threshold.to_string())
    builder.write_char(',')
    builder.write_string(metric.sensitivity.to_string())
    builder.write_char(',')
    builder.write_string(metric.specificity.to_string())
    builder.write_char(',')
    builder.write_string(metric.positive_predictive_value.to_string())
    builder.write_char(',')
    builder.write_string(metric.negative_predictive_value.to_string())
    builder.write_char(',')
    builder.write_string(metric.net_benefit.to_string())
    builder.write_char(',')
    builder.write_string(metric.sample_size.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate a risk-score vector.
pub fn survival_risk_scores_is_valid(scores : Array[SurvivalRiskScore]) -> Bool {
  for score in scores {
    if score.risk < 0.0 ||
      score.risk > 1.0 ||
      score.survival < 0.0 ||
      score.survival > 1.0 ||
      score.horizon < 0.0 {
      return false
    }
  }
  true
}
