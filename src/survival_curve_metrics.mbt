///|
/// Metrics for comparing and validating predicted survival curves.
pub(all) struct SurvivalMetricPoint {
  metric_row : Int
  time : Double
  observed : Double
  predicted : Double
  residual : Double
  absolute_error : Double
  squared_error : Double
} derive(Eq, Debug)

///|
pub(all) struct SurvivalCurveMetricSummary {
  count : Int
  mae : Double
  rmse : Double
  integrated_absolute_error : Double
  integrated_squared_error : Double
  max_absolute_error : Double
  concordance : Double
  metric_row : Int
} derive(Eq, Debug)

///|
/// Pair two curves at common time points and compute pointwise errors.
pub fn survival_curve_metric_points(
  times : Array[Double],
  observed : Array[Double],
  predicted : Array[Double],
) -> Array[SurvivalMetricPoint] {
  let output : Array[SurvivalMetricPoint] = []
  let count = times.length().min(observed.length()).min(predicted.length())
  for i = 0; i < count; i = i + 1 {
    let residual = predicted[i] - observed[i]
    let absolute_error = residual.abs()
    output.push({
      metric_row: i,
      time: times[i],
      observed: clamp_probability(observed[i]),
      predicted: clamp_probability(predicted[i]),
      residual,
      absolute_error,
      squared_error: residual * residual,
    })
  }
  output
}

///|
/// Integrate a metric by the trapezoidal rule over sorted time points.
pub fn survival_curve_integral(
  times : Array[Double],
  values : Array[Double],
) -> Double {
  let count = times.length().min(values.length())
  if count < 2 {
    return 0.0
  }
  let mut area = 0.0
  for i = 1; i < count; i = i + 1 {
    let width = (times[i] - times[i - 1]).max(0.0)
    area += width * (values[i] + values[i - 1]) / 2.0
  }
  area
}

///|
/// Summarize error across a pair of survival curves.
pub fn survival_curve_metric_summary(
  points : Array[SurvivalMetricPoint],
) -> SurvivalCurveMetricSummary {
  if points.length() == 0 {
    return {
      count: 0,
      mae: 0.0,
      rmse: 0.0,
      integrated_absolute_error: 0.0,
      integrated_squared_error: 0.0,
      max_absolute_error: 0.0,
      concordance: 0.0,
      metric_row: 0,
    }
  }
  let times : Array[Double] = []
  let absolute : Array[Double] = []
  let squared : Array[Double] = []
  let mut sum_absolute = 0.0
  let mut sum_squared = 0.0
  let mut max_absolute = 0.0
  for point in points {
    times.push(point.time)
    absolute.push(point.absolute_error)
    squared.push(point.squared_error)
    sum_absolute += point.absolute_error
    sum_squared += point.squared_error
    if point.absolute_error > max_absolute {
      max_absolute = point.absolute_error
    }
  }
  {
    count: points.length(),
    mae: sum_absolute / points.length().to_double(),
    rmse: double_sqrt(sum_squared / points.length().to_double()),
    integrated_absolute_error: survival_curve_integral(times, absolute),
    integrated_squared_error: survival_curve_integral(times, squared),
    max_absolute_error: max_absolute,
    concordance: survival_curve_concordance(points),
    metric_row: points.length(),
  }
}

///|
/// Fraction of ordered pairs where the predicted curve orders observations correctly.
pub fn survival_curve_concordance(
  points : Array[SurvivalMetricPoint],
) -> Double {
  let mut comparable = 0
  let mut concordant = 0
  for i = 0; i < points.length(); i = i + 1 {
    for j = i + 1; j < points.length(); j = j + 1 {
      let observed_order = points[i].observed - points[j].observed
      let predicted_order = points[i].predicted - points[j].predicted
      if observed_order.abs() > survival_epsilon {
        comparable += 1
        if observed_order * predicted_order > 0.0 {
          concordant += 1
        }
        if predicted_order.abs() <= survival_epsilon {
          concordant += 1
        }
      }
    }
  }
  if comparable == 0 {
    0.0
  } else {
    concordant.to_double() / comparable.to_double()
  }
}

///|
/// Add a constant calibration offset and clamp the resulting probability.
pub fn survival_curve_shift(
  predicted : Array[Double],
  offset : Double,
) -> Array[Double] {
  let output : Array[Double] = []
  for value in predicted {
    output.push(clamp_probability(value + offset))
  }
  output
}

///|
pub fn survival_curve_metrics_to_csv(
  points : Array[SurvivalMetricPoint],
) -> String {
  let mut out = "metric_row,time,observed,predicted,residual,absolute_error,squared_error\n"
  for point in points {
    out = out +
      point.metric_row.to_string() +
      "," +
      point.time.to_string() +
      "," +
      point.observed.to_string() +
      "," +
      point.predicted.to_string() +
      "," +
      point.residual.to_string() +
      "," +
      point.absolute_error.to_string() +
      "," +
      point.squared_error.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_curve_metric_points_are_valid(
  points : Array[SurvivalMetricPoint],
) -> Bool {
  for point in points {
    if point.time.is_nan() ||
      point.observed < 0.0 ||
      point.observed > 1.0 ||
      point.predicted < 0.0 ||
      point.predicted > 1.0 ||
      point.absolute_error < 0.0 ||
      point.squared_error < 0.0 {
      return false
    }
  }
  true
}

///|
pub fn survival_curve_metric_summary_is_valid(
  summary : SurvivalCurveMetricSummary,
) -> Bool {
  summary.count >= 0 &&
  summary.mae >= 0.0 &&
  summary.rmse >= 0.0 &&
  summary.integrated_absolute_error >= 0.0 &&
  summary.integrated_squared_error >= 0.0 &&
  summary.max_absolute_error >= 0.0 &&
  summary.concordance >= 0.0 &&
  summary.concordance <= 1.0
}
