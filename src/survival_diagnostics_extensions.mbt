///|
/// A residual summary used for model audit.
pub(all) struct SurvivalResidualSummary {
  count : Int
  mean : Double
  standard_deviation : Double
  minimum : Double
  maximum : Double
  positive_fraction : Double
} derive(Eq, Debug)

///|
/// A proportional-hazards diagnostic result.
pub(all) struct SurvivalProportionalHazardsCheck {
  covariate : Int
  correlation : Double
  statistic : Double
  p_value : Double
  violated : Bool
  diagnostic_row : Bool
} derive(Eq, Debug)

///|
/// Summarize residuals.
pub fn survival_residual_summary(
  residuals : Array[Double],
) -> SurvivalResidualSummary {
  let mut positive = 0
  for value in residuals {
    if value > 0.0 {
      positive += 1
    }
  }
  {
    count: residuals.length(),
    mean: survival_ext_mean(residuals),
    standard_deviation: survival_ext_standard_deviation(residuals),
    minimum: survival_ext_min(residuals),
    maximum: survival_ext_max(residuals),
    positive_fraction: if residuals.length() == 0 {
      0.0
    } else {
      positive.to_double() / residuals.length().to_double()
    },
  }
}

///|
/// Compute martingale residuals from an observed event and cumulative hazard.
pub fn survival_martingale_residuals(
  data : SurvivalData,
  cumulative_hazards : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  let count = if data.length() < cumulative_hazards.length() {
    data.length()
  } else {
    cumulative_hazards.length()
  }
  for i in 0..<count {
    result.push(data.get(i).event.to_double() - cumulative_hazards[i])
  }
  result
}

///|
/// Convert martingale residuals to deviance residuals.
pub fn survival_deviance_residuals(
  residuals : Array[Double],
  cumulative_hazards : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  let count = if residuals.length() < cumulative_hazards.length() {
    residuals.length()
  } else {
    cumulative_hazards.length()
  }
  for i in 0..<count {
    let hazard = clamp_positive(cumulative_hazards[i])
    let martingale = residuals[i]
    let deviance = if martingale >= 1.0 {
      0.0
    } else {
      2.0 *
      (
        martingale * double_log(clamp_positive((martingale + hazard) / hazard)) -
        (martingale + hazard)
      )
    }
    result.push(
      if deviance < 0.0 {
        -double_sqrt(-deviance)
      } else {
        double_sqrt(deviance)
      },
    )
  }
  result
}

///|
/// Compute a DFBETA-like influence approximation.
pub fn survival_dfbeta_approximation(
  residuals : Array[Double],
  covariate : Array[Double],
  information : Double,
) -> Array[Double] {
  let result : Array[Double] = []
  let denominator = clamp_positive(information)
  let count = if residuals.length() < covariate.length() {
    residuals.length()
  } else {
    covariate.length()
  }
  for i in 0..<count {
    result.push(residuals[i] * covariate[i] / denominator)
  }
  result
}

///|
/// Compute correlation between ordered time and a residual vector.
pub fn survival_time_residual_correlation(
  data : SurvivalData,
  residuals : Array[Double],
) -> Double {
  let count = if data.length() < residuals.length() {
    data.length()
  } else {
    residuals.length()
  }
  if count < 2 {
    return 0.0
  }
  let times : Array[Double] = []
  let values : Array[Double] = []
  for i in 0..<count {
    times.push(data.get(i).time)
    values.push(residuals[i])
  }
  let time_mean = survival_ext_mean(times)
  let value_mean = survival_ext_mean(values)
  let mut numerator = 0.0
  let mut time_variance = 0.0
  let mut value_variance = 0.0
  for i in 0..<count {
    let dt = times[i] - time_mean
    let dv = values[i] - value_mean
    numerator += dt * dv
    time_variance += dt * dt
    value_variance += dv * dv
  }
  survival_ext_safe_ratio(
    numerator,
    double_sqrt(time_variance * value_variance),
  )
}

///|
/// Build a proportional-hazards diagnostic from Schoenfeld-like residuals.
pub fn survival_ph_check(
  data : SurvivalData,
  residuals : Array[Array[Double]],
  alpha : Double,
) -> Array[SurvivalProportionalHazardsCheck] {
  if residuals.length() == 0 {
    return []
  }
  let width = residuals[0].length()
  let result : Array[SurvivalProportionalHazardsCheck] = []
  for covariate in 0..<width {
    let values : Array[Double] = []
    let times : Array[Double] = []
    for i in 0..<residuals.length() {
      if covariate < residuals[i].length() && i < data.length() {
        values.push(residuals[i][covariate])
        times.push(data.get(i).time)
      }
    }
    let correlation = survival_ext_safe_ratio(
      survival_ext_dot(
        survival_ext_add(
          times,
          survival_ext_scale(values, -survival_ext_mean(values)),
        ),
        survival_ext_ones(times.length()),
      ),
      survival_ext_norm(times),
    )
    let statistic = correlation * double_sqrt(values.length().to_double())
    let p_value = 2.0 * (1.0 - normal_cdf(statistic.abs()))
    result.push({
      covariate,
      correlation,
      statistic,
      p_value,
      violated: p_value < alpha,
      diagnostic_row: true,
    })
  }
  result
}

///|
/// Compute influence cutoff based on effective sample size.
pub fn survival_influence_cutoff(sample_size : Int) -> Double {
  if sample_size <= 0 {
    0.0
  } else {
    2.0 / sample_size.to_double()
  }
}

///|
/// Return indices exceeding an influence cutoff.
pub fn survival_influential_indices(
  values : Array[Double],
  cutoff : Double,
) -> Array[Int] {
  let result : Array[Int] = []
  for i in 0..<values.length() {
    if values[i].abs() > cutoff {
      result.push(i)
    }
  }
  result
}

///|
/// Export residual summary.
pub fn survival_residual_summary_to_csv(
  summary : SurvivalResidualSummary,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("count,")
  builder.write_string(summary.count.to_string())
  builder.write_char('\n')
  builder.write_string("mean,")
  builder.write_string(summary.mean.to_string())
  builder.write_char('\n')
  builder.write_string("standard_deviation,")
  builder.write_string(summary.standard_deviation.to_string())
  builder.write_char('\n')
  builder.write_string("minimum,")
  builder.write_string(summary.minimum.to_string())
  builder.write_char('\n')
  builder.write_string("maximum,")
  builder.write_string(summary.maximum.to_string())
  builder.write_char('\n')
  builder.write_string("positive_fraction,")
  builder.write_string(summary.positive_fraction.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Export proportional-hazards checks.
pub fn survival_ph_checks_to_csv(
  checks : Array[SurvivalProportionalHazardsCheck],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("covariate,correlation,statistic,p_value,violated\n")
  for check in checks {
    builder.write_string(check.covariate.to_string())
    builder.write_char(',')
    builder.write_string(check.correlation.to_string())
    builder.write_char(',')
    builder.write_string(check.statistic.to_string())
    builder.write_char(',')
    builder.write_string(check.p_value.to_string())
    builder.write_char(',')
    builder.write_string(check.violated.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
