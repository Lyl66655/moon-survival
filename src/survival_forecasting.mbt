///|
/// A forecast point with uncertainty and scenario metadata.
pub(all) struct SurvivalForecastPoint {
  horizon : Double
  survival : Double
  risk : Double
  lower : Double
  upper : Double
  scenario : String
} derive(Eq, Debug)

///|
/// A future forecast profile.
pub(all) struct SurvivalForecastProfile {
  points : Array[SurvivalForecastPoint]
  last_observed_time : Double
  forecast_horizon : Double
  forecast_row : Bool
} derive(Eq, Debug)

///|
/// Forecast a KM curve with a flat-tail scenario.
pub fn survival_forecast_km(
  data : SurvivalData,
  horizons : Array[Double],
  scenario : String,
) -> SurvivalForecastProfile {
  let fit = kaplan_meier(data)
  let last_time = if data.length() == 0 {
    0.0
  } else {
    data.get(data.length() - 1).time
  }
  let points : Array[SurvivalForecastPoint] = []
  let maximum = if horizons.length() == 0 {
    0.0
  } else {
    survival_ext_max(horizons)
  }
  for horizon in horizons {
    let survival = fit.survival_at(horizon)
    let variance = match fit.step_at(horizon) {
      Some(step) => step.std_err
      None => 0.0
    }
    points.push({
      horizon,
      survival,
      risk: 1.0 - survival,
      lower: clamp_probability(survival - 1.96 * variance),
      upper: clamp_probability(survival + 1.96 * variance),
      scenario,
    })
  }
  {
    points,
    last_observed_time: last_time,
    forecast_horizon: maximum,
    forecast_row: true,
  }
}

///|
/// Extrapolate with a constant hazard after the last observed time.
pub fn survival_forecast_constant_hazard(
  data : SurvivalData,
  horizons : Array[Double],
  hazard : Double,
) -> SurvivalForecastProfile {
  let fit = kaplan_meier(data)
  let last_time = if data.length() == 0 {
    0.0
  } else {
    data.get(data.length() - 1).time
  }
  let base_survival = fit.survival_at(last_time)
  let points : Array[SurvivalForecastPoint] = []
  for horizon in horizons {
    let survival = if horizon <= last_time {
      fit.survival_at(horizon)
    } else {
      base_survival *
      double_exp(-clamp_positive(hazard) * (horizon - last_time))
    }
    points.push({
      horizon,
      survival: clamp_probability(survival),
      risk: clamp_probability(1.0 - survival),
      lower: clamp_probability(survival - 0.05),
      upper: clamp_probability(survival + 0.05),
      scenario: "constant_hazard",
    })
  }
  {
    points,
    last_observed_time: last_time,
    forecast_horizon: if horizons.length() == 0 {
      0.0
    } else {
      survival_ext_max(horizons)
    },
    forecast_row: true,
  }
}

///|
/// Extrapolate a parametric fit to a grid.
pub fn survival_forecast_parametric(
  fit : ParametricFit,
  horizons : Array[Double],
  _covariates : Array[Double],
) -> SurvivalForecastProfile {
  let points : Array[SurvivalForecastPoint] = []
  for horizon in horizons {
    let survival = fit.survival_at(horizon)
    points.push({
      horizon,
      survival,
      risk: 1.0 - survival,
      lower: clamp_probability(survival - 0.05),
      upper: clamp_probability(survival + 0.05),
      scenario: parametric_family_name(fit.family),
    })
  }
  {
    points,
    last_observed_time: 0.0,
    forecast_horizon: if horizons.length() == 0 {
      0.0
    } else {
      survival_ext_max(horizons)
    },
    forecast_row: true,
  }
}

///|
/// Compute a forecast horizon from a target survival.
pub fn survival_forecast_target_time(
  profile : SurvivalForecastProfile,
  target_survival : Double,
) -> Double? {
  for point in profile.points {
    if point.survival <= target_survival {
      return Some(point.horizon)
    }
  }
  None
}

///|
/// Compare two forecast profiles on a shared grid.
pub fn survival_forecast_difference(
  first : SurvivalForecastProfile,
  second : SurvivalForecastProfile,
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  let count = if first.points.length() < second.points.length() {
    first.points.length()
  } else {
    second.points.length()
  }
  for index in 0..<count {
    result.push({
      time: first.points[index].horizon,
      value: first.points[index].survival - second.points[index].survival,
    })
  }
  result
}

///|
/// Integrate a forecast profile.
pub fn survival_forecast_rmst(profile : SurvivalForecastProfile) -> Double {
  survival_curve_area(
    survival_forecast_to_curve(profile),
    profile.forecast_horizon,
  )
}

///|
/// Convert a profile to a generic survival curve.
pub fn survival_forecast_to_curve(
  profile : SurvivalForecastProfile,
) -> Array[CurvePoint] {
  let result : Array[CurvePoint] = []
  for point in profile.points {
    result.push({ time: point.horizon, value: point.survival })
  }
  result
}

///|
/// Serialize a forecast profile.
pub fn survival_forecast_to_csv(profile : SurvivalForecastProfile) -> String {
  let builder = StringBuilder::new()
  builder.write_string("horizon,survival,risk,lower,upper,scenario\n")
  for point in profile.points {
    builder.write_string(point.horizon.to_string())
    builder.write_char(',')
    builder.write_string(point.survival.to_string())
    builder.write_char(',')
    builder.write_string(point.risk.to_string())
    builder.write_char(',')
    builder.write_string(point.lower.to_string())
    builder.write_char(',')
    builder.write_string(point.upper.to_string())
    builder.write_char(',')
    builder.write_string(point.scenario)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate forecast values and ordering.
pub fn survival_forecast_is_valid(profile : SurvivalForecastProfile) -> Bool {
  for index in 0..<profile.points.length() {
    let point = profile.points[index]
    if point.horizon < 0.0 ||
      point.survival < 0.0 ||
      point.survival > 1.0 ||
      point.lower < 0.0 ||
      point.upper > 1.0 ||
      (index > 0 && point.horizon < profile.points[index - 1].horizon) {
      return false
    }
  }
  true
}
