///|
/// A clinically interpretable risk stratification table.
///
/// The implementation uses fixed score cut points and a common prediction
/// horizon. It is intentionally small and deterministic so that the result
/// can be embedded in reports, dashboards, or downstream validation jobs.
pub(all) struct SurvivalRiskStratum {
  name : String
  lower : Double
  upper : Double
  count : Int
  events : Int
  observed_risk : Double
  predicted_risk : Double
  stratum_row : Int
} derive(Eq, Debug)

///|
fn survival_risk_stratum_name(
  row : Int,
  lower : Double,
  upper : Double,
) -> String {
  if row == 0 {
    "low [" + lower.to_string() + ", " + upper.to_string() + ")"
  } else {
    "stratum_" +
    row.to_string() +
    " [" +
    lower.to_string() +
    ", " +
    upper.to_string() +
    ")"
  }
}

///|
/// Build observed-versus-predicted risk summaries for score strata.
pub fn survival_risk_stratify(
  risks : Array[Double],
  data : Array[SurvivalData],
  cuts : Array[Double],
  horizon : Double,
) -> Array[SurvivalRiskStratum] {
  if risks.length() != data.length() || cuts.length() == 0 || horizon <= 0.0 {
    return []
  }
  let sorted = survival_ext_sorted(cuts)
  let strata = Array::make(sorted.length() + 1, 0)
  let counts = Array::make(sorted.length() + 1, 0)
  let events = Array::make(sorted.length() + 1, 0)
  let risk_sum = Array::make(sorted.length() + 1, 0.0)
  for i = 0; i < risks.length(); i = i + 1 {
    let score = risks[i]
    let mut group = 0
    while group < sorted.length() && score >= sorted[group] {
      group = group + 1
    }
    counts[group] = counts[group] + 1
    let row = data[i].observations[0]
    let event_at_horizon = row.event == 1 && row.time <= horizon
    if event_at_horizon {
      events[group] = events[group] + 1
    }
    risk_sum[group] = risk_sum[group] + score
    strata[group] = strata[group] + 1
  }
  let output = Array::make(sorted.length() + 1, SurvivalRiskStratum::{
    name: "",
    lower: 0.0,
    upper: 0.0,
    count: 0,
    events: 0,
    observed_risk: 0.0,
    predicted_risk: 0.0,
    stratum_row: 0,
  })
  for group = 0; group < output.length(); group = group + 1 {
    let lower = if group == 0 { 0.0 } else { sorted[group - 1] }
    let upper = if group < sorted.length() { sorted[group] } else { 1.0 }
    let n = counts[group]
    let observed = if n == 0 {
      0.0
    } else {
      events[group].to_double() / n.to_double()
    }
    let predicted = if n == 0 { 0.0 } else { risk_sum[group] / n.to_double() }
    output[group] = SurvivalRiskStratum::{
      name: survival_risk_stratum_name(group, lower, upper),
      lower,
      upper,
      count: n,
      events: events[group],
      observed_risk: observed,
      predicted_risk: predicted,
      stratum_row: group,
    }
  }
  output
}

///|
/// Sum of absolute observed-predicted risk differences across strata.
pub fn survival_risk_stratification_error(
  strata : Array[SurvivalRiskStratum],
) -> Double {
  let mut total = 0.0
  for row in strata {
    total = total + (row.observed_risk - row.predicted_risk).abs()
  }
  total
}

///|
/// Return the stratum with the largest observed event risk.
pub fn survival_risk_worst_stratum(
  strata : Array[SurvivalRiskStratum],
) -> SurvivalRiskStratum? {
  if strata.length() == 0 {
    return None
  }
  let mut best = strata[0]
  for row in strata[1:] {
    if row.observed_risk > best.observed_risk {
      best = row
    }
  }
  Some(best)
}

///|
pub fn survival_risk_strata_to_csv(
  strata : Array[SurvivalRiskStratum],
) -> String {
  let mut out = "name,lower,upper,count,events,observed_risk,predicted_risk\n"
  for row in strata {
    out = out +
      row.name +
      "," +
      row.lower.to_string() +
      "," +
      row.upper.to_string() +
      "," +
      row.count.to_string() +
      "," +
      row.events.to_string() +
      "," +
      row.observed_risk.to_string() +
      "," +
      row.predicted_risk.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_risk_strata_is_valid(
  strata : Array[SurvivalRiskStratum],
) -> Bool {
  for row in strata {
    if row.count < 0 ||
      row.events < 0 ||
      row.events > row.count ||
      row.lower > row.upper ||
      row.observed_risk < 0.0 ||
      row.observed_risk > 1.0 {
      return false
    }
  }
  true
}
