///|
/// Group-level Kaplan-Meier operations for cohort reporting.
pub(all) struct SurvivalGroupCurve {
  group_row : Int
  label : String
  sample_size : Int
  events : Int
  median : Double
  rmst : Double
  horizon : Double
  survival_at_horizon : Double
  curve : KMResult
} derive(Eq, Debug)

///|
pub(all) struct SurvivalGroupKMContrast {
  left_row : Int
  right_row : Int
  horizon : Double
  left_survival : Double
  right_survival : Double
  difference : Double
  ratio : Double
  left_rmst : Double
  right_rmst : Double
  rmst_difference : Double
} derive(Eq, Debug)

///|
fn survival_group_median(curve : KMResult) -> Double {
  match curve.median() {
    Some(value) => value
    None => -1.0
  }
}

///|
/// Estimate one KM curve per labeled group.
pub fn survival_group_km(
  groups : Array[SurvivalData],
  labels : Array[String],
  horizon : Double,
) -> Array[SurvivalGroupCurve] {
  let output = []
  for i = 0; i < groups.length(); i = i + 1 {
    let label = if i < labels.length() {
      labels[i]
    } else {
      "group_" + i.to_string()
    }
    let curve = kaplan_meier(groups[i])
    output.push({
      group_row: i,
      label,
      sample_size: groups[i].length(),
      events: groups[i].event_count(),
      median: survival_group_median(curve),
      rmst: curve.restricted_mean(horizon),
      horizon,
      survival_at_horizon: curve.predict(horizon),
      curve,
    })
  }
  output
}

///|
/// Compare every pair of group curves at a common horizon.
pub fn survival_group_km_contrasts(
  curves : Array[SurvivalGroupCurve],
) -> Array[SurvivalGroupKMContrast] {
  let output = []
  for left = 0; left < curves.length(); left = left + 1 {
    for right = left + 1; right < curves.length(); right = right + 1 {
      let a = curves[left]
      let b = curves[right]
      let difference = a.survival_at_horizon - b.survival_at_horizon
      let ratio = if b.survival_at_horizon == 0.0 {
        0.0
      } else {
        a.survival_at_horizon / b.survival_at_horizon
      }
      output.push({
        left_row: a.group_row,
        right_row: b.group_row,
        horizon: a.horizon,
        left_survival: a.survival_at_horizon,
        right_survival: b.survival_at_horizon,
        difference,
        ratio,
        left_rmst: a.rmst,
        right_rmst: b.rmst,
        rmst_difference: a.rmst - b.rmst,
      })
    }
  }
  output
}

///|
/// Select the group with the greatest survival at the requested horizon.
pub fn survival_group_best(
  curves : Array[SurvivalGroupCurve],
) -> SurvivalGroupCurve? {
  if curves.length() == 0 {
    return None
  }
  let mut best = curves[0]
  for curve in curves[1:] {
    if curve.survival_at_horizon > best.survival_at_horizon {
      best = curve
    }
  }
  Some(best)
}

///|
/// Select the group with the lowest survival at the requested horizon.
pub fn survival_group_worst(
  curves : Array[SurvivalGroupCurve],
) -> SurvivalGroupCurve? {
  if curves.length() == 0 {
    return None
  }
  let mut worst = curves[0]
  for curve in curves[1:] {
    if curve.survival_at_horizon < worst.survival_at_horizon {
      worst = curve
    }
  }
  Some(worst)
}

///|
pub fn survival_group_curves_to_csv(
  curves : Array[SurvivalGroupCurve],
) -> String {
  let mut out = "group_row,label,sample_size,events,median,rmst,horizon,survival_at_horizon\n"
  for row in curves {
    out = out +
      row.group_row.to_string() +
      "," +
      row.label +
      "," +
      row.sample_size.to_string() +
      "," +
      row.events.to_string() +
      "," +
      row.median.to_string() +
      "," +
      row.rmst.to_string() +
      "," +
      row.horizon.to_string() +
      "," +
      row.survival_at_horizon.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_group_km_contrasts_to_csv(
  contrasts : Array[SurvivalGroupKMContrast],
) -> String {
  let mut out = "left_row,right_row,horizon,left_survival,right_survival,difference,ratio,left_rmst,right_rmst,rmst_difference\n"
  for row in contrasts {
    out = out +
      row.left_row.to_string() +
      "," +
      row.right_row.to_string() +
      "," +
      row.horizon.to_string() +
      "," +
      row.left_survival.to_string() +
      "," +
      row.right_survival.to_string() +
      "," +
      row.difference.to_string() +
      "," +
      row.ratio.to_string() +
      "," +
      row.left_rmst.to_string() +
      "," +
      row.right_rmst.to_string() +
      "," +
      row.rmst_difference.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_group_curves_are_valid(
  curves : Array[SurvivalGroupCurve],
) -> Bool {
  for curve in curves {
    if curve.sample_size < 0 ||
      curve.events < 0 ||
      curve.events > curve.sample_size ||
      curve.horizon <= 0.0 ||
      curve.survival_at_horizon < 0.0 ||
      curve.survival_at_horizon > 1.0 ||
      curve.rmst < 0.0 {
      return false
    }
  }
  true
}
