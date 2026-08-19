///|
/// A generic model comparison row.
pub(all) struct SurvivalModelComparisonRow {
  name : String
  criterion : Double
  delta : Double
  weight : Double
  rank : Int
  comparison_row : Bool
} derive(Eq, Debug)

///|
/// Compare named models with a lower-is-better criterion.
pub fn survival_compare_models(
  names : Array[String],
  criteria : Array[Double],
) -> Array[SurvivalModelComparisonRow] {
  let count = if names.length() < criteria.length() {
    names.length()
  } else {
    criteria.length()
  }
  let result : Array[SurvivalModelComparisonRow] = []
  if count == 0 {
    return result
  }
  let minimum = survival_ext_min(criteria)
  let mut denominator = 0.0
  for i in 0..<count {
    denominator += double_exp(-0.5 * (criteria[i] - minimum))
  }
  for i in 0..<count {
    result.push({
      name: names[i],
      criterion: criteria[i],
      delta: criteria[i] - minimum,
      weight: survival_ext_safe_ratio(
        double_exp(-0.5 * (criteria[i] - minimum)),
        denominator,
      ),
      rank: 0,
      comparison_row: true,
    })
  }
  result.sort_by(fn(a, b) {
    if a.criterion < b.criterion {
      -1
    } else if a.criterion > b.criterion {
      1
    } else {
      a.name.compare(b.name)
    }
  })
  for i in 0..<result.length() {
    let row = result[i]
    result[i] = {
      name: row.name,
      criterion: row.criterion,
      delta: row.delta,
      weight: row.weight,
      rank: i + 1,
      comparison_row: true,
    }
  }
  result
}

///|
/// Select a model name by minimum criterion.
pub fn survival_select_model_name(
  rows : Array[SurvivalModelComparisonRow],
) -> String {
  if rows.length() == 0 {
    ""
  } else {
    rows[0].name
  }
}

///|
/// Compute a model-averaged criterion.
pub fn survival_model_average(
  rows : Array[SurvivalModelComparisonRow],
) -> Double {
  let mut total = 0.0
  for row in rows {
    total += row.weight * row.criterion
  }
  total
}

///|
/// Return models within a delta threshold.
pub fn survival_models_within_delta(
  rows : Array[SurvivalModelComparisonRow],
  delta : Double,
) -> Array[SurvivalModelComparisonRow] {
  let result : Array[SurvivalModelComparisonRow] = []
  for row in rows {
    if row.delta <= delta {
      result.push(row)
    }
  }
  result
}

///|
/// Compute the evidence ratio of two model weights.
pub fn survival_model_evidence_ratio(
  first : SurvivalModelComparisonRow,
  second : SurvivalModelComparisonRow,
) -> Double {
  survival_ext_safe_ratio(first.weight, second.weight)
}

///|
/// Export model comparisons.
pub fn survival_model_comparison_rows_to_csv(
  rows : Array[SurvivalModelComparisonRow],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("name,criterion,delta,weight,rank\n")
  for row in rows {
    builder.write_string(row.name)
    builder.write_char(',')
    builder.write_string(row.criterion.to_string())
    builder.write_char(',')
    builder.write_string(row.delta.to_string())
    builder.write_char(',')
    builder.write_string(row.weight.to_string())
    builder.write_char(',')
    builder.write_string(row.rank.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate comparison weights and ranks.
pub fn survival_model_comparisons_are_valid(
  rows : Array[SurvivalModelComparisonRow],
) -> Bool {
  let mut total = 0.0
  for row in rows {
    if row.rank <= 0 || row.criterion.is_nan() || row.weight < 0.0 {
      return false
    }
    total += row.weight
  }
  (total - 1.0).abs() < 0.000001 || rows.length() == 0
}
