///|
/// A named term in a survival-model design matrix.
pub(all) struct SurvivalFormulaTerm {
  name : String
  source_column : Int
  level : Int
  interaction_column : Int
  interaction_term : Bool
} derive(Eq, Debug)

///|
/// A validated design specification.
pub(all) struct SurvivalDesignSpec {
  terms : Array[SurvivalFormulaTerm]
  include_intercept : Bool
  reference_levels : Array[Int]
  design_version : String
} derive(Eq, Debug)

///|
/// Design matrix output with row and column audit information.
pub(all) struct SurvivalDesignMatrix {
  values : Array[Array[Double]]
  terms : Array[SurvivalFormulaTerm]
  row_count : Int
  column_count : Int
} derive(Eq, Debug)

///|
/// A compact formula parsing error.
pub(all) enum FormulaError {
  EmptyFormula
  InvalidTerm
  InvalidColumn
  InvalidLevel
  DuplicateTerm
} derive(Eq, Debug)

///|
/// Create a numeric design specification from column count.
pub fn survival_design_spec(
  column_count : Int,
  include_intercept : Bool,
) -> SurvivalDesignSpec {
  let terms : Array[SurvivalFormulaTerm] = []
  if include_intercept {
    terms.push({
      name: "(Intercept)",
      source_column: -1,
      level: 0,
      interaction_column: -1,
      interaction_term: false,
    })
  }
  for column in 0..<column_count {
    terms.push({
      name: "x" + column.to_string(),
      source_column: column,
      level: 0,
      interaction_column: -1,
      interaction_term: false,
    })
  }
  let references = Array::make(column_count, 0)
  {
    terms,
    include_intercept,
    reference_levels: references,
    design_version: "1",
  }
}

///|
/// Add a single interaction term to a specification.
pub fn survival_design_add_interaction(
  spec : SurvivalDesignSpec,
  left : Int,
  right : Int,
) -> Result[SurvivalDesignSpec, FormulaError] {
  if left < 0 || right < 0 {
    return Err(InvalidColumn)
  }
  let terms : Array[SurvivalFormulaTerm] = []
  for term in spec.terms {
    terms.push(term)
  }
  terms.push({
    name: "x" + left.to_string() + ":x" + right.to_string(),
    source_column: left,
    level: right,
    interaction_column: right,
    interaction_term: true,
  })
  Ok({
    terms,
    include_intercept: spec.include_intercept,
    reference_levels: spec.reference_levels,
    design_version: spec.design_version,
  })
}

///|
/// Add one-hot terms for a categorical column.
pub fn survival_design_add_factor(
  spec : SurvivalDesignSpec,
  column : Int,
  level_count : Int,
  reference : Int,
) -> Result[SurvivalDesignSpec, FormulaError] {
  if column < 0 || level_count <= 1 || reference < 0 || reference >= level_count {
    return Err(InvalidLevel)
  }
  let terms : Array[SurvivalFormulaTerm] = []
  for term in spec.terms {
    if term.source_column != column || term.level == 0 {
      terms.push(term)
    }
  }
  for level in 0..<level_count {
    if level != reference {
      terms.push({
        name: "x" + column.to_string() + "[" + level.to_string() + "]",
        source_column: column,
        level,
        interaction_column: -1,
        interaction_term: false,
      })
    }
  }
  let references : Array[Int] = []
  for value in spec.reference_levels {
    references.push(value)
  }
  while references.length() <= column {
    references.push(0)
  }
  references[column] = reference
  Ok({
    terms,
    include_intercept: spec.include_intercept,
    reference_levels: references,
    design_version: spec.design_version,
  })
}

///|
/// Evaluate a numeric design specification.
pub fn survival_design_build(
  spec : SurvivalDesignSpec,
  values : Array[Array[Double]],
) -> SurvivalDesignMatrix {
  let output : Array[Array[Double]] = []
  for row in values {
    let design_row : Array[Double] = []
    for term in spec.terms {
      if term.source_column < 0 {
        design_row.push(1.0)
      } else if term.interaction_term {
        if term.source_column < row.length() &&
          term.interaction_column < row.length() {
          design_row.push(
            row[term.source_column] * row[term.interaction_column],
          )
        } else {
          design_row.push(0.0)
        }
      } else if term.level > 0 {
        let value = if term.source_column < row.length() {
          row[term.source_column]
        } else {
          0.0
        }
        design_row.push(if value == term.level.to_double() { 1.0 } else { 0.0 })
      } else if term.source_column < row.length() {
        design_row.push(row[term.source_column])
      } else {
        design_row.push(0.0)
      }
    }
    output.push(design_row)
  }
  {
    values: output,
    terms: spec.terms,
    row_count: output.length(),
    column_count: spec.terms.length(),
  }
}

///|
/// Center and scale a design matrix column-wise.
pub fn survival_design_standardize(
  matrix : SurvivalDesignMatrix,
) -> SurvivalDesignMatrix {
  let means : Array[Double] = []
  let scales : Array[Double] = []
  for column in 0..<matrix.column_count {
    let column_values : Array[Double] = []
    for row in matrix.values {
      column_values.push(row[column])
    }
    let mean = if column_values.length() == 0 {
      0.0
    } else {
      mean_doubles(column_values)
    }
    let mut variance = 0.0
    for value in column_values {
      let delta = value - mean
      variance += delta * delta
    }
    means.push(mean)
    scales.push(
      double_sqrt(variance / (column_values.length() + 1).to_double()),
    )
  }
  let output : Array[Array[Double]] = []
  for row in matrix.values {
    let scaled : Array[Double] = []
    for column in 0..<matrix.column_count {
      scaled.push(
        (row[column] - means[column]) / clamp_positive(scales[column]),
      )
    }
    output.push(scaled)
  }
  {
    values: output,
    terms: matrix.terms,
    row_count: matrix.row_count,
    column_count: matrix.column_count,
  }
}

///|
/// Return a detached column from a design matrix.
pub fn SurvivalDesignMatrix::column(
  self : SurvivalDesignMatrix,
  index : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  if index < 0 || index >= self.column_count {
    return result
  }
  for row in self.values {
    result.push(row[index])
  }
  result
}

///|
/// Return the linear predictor for a design row.
pub fn SurvivalDesignMatrix::linear_predictor(
  self : SurvivalDesignMatrix,
  row : Array[Double],
  coefficients : Array[Double],
) -> Double {
  if row.length() != coefficients.length() || row.length() != self.column_count {
    return 0.0
  }
  dot_product(row, coefficients)
}

///|
/// Add two design matrices with identical terms.
pub fn survival_design_append(
  left : SurvivalDesignMatrix,
  right : SurvivalDesignMatrix,
) -> SurvivalDesignMatrix {
  let values : Array[Array[Double]] = []
  for row in left.values {
    values.push(row)
  }
  for row in right.values {
    values.push(row)
  }
  {
    values,
    terms: left.terms,
    row_count: values.length(),
    column_count: left.column_count,
  }
}

///|
/// Remove rows containing non-finite design values.
pub fn survival_design_complete_cases(
  matrix : SurvivalDesignMatrix,
) -> SurvivalDesignMatrix {
  let values : Array[Array[Double]] = []
  for row in matrix.values {
    let mut complete = true
    for value in row {
      if value.is_nan() || value.is_inf() {
        complete = false
      }
    }
    if complete {
      values.push(row)
    }
  }
  {
    values,
    terms: matrix.terms,
    row_count: values.length(),
    column_count: matrix.column_count,
  }
}

///|
/// Test whether a specification has unique term names.
pub fn survival_design_is_valid(spec : SurvivalDesignSpec) -> Bool {
  if spec.terms.length() == 0 {
    return false
  }
  for i in 0..<spec.terms.length() {
    for j in (i + 1)..<spec.terms.length() {
      if spec.terms[i].name == spec.terms[j].name {
        return false
      }
    }
  }
  true
}

///|
/// Export design metadata.
pub fn SurvivalDesignMatrix::to_csv(self : SurvivalDesignMatrix) -> String {
  let builder = StringBuilder::new()
  for index in 0..<self.terms.length() {
    if index > 0 {
      builder.write_char(',')
    }
    builder.write_string(self.terms[index].name)
  }
  builder.write_char('\n')
  for row in self.values {
    for index in 0..<row.length() {
      if index > 0 {
        builder.write_char(',')
      }
      builder.write_string(row[index].to_string())
    }
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Create an interaction-only matrix for two columns.
pub fn survival_interaction_column(
  values : Array[Array[Double]],
  left : Int,
  right : Int,
) -> Array[Double] {
  let result : Array[Double] = []
  for row in values {
    if left < row.length() && right < row.length() {
      result.push(row[left] * row[right])
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Append one numeric column to a row matrix.
pub fn survival_design_append_column(
  values : Array[Array[Double]],
  column : Array[Double],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for i in 0..<values.length() {
    let row : Array[Double] = []
    for value in values[i] {
      row.push(value)
    }
    if i < column.length() {
      row.push(column[i])
    }
    result.push(row)
  }
  result
}

///|
/// Export formula terms.
pub fn survival_formula_terms_to_csv(
  terms : Array[SurvivalFormulaTerm],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string(
    "term,source_column,level,interaction_column,is_interaction\n",
  )
  for term in terms {
    builder.write_string(term.name)
    builder.write_char(',')
    builder.write_string(term.source_column.to_string())
    builder.write_char(',')
    builder.write_string(term.level.to_string())
    builder.write_char(',')
    builder.write_string(term.interaction_column.to_string())
    builder.write_char(',')
    builder.write_string(term.interaction_term.to_string())
    builder.write_char('\n')
  }
  builder.to_string()
}
