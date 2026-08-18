///|
/// Covariance matrix for a rectangular observation-by-feature matrix.
pub fn covariance_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  if matrix.length() < 2 {
    return []
  }
  let means = column_means(matrix)
  let columns = means.length()
  let result = []
  for i in 0..<columns {
    let row = []
    for j in 0..<columns {
      let mut total = 0.0
      for observation in matrix {
        total += (observation[i] - means[i]) * (observation[j] - means[j])
      }
      row.push(total / (matrix.length() - 1).to_double())
    }
    result.push(row)
  }
  result
}

///|
/// Correlation matrix for a rectangular design matrix.
pub fn correlation_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let covariance = covariance_matrix(matrix)
  if covariance.length() == 0 {
    []
  } else {
    let scales = []
    for i in 0..<covariance.length() {
      scales.push(double_sqrt(covariance[i][i].abs()))
    }
    let result = []
    for i in 0..<covariance.length() {
      let row = []
      for j in 0..<covariance.length() {
        let denominator = scales[i] * scales[j]
        row.push(
          if denominator <= survival_epsilon {
            0.0
          } else {
            covariance[i][j] / denominator
          },
        )
      }
      result.push(row)
    }
    result
  }
}

///|
/// Variance inflation factor estimated from a correlation matrix inverse.
pub fn variance_inflation_factors(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let correlation = correlation_matrix(matrix)
  if correlation.length() == 0 {
    []
  } else {
    match matrix_inverse(matrix_ridge(correlation, 0.000000001)) {
      Some(inverse) => matrix_diagonal(inverse)
      None => {
        let result = []
        for _ in 0..<correlation.length() {
          result.push(0.0)
        }
        result
      }
    }
  }
}

///|
/// Standardize a covariate matrix and report its condition proxy.
pub(all) struct DesignDiagnostics {
  rows : Int
  columns : Int
  rank_proxy : Int
  condition_proxy : Double
  means : Array[Double]
  scales : Array[Double]
  variance_inflation : Array[Double]
} derive(Eq, Debug)

///|
/// Diagnose collinearity before Cox or regression fitting.
pub fn design_diagnostics(matrix : Array[Array[Double]]) -> DesignDiagnostics {
  if matrix.length() == 0 {
    return {
      rows: 0,
      columns: 0,
      rank_proxy: 0,
      condition_proxy: 0.0,
      means: [],
      scales: [],
      variance_inflation: [],
    }
  }
  let standardization = standardize_matrix(matrix)
  let scales = standardization.scales
  let mut minimum = 1.0e300
  let mut maximum = 0.0
  let mut rank = 0
  for scale in scales {
    if scale > survival_epsilon {
      rank += 1
      if scale < minimum {
        minimum = scale
      }
      if scale > maximum {
        maximum = scale
      }
    }
  }
  {
    rows: matrix.length(),
    columns: scales.length(),
    rank_proxy: rank,
    condition_proxy: if minimum <= survival_epsilon {
      0.0
    } else {
      maximum / minimum
    },
    means: standardization.means,
    scales,
    variance_inflation: variance_inflation_factors(matrix),
  }
}

///|
/// Return the upper-triangular correlation pairs above a threshold.
pub(all) struct CorrelationPair {
  left : Int
  right : Int
  correlation : Double
} derive(Eq, Debug)

///|
/// Find potentially collinear feature pairs.
pub fn high_correlation_pairs(
  matrix : Array[Array[Double]],
  threshold? : Double = 0.8,
) -> Array[CorrelationPair] {
  let correlation = correlation_matrix(matrix)
  let result = []
  for i in 0..<correlation.length() {
    for j in (i + 1)..<correlation.length() {
      if correlation[i][j].abs() >= threshold {
        result.push({ left: i, right: j, correlation: correlation[i][j] })
      }
    }
  }
  result
}
