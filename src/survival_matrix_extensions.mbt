///|
/// Validate a rectangular numeric matrix.
pub fn survival_matrix_is_rectangular(matrix : Array[Array[Double]]) -> Bool {
  if matrix.length() == 0 {
    return true
  }
  let width = matrix[0].length()
  for row in matrix {
    if row.length() != width {
      return false
    }
    for value in row {
      if value.is_nan() || value.is_inf() {
        return false
      }
    }
  }
  true
}

///|
/// Return a matrix transpose.
pub fn survival_matrix_transpose(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  if !survival_matrix_is_rectangular(matrix) || matrix.length() == 0 {
    return []
  }
  let width = matrix[0].length()
  let result : Array[Array[Double]] = []
  for column in 0..<width {
    let row : Array[Double] = []
    for source in matrix {
      row.push(source[column])
    }
    result.push(row)
  }
  result
}

///|
/// Matrix multiplication.
pub fn survival_matrix_multiply(
  left : Array[Array[Double]],
  right : Array[Array[Double]],
) -> Array[Array[Double]] {
  if !survival_matrix_is_rectangular(left) ||
    !survival_matrix_is_rectangular(right) ||
    left.length() == 0 ||
    right.length() == 0 ||
    left[0].length() != right.length() {
    return []
  }
  let result : Array[Array[Double]] = []
  for i in 0..<left.length() {
    let row : Array[Double] = []
    for j in 0..<right[0].length() {
      let mut value = 0.0
      for k in 0..<right.length() {
        value += left[i][k] * right[k][j]
      }
      row.push(value)
    }
    result.push(row)
  }
  result
}

///|
/// Matrix-vector multiplication.
pub fn survival_matrix_vector_multiply(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  if !survival_matrix_is_rectangular(matrix) ||
    matrix.length() == 0 ||
    matrix[0].length() != vector.length() {
    return []
  }
  let result : Array[Double] = []
  for row in matrix {
    result.push(dot_product(row, vector))
  }
  result
}

///|
/// Add two matrices.
pub fn survival_matrix_add(
  left : Array[Array[Double]],
  right : Array[Array[Double]],
) -> Array[Array[Double]] {
  if !survival_matrix_is_rectangular(left) ||
    !survival_matrix_is_rectangular(right) ||
    left.length() != right.length() ||
    left.length() == 0 ||
    left[0].length() != right[0].length() {
    return []
  }
  let result : Array[Array[Double]] = []
  for i in 0..<left.length() {
    let row : Array[Double] = []
    for j in 0..<left[i].length() {
      row.push(left[i][j] + right[i][j])
    }
    result.push(row)
  }
  result
}

///|
/// Scale a matrix.
pub fn survival_matrix_scale(
  matrix : Array[Array[Double]],
  factor : Double,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for source in matrix {
    let row : Array[Double] = []
    for value in source {
      row.push(value * factor)
    }
    result.push(row)
  }
  result
}

///|
/// Construct an identity matrix.
pub fn survival_matrix_identity(size : Int) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  if size <= 0 {
    return result
  }
  for i in 0..<size {
    let row = Array::make(size, 0.0)
    row[i] = 1.0
    result.push(row)
  }
  result
}

///|
/// Return the diagonal.
pub fn survival_matrix_diagonal(matrix : Array[Array[Double]]) -> Array[Double] {
  let result : Array[Double] = []
  let size = if matrix.length() <
    (if matrix.length() == 0 { 0 } else { matrix[0].length() }) {
    matrix.length()
  } else if matrix.length() == 0 {
    0
  } else {
    matrix[0].length()
  }
  for i in 0..<size {
    result.push(matrix[i][i])
  }
  result
}

///|
/// Matrix trace.
pub fn survival_matrix_trace(matrix : Array[Array[Double]]) -> Double {
  survival_ext_sum(survival_matrix_diagonal(matrix))
}

///|
/// Frobenius norm.
pub fn survival_matrix_frobenius(matrix : Array[Array[Double]]) -> Double {
  let mut total = 0.0
  for row in matrix {
    for value in row {
      total += value * value
    }
  }
  double_sqrt(total)
}

///|
/// Symmetrize a square matrix.
pub fn survival_matrix_symmetrize(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  if !survival_matrix_is_rectangular(matrix) ||
    matrix.length() == 0 ||
    matrix.length() != matrix[0].length() {
    return []
  }
  let result : Array[Array[Double]] = []
  for i in 0..<matrix.length() {
    let row : Array[Double] = []
    for j in 0..<matrix.length() {
      row.push((matrix[i][j] + matrix[j][i]) / 2.0)
    }
    result.push(row)
  }
  result
}

///|
/// Check symmetry within tolerance.
pub fn survival_matrix_is_symmetric(
  matrix : Array[Array[Double]],
  tolerance : Double,
) -> Bool {
  if !survival_matrix_is_rectangular(matrix) ||
    matrix.length() == 0 ||
    matrix.length() != matrix[0].length() {
    return false
  }
  for i in 0..<matrix.length() {
    for j in 0..<matrix.length() {
      if (matrix[i][j] - matrix[j][i]).abs() > tolerance {
        return false
      }
    }
  }
  true
}

///|
/// Return a positive diagonal ridge regularization.
pub fn survival_matrix_ridge(
  matrix : Array[Array[Double]],
  ridge : Double,
) -> Array[Array[Double]] {
  let result = survival_matrix_copy(matrix)
  let amount = if ridge <= 0.0 { 0.000001 } else { ridge }
  if result.length() == 0 || result.length() != result[0].length() {
    return result
  }
  for i in 0..<result.length() {
    result[i][i] += amount
  }
  result
}

///|
/// Defensive matrix copy.
pub fn survival_matrix_copy(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for source in matrix {
    let row : Array[Double] = []
    for value in source {
      row.push(value)
    }
    result.push(row)
  }
  result
}

///|
/// Return a matrix determinant for small square matrices.
pub fn survival_matrix_determinant(matrix : Array[Array[Double]]) -> Double {
  if !survival_matrix_is_rectangular(matrix) ||
    matrix.length() == 0 ||
    matrix.length() != matrix[0].length() {
    return 0.0
  }
  if matrix.length() == 1 {
    return matrix[0][0]
  }
  if matrix.length() == 2 {
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
  }
  let mut determinant = 0.0
  for column in 0..<matrix.length() {
    let minor : Array[Array[Double]] = []
    for i in 1..<matrix.length() {
      let row : Array[Double] = []
      for j in 0..<matrix.length() {
        if j != column {
          row.push(matrix[i][j])
        }
      }
      minor.push(row)
    }
    let sign = if column % 2 == 0 { 1.0 } else { -1.0 }
    determinant += sign * matrix[0][column] * survival_matrix_determinant(minor)
  }
  determinant
}

///|
/// Return condition-friendly inverse with ridge fallback.
pub fn survival_matrix_safe_inverse(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  match matrix_inverse(matrix) {
    Some(inverse) => inverse
    None => matrix_inverse_or_identity(survival_matrix_ridge(matrix, 0.000001))
  }
}

///|
/// Return inverse if available or an identity fallback.
fn matrix_inverse_or_identity(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  match matrix_inverse(matrix) {
    Some(value) => value
    None => survival_matrix_identity(matrix.length())
  }
}

///|
/// Compute a covariance-like crossproduct matrix.
pub fn survival_matrix_crossproduct(
  values : Array[Array[Double]],
) -> Array[Array[Double]] {
  survival_matrix_multiply(survival_matrix_transpose(values), values)
}

///|
/// Center matrix columns.
pub fn survival_matrix_center(
  values : Array[Array[Double]],
) -> Array[Array[Double]] {
  if values.length() == 0 || !survival_matrix_is_rectangular(values) {
    return []
  }
  let means : Array[Double] = []
  for column in 0..<values[0].length() {
    let column_values : Array[Double] = []
    for row in values {
      column_values.push(row[column])
    }
    means.push(survival_ext_mean(column_values))
  }
  let result : Array[Array[Double]] = []
  for row in values {
    let centered : Array[Double] = []
    for column in 0..<row.length() {
      centered.push(row[column] - means[column])
    }
    result.push(centered)
  }
  result
}

///|
/// Standardize matrix columns.
pub fn survival_matrix_standardize(
  values : Array[Array[Double]],
) -> Array[Array[Double]] {
  let centered = survival_matrix_center(values)
  if centered.length() == 0 {
    return centered
  }
  let result : Array[Array[Double]] = []
  for row in centered {
    let output : Array[Double] = []
    for column in 0..<row.length() {
      let column_values : Array[Double] = []
      for source in centered {
        column_values.push(source[column])
      }
      output.push(
        survival_ext_safe_ratio(
          row[column],
          survival_ext_standard_deviation(column_values),
        ),
      )
    }
    result.push(output)
  }
  result
}

///|
/// Export a matrix in CSV form.
pub fn survival_matrix_to_csv(matrix : Array[Array[Double]]) -> String {
  let builder = StringBuilder::new()
  for row in matrix {
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
