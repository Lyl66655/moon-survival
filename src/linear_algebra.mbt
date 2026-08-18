///|
/// Matrix shape validation shared by regression and calibration routines.
fn matrix_is_square(matrix : Array[Array[Double]]) -> Bool {
  let size = matrix.length()
  for row in matrix {
    if row.length() != size {
      return false
    }
  }
  true
}

///|
/// Return a detached copy of a matrix.
pub fn matrix_copy(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let result = []
  for row in matrix {
    let copied = []
    for value in row {
      copied.push(value)
    }
    result.push(copied)
  }
  result
}

///|
/// Matrix transpose.
pub fn matrix_transpose(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let columns = matrix[0].length()
  for row in matrix {
    if row.length() != columns {
      panic()
    }
  }
  let result = []
  for j in 0..<columns {
    let row = []
    for i in 0..<matrix.length() {
      row.push(matrix[i][j])
    }
    result.push(row)
  }
  result
}

///|
/// Matrix multiplication.
pub fn matrix_multiply(
  left : Array[Array[Double]],
  right : Array[Array[Double]],
) -> Array[Array[Double]] {
  if left.length() == 0 || right.length() == 0 {
    return []
  }
  let right_rows = right.length()
  let right_columns = right[0].length()
  for row in right {
    if row.length() != right_columns {
      panic()
    }
  }
  for row in left {
    if row.length() != right_rows {
      panic()
    }
  }
  let result = []
  for i in 0..<left.length() {
    let row = []
    for j in 0..<right_columns {
      let mut value = 0.0
      for k in 0..<right_rows {
        value += left[i][k] * right[k][j]
      }
      row.push(value)
    }
    result.push(row)
  }
  result
}

///|
/// Multiply a matrix by a vector.
pub fn matrix_vector_multiply(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let result = []
  for row in matrix {
    if row.length() != vector.length() {
      panic()
    }
    result.push(dot_product(row, vector))
  }
  result
}

///|
/// Add matrices of equal shape.
pub fn matrix_add(
  left : Array[Array[Double]],
  right : Array[Array[Double]],
) -> Array[Array[Double]] {
  if left.length() != right.length() {
    panic()
  }
  let result = []
  for i in 0..<left.length() {
    if left[i].length() != right[i].length() {
      panic()
    }
    let row = []
    for j in 0..<left[i].length() {
      row.push(left[i][j] + right[i][j])
    }
    result.push(row)
  }
  result
}

///|
/// Scale a matrix.
pub fn matrix_scale(
  matrix : Array[Array[Double]],
  scale : Double,
) -> Array[Array[Double]] {
  let result = []
  for row in matrix {
    let output = []
    for value in row {
      output.push(value * scale)
    }
    result.push(output)
  }
  result
}

///|
/// Trace of a square matrix.
pub fn matrix_trace(matrix : Array[Array[Double]]) -> Double {
  if !matrix_is_square(matrix) {
    panic()
  }
  let mut trace = 0.0
  for i in 0..<matrix.length() {
    trace += matrix[i][i]
  }
  trace
}

///|
/// Solve a linear system with partial pivoting.
pub fn solve_linear_system(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double]? {
  let n = matrix.length()
  if n == 0 || !matrix_is_square(matrix) || vector.length() != n {
    return None
  }
  let augmented = []
  for i in 0..<n {
    let row = []
    for j in 0..<n {
      row.push(matrix[i][j])
    }
    row.push(vector[i])
    augmented.push(row)
  }
  for column in 0..<n {
    let mut pivot = column
    for row in (column + 1)..<n {
      if augmented[row][column].abs() > augmented[pivot][column].abs() {
        pivot = row
      }
    }
    if augmented[pivot][column].abs() <= survival_epsilon {
      return None
    }
    if pivot != column {
      let temporary = augmented[column]
      augmented[column] = augmented[pivot]
      augmented[pivot] = temporary
    }
    let pivot_value = augmented[column][column]
    for j in column..<(n + 1) {
      augmented[column][j] /= pivot_value
    }
    for row in 0..<n {
      if row != column {
        let multiplier = augmented[row][column]
        for j in column..<(n + 1) {
          augmented[row][j] -= multiplier * augmented[column][j]
        }
      }
    }
  }
  let result = []
  for i in 0..<n {
    result.push(augmented[i][n])
  }
  Some(result)
}

///|
/// Invert a square matrix using the same pivoted elimination routine.
pub fn matrix_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let n = matrix.length()
  if n == 0 || !matrix_is_square(matrix) {
    return None
  }
  let augmented = []
  for i in 0..<n {
    let row = []
    for j in 0..<n {
      row.push(matrix[i][j])
    }
    for j in 0..<n {
      row.push(if i == j { 1.0 } else { 0.0 })
    }
    augmented.push(row)
  }
  for column in 0..<n {
    let mut pivot = column
    for row in (column + 1)..<n {
      if augmented[row][column].abs() > augmented[pivot][column].abs() {
        pivot = row
      }
    }
    if augmented[pivot][column].abs() <= survival_epsilon {
      return None
    }
    if pivot != column {
      let temporary = augmented[column]
      augmented[column] = augmented[pivot]
      augmented[pivot] = temporary
    }
    let pivot_value = augmented[column][column]
    for j in 0..<(2 * n) {
      augmented[column][j] /= pivot_value
    }
    for row in 0..<n {
      if row != column {
        let multiplier = augmented[row][column]
        for j in 0..<(2 * n) {
          augmented[row][j] -= multiplier * augmented[column][j]
        }
      }
    }
  }
  let inverse = []
  for i in 0..<n {
    let row = []
    for j in 0..<n {
      row.push(augmented[i][n + j])
    }
    inverse.push(row)
  }
  Some(inverse)
}

///|
/// Determinant computed with pivoting, useful for diagnostics.
pub fn matrix_determinant(matrix : Array[Array[Double]]) -> Double {
  let n = matrix.length()
  if n == 0 || !matrix_is_square(matrix) {
    panic()
  }
  let working = matrix_copy(matrix)
  let mut determinant = 1.0
  let mut sign = 1.0
  for column in 0..<n {
    let mut pivot = column
    for row in (column + 1)..<n {
      if working[row][column].abs() > working[pivot][column].abs() {
        pivot = row
      }
    }
    if working[pivot][column].abs() <= survival_epsilon {
      return 0.0
    }
    if pivot != column {
      let temporary = working[column]
      working[column] = working[pivot]
      working[pivot] = temporary
      sign = -sign
    }
    let pivot_value = working[column][column]
    determinant *= pivot_value
    for row in (column + 1)..<n {
      let multiplier = working[row][column] / pivot_value
      for j in (column + 1)..<n {
        working[row][j] -= multiplier * working[column][j]
      }
    }
  }
  sign * determinant
}

///|
/// Symmetrize a numerical Hessian or covariance matrix.
pub fn matrix_symmetrize(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  if !matrix_is_square(matrix) {
    panic()
  }
  let result = matrix_copy(matrix)
  for i in 0..<result.length() {
    for j in (i + 1)..<result.length() {
      let value = 0.5 * (result[i][j] + result[j][i])
      result[i][j] = value
      result[j][i] = value
    }
  }
  result
}

///|
/// Add a small ridge to the diagonal to make near-singular model matrices
/// diagnosable and solvable without changing well-conditioned estimates.
pub fn matrix_ridge(
  matrix : Array[Array[Double]],
  ridge : Double,
) -> Array[Array[Double]] {
  if !matrix_is_square(matrix) {
    panic()
  }
  let result = matrix_copy(matrix)
  for i in 0..<result.length() {
    result[i][i] += ridge
  }
  result
}

///|
/// Return the diagonal of a square matrix.
pub fn matrix_diagonal(matrix : Array[Array[Double]]) -> Array[Double] {
  if !matrix_is_square(matrix) {
    panic()
  }
  let result = []
  for i in 0..<matrix.length() {
    result.push(matrix[i][i])
  }
  result
}

///|
/// Format a matrix as rows separated by newlines for diagnostic output.
pub fn matrix_to_string(matrix : Array[Array[Double]]) -> String {
  let builder = StringBuilder::new()
  for row in matrix {
    for j in 0..<row.length() {
      if j > 0 {
        builder.write_string(",")
      }
      builder.write_string("\{row[j]}")
    }
    builder.write_string("\n")
  }
  builder.to_string()
}
