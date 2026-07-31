///|
pub fn double_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut z = x
  for _ in 0..<20 {
    let next_z = z - (z * z - x) / (2.0 * z)
    let diff = next_z - z
    z = next_z
    if diff < 0.000000001 && diff > -0.000000001 {
      break
    }
  }
  z
}

///|
pub fn double_exp(x : Double) -> Double {
  let mut sum = 1.0
  let mut term = 1.0
  for i in 1..=30 {
    term = term * x / i.to_double()
    sum = sum + term
    if term < 0.000000001 && term > -0.000000001 {
      break
    }
  }
  sum
}
