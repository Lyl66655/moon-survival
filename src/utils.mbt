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
