///|
/// A structured result for public API input validation.
pub(all) struct SurvivalApiGuardResult {
  valid : Bool
  code : String
  message : String
  guard_row : Bool
} derive(Eq, Debug)

///|
/// Return a successful guard result.
pub fn survival_guard_ok() -> SurvivalApiGuardResult {
  { valid: true, code: "OK", message: "input accepted", guard_row: true }
}

///|
/// Return a failed guard result.
pub fn survival_guard_error(
  code : String,
  message : String,
) -> SurvivalApiGuardResult {
  { valid: false, code, message, guard_row: true }
}

///|
/// Validate parallel array lengths.
pub fn survival_guard_lengths(
  first : Int,
  second : Int,
  label : String,
) -> SurvivalApiGuardResult {
  if first == second {
    survival_guard_ok()
  } else {
    survival_guard_error("LENGTH_MISMATCH", label + " lengths differ")
  }
}

///|
/// Validate a non-empty sample size.
pub fn survival_guard_nonempty(size : Int) -> SurvivalApiGuardResult {
  if size > 0 {
    survival_guard_ok()
  } else {
    survival_guard_error("EMPTY_DATA", "data must contain at least one row")
  }
}

///|
/// Validate a probability.
pub fn survival_guard_probability(
  value : Double,
  label : String,
) -> SurvivalApiGuardResult {
  if value.is_nan() || value.is_inf() || value < 0.0 || value > 1.0 {
    survival_guard_error("INVALID_PROBABILITY", label + " must be in [0,1]")
  } else {
    survival_guard_ok()
  }
}

///|
/// Validate a positive parameter.
pub fn survival_guard_positive(
  value : Double,
  label : String,
) -> SurvivalApiGuardResult {
  if value.is_nan() || value.is_inf() || value <= 0.0 {
    survival_guard_error("INVALID_POSITIVE", label + " must be positive")
  } else {
    survival_guard_ok()
  }
}

///|
/// Validate a non-negative duration.
pub fn survival_guard_time(value : Double) -> SurvivalApiGuardResult {
  if value.is_nan() || value.is_inf() || value < 0.0 {
    survival_guard_error("INVALID_TIME", "time must be finite and non-negative")
  } else {
    survival_guard_ok()
  }
}

///|
/// Validate a binary event code.
pub fn survival_guard_event(value : Int) -> SurvivalApiGuardResult {
  if value == 0 || value == 1 {
    survival_guard_ok()
  } else {
    survival_guard_error("INVALID_EVENT", "event must be 0 or 1")
  }
}

///|
/// Validate a tolerance.
pub fn survival_guard_tolerance(value : Double) -> SurvivalApiGuardResult {
  if value.is_nan() || value.is_inf() || value <= 0.0 {
    survival_guard_error("INVALID_TOLERANCE", "tolerance must be positive")
  } else {
    survival_guard_ok()
  }
}

///|
/// Validate a finite vector.
pub fn survival_guard_vector(
  values : Array[Double],
  label : String,
) -> SurvivalApiGuardResult {
  for value in values {
    if value.is_nan() || value.is_inf() {
      return survival_guard_error(
        "NONFINITE_VALUE",
        label + " contains a non-finite value",
      )
    }
  }
  survival_guard_ok()
}

///|
/// Validate a matrix shape.
pub fn survival_guard_matrix(
  matrix : Array[Array[Double]],
  rows : Int,
  columns : Int,
) -> SurvivalApiGuardResult {
  if matrix.length() != rows {
    return survival_guard_error("ROW_MISMATCH", "matrix row count differs")
  }
  for row in matrix {
    if row.length() != columns {
      return survival_guard_error(
        "COLUMN_MISMATCH", "matrix column count differs",
      )
    }
  }
  survival_guard_ok()
}

///|
/// Validate a confidence level.
pub fn survival_guard_confidence(value : Double) -> SurvivalApiGuardResult {
  if value > 0.0 && value < 1.0 {
    survival_guard_ok()
  } else {
    survival_guard_error(
      "INVALID_CONFIDENCE", "confidence must be between zero and one",
    )
  }
}

///|
/// Validate an iteration budget.
pub fn survival_guard_iterations(value : Int) -> SurvivalApiGuardResult {
  if value > 0 {
    survival_guard_ok()
  } else {
    survival_guard_error(
      "INVALID_ITERATIONS", "iteration budget must be positive",
    )
  }
}

///|
/// Return the first failed guard.
pub fn survival_guard_first(
  results : Array[SurvivalApiGuardResult],
) -> SurvivalApiGuardResult {
  for result in results {
    if !result.valid {
      return result
    }
  }
  survival_guard_ok()
}

///|
/// Export a guard result.
pub fn survival_guard_to_csv(result : SurvivalApiGuardResult) -> String {
  "valid,code,message\n" +
  result.valid.to_string() +
  "," +
  result.code +
  "," +
  result.message +
  "\n"
}

///|
/// Validate survival parallel arrays.
pub fn survival_guard_survival_arrays(
  times : Array[Double],
  events : Array[Int],
) -> SurvivalApiGuardResult {
  let length = survival_guard_lengths(
    times.length(),
    events.length(),
    "time/event",
  )
  if !length.valid {
    return length
  }
  let nonempty = survival_guard_nonempty(times.length())
  if !nonempty.valid {
    return nonempty
  }
  for time in times {
    let check_result = survival_guard_time(time)
    if !check_result.valid {
      return check_result
    }
  }
  for event in events {
    let check_result = survival_guard_event(event)
    if !check_result.valid {
      return check_result
    }
  }
  survival_guard_ok()
}

///|
/// Validate a target horizon relative to an observed maximum.
pub fn survival_guard_horizon(
  horizon : Double,
  maximum_time : Double,
) -> SurvivalApiGuardResult {
  let time = survival_guard_time(horizon)
  if !time.valid {
    time
  } else if horizon > maximum_time {
    survival_guard_error(
      "HORIZON_AFTER_DATA", "horizon exceeds observed follow-up",
    )
  } else {
    survival_guard_ok()
  }
}
