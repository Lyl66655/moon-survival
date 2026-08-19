///|
/// Severity for a data-quality issue.
pub(all) enum SurvivalQualitySeverity {
  Info
  Warning
  Error
} derive(Eq, Debug)

///|
/// A machine-readable data-quality issue.
pub(all) struct SurvivalQualityIssue {
  code : String
  severity : SurvivalQualitySeverity
  row : Int
  column : Int
  message : String
  quality_issue : Bool
} derive(Eq, Debug)

///|
/// A validation rule applied to a survival table.
pub(all) struct SurvivalQualityRule {
  code : String
  description : String
  severity : SurvivalQualitySeverity
  enabled : Bool
} derive(Eq, Debug)

///|
/// Full quality-control report.
pub(all) struct SurvivalQualityReport {
  row_count : Int
  issue_count : Int
  error_count : Int
  warning_count : Int
  issues : Array[SurvivalQualityIssue]
  passed : Bool
} derive(Eq, Debug)

///|
/// Default survival-data quality rules.
pub fn survival_default_quality_rules() -> Array[SurvivalQualityRule] {
  [
    {
      code: "TIME_FINITE",
      description: "time must be finite and non-negative",
      severity: Error,
      enabled: true,
    },
    {
      code: "EVENT_BINARY",
      description: "event must be zero or one",
      severity: Error,
      enabled: true,
    },
    {
      code: "DUPLICATE_ROW",
      description: "duplicate subject rows should be reviewed",
      severity: Warning,
      enabled: true,
    },
    {
      code: "LATE_FOLLOWUP",
      description: "extreme follow-up values should be reviewed",
      severity: Info,
      enabled: true,
    },
  ]
}

///|
/// Check core survival data against basic rules.
pub fn survival_quality_check(
  data : SurvivalData,
  rules : Array[SurvivalQualityRule],
) -> SurvivalQualityReport {
  let has_rules = rules.length() > 0
  let issues : Array[SurvivalQualityIssue] = []
  for i in 0..<data.length() {
    let observation = data.get(i)
    if observation.time.is_nan() ||
      observation.time.is_inf() ||
      observation.time < 0.0 {
      issues.push({
        code: "TIME_FINITE",
        severity: Error,
        row: i,
        column: 0,
        message: "invalid time",
        quality_issue: true,
      })
    }
    if observation.event != 0 && observation.event != 1 {
      issues.push({
        code: "EVENT_BINARY",
        severity: Error,
        row: i,
        column: 1,
        message: "invalid event indicator",
        quality_issue: true,
      })
    }
    if i > 0 &&
      observation.time == data.get(i - 1).time &&
      observation.event == data.get(i - 1).event {
      issues.push({
        code: "DUPLICATE_ROW",
        severity: Warning,
        row: i,
        column: -1,
        message: "duplicate time/event row",
        quality_issue: true,
      })
    }
  }
  let mut errors = 0
  let mut warnings = 0
  for issue in issues {
    match issue.severity {
      Error => errors += 1
      Warning => warnings += 1
      Info => ()
    }
  }
  {
    row_count: data.length(),
    issue_count: issues.length(),
    error_count: errors,
    warning_count: warnings,
    issues,
    passed: if has_rules {
      errors == 0
    } else {
      errors == 0
    },
  }
}

///|
/// Detect outliers using a robust interquartile fence.
pub fn survival_time_outliers(
  data : SurvivalData,
  multiplier : Double,
) -> Array[Int] {
  let values : Array[Double] = []
  for observation in data.to_array() {
    values.push(observation.time)
  }
  if values.length() < 4 {
    return []
  }
  values.sort()
  let q1 = values[values.length() / 4]
  let q3 = values[3 * values.length() / 4]
  let fence = if multiplier <= 0.0 { 1.5 } else { multiplier }
  let lower = q1 - fence * (q3 - q1)
  let upper = q3 + fence * (q3 - q1)
  let result : Array[Int] = []
  for i in 0..<data.length() {
    if data.get(i).time < lower || data.get(i).time > upper {
      result.push(i)
    }
  }
  result
}

///|
/// Check whether event times respect a fixed administrative horizon.
pub fn survival_horizon_issues(
  data : SurvivalData,
  horizon : Double,
) -> Array[SurvivalQualityIssue] {
  let result : Array[SurvivalQualityIssue] = []
  for i in 0..<data.length() {
    if horizon > 0.0 && data.get(i).time > horizon {
      result.push({
        code: "LATE_FOLLOWUP",
        severity: Info,
        row: i,
        column: 0,
        message: "time exceeds administrative horizon",
        quality_issue: true,
      })
    }
  }
  result
}

///|
/// Check that a treatment/group vector is aligned and binary.
pub fn survival_group_quality(
  groups : Array[Int],
  expected_length : Int,
) -> Array[SurvivalQualityIssue] {
  let result : Array[SurvivalQualityIssue] = []
  if groups.length() != expected_length {
    result.push({
      code: "GROUP_LENGTH",
      severity: Error,
      row: -1,
      column: -1,
      message: "group vector length mismatch",
      quality_issue: true,
    })
  }
  for i in 0..<groups.length() {
    if groups[i] < 0 {
      result.push({
        code: "GROUP_VALUE",
        severity: Error,
        row: i,
        column: 0,
        message: "group must be non-negative",
        quality_issue: true,
      })
    }
  }
  result
}

///|
/// Combine quality issues from independent checks.
pub fn survival_quality_merge(
  row_count : Int,
  groups : Array[Array[SurvivalQualityIssue]],
) -> SurvivalQualityReport {
  let issues : Array[SurvivalQualityIssue] = []
  for issue_group in groups {
    for issue in issue_group {
      issues.push(issue)
    }
  }
  let mut errors = 0
  let mut warnings = 0
  for issue in issues {
    match issue.severity {
      Error => errors += 1
      Warning => warnings += 1
      Info => ()
    }
  }
  {
    row_count,
    issue_count: issues.length(),
    error_count: errors,
    warning_count: warnings,
    issues,
    passed: errors == 0,
  }
}

///|
/// Compute a completeness fraction for a numeric matrix.
pub fn survival_quality_completeness(values : Array[Array[Double]]) -> Double {
  let mut total = 0
  let mut observed = 0
  for row in values {
    for value in row {
      total += 1
      if !value.is_nan() && !value.is_inf() {
        observed += 1
      }
    }
  }
  if total == 0 {
    1.0
  } else {
    observed.to_double() / total.to_double()
  }
}

///|
/// Detect impossible event/censoring combinations.
pub fn survival_event_consistency(
  data : SurvivalData,
) -> Array[SurvivalQualityIssue] {
  let result : Array[SurvivalQualityIssue] = []
  for i in 0..<data.length() {
    if data.get(i).time == 0.0 && data.get(i).event == 0 {
      result.push({
        code: "ZERO_CENSORED_TIME",
        severity: Warning,
        row: i,
        column: -1,
        message: "zero-time censoring",
        quality_issue: true,
      })
    }
  }
  result
}

///|
/// Return a human-readable severity label.
pub fn survival_quality_severity_name(
  severity : SurvivalQualitySeverity,
) -> String {
  match severity {
    Info => "info"
    Warning => "warning"
    Error => "error"
  }
}

///|
/// Export quality issues.
pub fn survival_quality_issues_to_csv(
  issues : Array[SurvivalQualityIssue],
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("code,severity,row,column,message\n")
  for issue in issues {
    builder.write_string(issue.code)
    builder.write_char(',')
    builder.write_string(survival_quality_severity_name(issue.severity))
    builder.write_char(',')
    builder.write_string(issue.row.to_string())
    builder.write_char(',')
    builder.write_string(issue.column.to_string())
    builder.write_char(',')
    builder.write_string(issue.message)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Export quality summary.
pub fn survival_quality_report_to_csv(report : SurvivalQualityReport) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("row_count,")
  builder.write_string(report.row_count.to_string())
  builder.write_char('\n')
  builder.write_string("issue_count,")
  builder.write_string(report.issue_count.to_string())
  builder.write_char('\n')
  builder.write_string("error_count,")
  builder.write_string(report.error_count.to_string())
  builder.write_char('\n')
  builder.write_string("warning_count,")
  builder.write_string(report.warning_count.to_string())
  builder.write_char('\n')
  builder.write_string("passed,")
  builder.write_string(report.passed.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Verify quality-report invariants.
pub fn survival_quality_report_is_valid(report : SurvivalQualityReport) -> Bool {
  report.row_count >= 0 &&
  report.issue_count == report.issues.length() &&
  report.error_count >= 0 &&
  report.warning_count >= 0 &&
  report.error_count <= report.issue_count
}
