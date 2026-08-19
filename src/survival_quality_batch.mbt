///|
/// Batch quality gates for multiple survival cohorts.
pub(all) struct SurvivalQualityBatchRow {
  batch_row : Int
  dataset_name : String
  sample_size : Int
  events : Int
  censored : Int
  issue_count : Int
  passed : Bool
  completeness : Double
  quality_row : Int
} derive(Eq, Debug)

///|
pub(all) struct SurvivalQualityGate {
  minimum_sample_size : Int
  minimum_event_count : Int
  maximum_censoring_fraction : Double
  require_nonnegative_time : Bool
  gate_row : Int
} derive(Eq, Debug)

///|
fn survival_quality_issue_count(report : SurvivalQualityReport) -> Int {
  report.issues.length()
}

///|
/// Evaluate named datasets against a common quality gate.
pub fn survival_quality_batch(
  names : Array[String],
  datasets : Array[SurvivalData],
  gate : SurvivalQualityGate,
) -> Array[SurvivalQualityBatchRow] {
  let output : Array[SurvivalQualityBatchRow] = []
  let count = names.length().min(datasets.length())
  for i = 0; i < count; i = i + 1 {
    let data = datasets[i]
    let report = survival_quality_check(data, survival_default_quality_rules())
    let sample = data.length()
    let events = data.event_count()
    let censored = data.censored_count()
    let fraction = if sample == 0 {
      1.0
    } else {
      censored.to_double() / sample.to_double()
    }
    let passed = sample >= gate.minimum_sample_size &&
      events >= gate.minimum_event_count &&
      fraction <= gate.maximum_censoring_fraction &&
      survival_quality_issue_count(report) == 0
    output.push({
      batch_row: i,
      dataset_name: names[i],
      sample_size: sample,
      events,
      censored,
      issue_count: survival_quality_issue_count(report),
      passed,
      completeness: if sample == 0 {
        0.0
      } else {
        1.0
      },
      quality_row: i,
    })
  }
  output
}

///|
/// Create a conservative gate suitable for exploratory analyses.
pub fn survival_quality_gate_default() -> SurvivalQualityGate {
  {
    minimum_sample_size: 1,
    minimum_event_count: 1,
    maximum_censoring_fraction: 0.99,
    require_nonnegative_time: true,
    gate_row: 0,
  }
}

///|
/// Return a stricter gate for production reporting.
pub fn survival_quality_gate_production(
  minimum_sample_size : Int,
  minimum_event_count : Int,
  maximum_censoring_fraction : Double,
) -> SurvivalQualityGate {
  {
    minimum_sample_size: minimum_sample_size.max(1),
    minimum_event_count: minimum_event_count.max(1),
    maximum_censoring_fraction: clamp_probability(maximum_censoring_fraction),
    require_nonnegative_time: true,
    gate_row: 1,
  }
}

///|
/// Number of datasets that pass the gate.
pub fn survival_quality_batch_pass_count(
  rows : Array[SurvivalQualityBatchRow],
) -> Int {
  let mut count = 0
  for row in rows {
    if row.passed {
      count += 1
    }
  }
  count
}

///|
pub fn survival_quality_batch_to_csv(
  rows : Array[SurvivalQualityBatchRow],
) -> String {
  let mut out = "batch_row,dataset_name,sample_size,events,censored,issue_count,passed,completeness,quality_row\n"
  for row in rows {
    out = out +
      row.batch_row.to_string() +
      "," +
      row.dataset_name +
      "," +
      row.sample_size.to_string() +
      "," +
      row.events.to_string() +
      "," +
      row.censored.to_string() +
      "," +
      row.issue_count.to_string() +
      "," +
      row.passed.to_string() +
      "," +
      row.completeness.to_string() +
      "," +
      row.quality_row.to_string() +
      "\n"
  }
  out
}

///|
pub fn survival_quality_gate_to_csv(gate : SurvivalQualityGate) -> String {
  "minimum_sample_size,minimum_event_count,maximum_censoring_fraction,require_nonnegative_time,gate_row\n" +
  gate.minimum_sample_size.to_string() +
  "," +
  gate.minimum_event_count.to_string() +
  "," +
  gate.maximum_censoring_fraction.to_string() +
  "," +
  gate.require_nonnegative_time.to_string() +
  "," +
  gate.gate_row.to_string() +
  "\n"
}

///|
pub fn survival_quality_batch_is_valid(
  rows : Array[SurvivalQualityBatchRow],
) -> Bool {
  for row in rows {
    if row.sample_size < 0 ||
      row.events < 0 ||
      row.censored < 0 ||
      row.events + row.censored != row.sample_size ||
      row.issue_count < 0 ||
      row.completeness < 0.0 ||
      row.completeness > 1.0 {
      return false
    }
  }
  true
}

///|
pub fn survival_quality_gate_is_valid(gate : SurvivalQualityGate) -> Bool {
  gate.minimum_sample_size >= 0 &&
  gate.minimum_event_count >= 0 &&
  gate.maximum_censoring_fraction >= 0.0 &&
  gate.maximum_censoring_fraction <= 1.0
}
