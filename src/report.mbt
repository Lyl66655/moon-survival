///|
/// Human-readable report assembled from the library's core estimators.
pub(all) struct SurvivalReport {
  title : String
  summary : SurvivalSummary
  median : Double?
  restricted_mean : Double
  exponential_rate : Double
  weibull_shape : Double
  weibull_scale : Double
  logrank : LogRankResult?
  quality_flags : Array[String]
} derive(Eq, Debug)

///|
/// Build a report for one dataset. The report is deliberately data-only so it
/// can be rendered by a CLI, notebook, or web application.
pub fn build_survival_report(
  title : String,
  data : SurvivalData,
  horizon : Double,
) -> SurvivalReport {
  let km = kaplan_meier(data)
  let exponential = exponential_fit(data)
  let weibull = weibull_fit(data)
  let flags : Array[String] = []
  if data.length() == 0 {
    flags.push("empty_dataset")
  }
  if data.length() > 0 && data.event_count() == 0 {
    flags.push("no_observed_events")
  }
  if data.length() > 0 &&
    data.censored_count().to_double() / data.length().to_double() > 0.5 {
    flags.push("high_censoring")
  }
  {
    title,
    summary: data.summary(),
    median: km.median(),
    restricted_mean: km.restricted_mean(horizon),
    exponential_rate: exponential.rate,
    weibull_shape: weibull.shape,
    weibull_scale: weibull.scale,
    logrank: None,
    quality_flags: flags,
  }
}

///|
/// Add a two-group comparison to a report.
pub fn SurvivalReport::with_logrank(
  self : SurvivalReport,
  result : LogRankResult,
) -> SurvivalReport {
  {
    title: self.title,
    summary: self.summary,
    median: self.median,
    restricted_mean: self.restricted_mean,
    exponential_rate: self.exponential_rate,
    weibull_shape: self.weibull_shape,
    weibull_scale: self.weibull_scale,
    logrank: Some(result),
    quality_flags: self.quality_flags,
  }
}

///|
/// Render a concise Markdown report suitable for a README artifact.
pub fn SurvivalReport::to_markdown(self : SurvivalReport) -> String {
  let builder = StringBuilder::new()
  builder.write_string("# ")
  builder.write_string(self.title)
  builder.write_string("\n\n")
  builder.write_string("| Metric | Value |\n|---|---:|\n")
  builder.write_string("| Sample size | \{self.summary.n} |\n")
  builder.write_string("| Events | \{self.summary.events} |\n")
  builder.write_string(
    "| Censoring fraction | \{1.0 - self.summary.event_fraction} |\n",
  )
  builder.write_string("| Restricted mean | \{self.restricted_mean} |\n")
  builder.write_string("| Exponential rate | \{self.exponential_rate} |\n")
  builder.write_string("| Weibull shape | \{self.weibull_shape} |\n")
  builder.write_string("| Weibull scale | \{self.weibull_scale} |\n")
  match self.median {
    Some(value) => builder.write_string("| Median | \{value} |\n")
    None => builder.write_string("| Median | not reached |\n")
  }
  match self.logrank {
    Some(value) =>
      builder.write_string("| Log-rank p-value | \{value.p_value} |\n")
    None => ()
  }
  if self.quality_flags.length() > 0 {
    builder.write_string("\nQuality flags: ")
    for i in 0..<self.quality_flags.length() {
      if i > 0 {
        builder.write_string(", ")
      }
      builder.write_string(self.quality_flags[i])
    }
    builder.write_string("\n")
  }
  builder.to_string()
}

///|
/// Render a machine-readable, stable CSV summary.
pub fn SurvivalReport::to_csv(self : SurvivalReport) -> String {
  let builder = StringBuilder::new()
  builder.write_string("metric,value\n")
  builder.write_string("sample_size,\{self.summary.n}\n")
  builder.write_string("events,\{self.summary.events}\n")
  builder.write_string("censored,\{self.summary.censored}\n")
  builder.write_string("event_fraction,\{self.summary.event_fraction}\n")
  builder.write_string("restricted_mean,\{self.restricted_mean}\n")
  builder.write_string("exponential_rate,\{self.exponential_rate}\n")
  builder.write_string("weibull_shape,\{self.weibull_shape}\n")
  builder.write_string("weibull_scale,\{self.weibull_scale}\n")
  builder.to_string()
}

///|
/// Return true when the report contains no high-risk quality flag.
pub fn SurvivalReport::is_analysis_ready(self : SurvivalReport) -> Bool {
  self.summary.n > 0 && self.summary.events > 0
}
