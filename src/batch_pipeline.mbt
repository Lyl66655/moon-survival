///|
/// Stage status in an analysis pipeline.
pub(all) enum SurvivalPipelineStatus {
  Pending
  Running
  Succeeded
  Failed
} derive(Eq, Debug)

///|
/// A pipeline stage result.
pub(all) struct SurvivalPipelineStage {
  name : String
  status : SurvivalPipelineStatus
  rows_in : Int
  rows_out : Int
  events_out : Int
  message : String
  pipeline_stage : Bool
} derive(Eq, Debug)

///|
/// Batch-analysis pipeline configuration.
pub(all) struct SurvivalPipelineConfig {
  name : String
  horizon : Double
  run_quality_checks : Bool
  run_bootstrap : Bool
  bootstrap_replicates : Int
  pipeline_config : Bool
} derive(Eq, Debug)

///|
/// Completed pipeline metadata and stages.
pub(all) struct SurvivalPipelineRun {
  config : SurvivalPipelineConfig
  stages : Array[SurvivalPipelineStage]
  status : SurvivalPipelineStatus
  manifest : SurvivalAuditManifest
  pipeline_run : Bool
} derive(Eq, Debug)

///|
/// Create default pipeline configuration.
pub fn survival_pipeline_config(
  name : String,
  horizon : Double,
) -> SurvivalPipelineConfig {
  {
    name,
    horizon: if horizon < 0.0 {
      0.0
    } else {
      horizon
    },
    run_quality_checks: true,
    run_bootstrap: true,
    bootstrap_replicates: 200,
    pipeline_config: true,
  }
}

///|
/// Add a pipeline stage.
pub fn SurvivalPipelineRun::add_stage(
  self : SurvivalPipelineRun,
  stage : SurvivalPipelineStage,
) -> SurvivalPipelineRun {
  let stages : Array[SurvivalPipelineStage] = []
  for value in self.stages {
    stages.push(value)
  }
  stages.push(stage)
  let status = match stage.status {
    Failed => Failed
    _ => self.status
  }
  {
    config: self.config,
    stages,
    status,
    manifest: self.manifest,
    pipeline_run: true,
  }
}

///|
/// Start a pipeline run.
pub fn survival_pipeline_start(
  config : SurvivalPipelineConfig,
  data : SurvivalData,
) -> SurvivalPipelineRun {
  {
    config,
    stages: [],
    status: Running,
    manifest: survival_audit_manifest(
      "moon-survival",
      config.name,
      "0.2.1",
      data,
    ),
    pipeline_run: true,
  }
}

///|
/// Record validation and quality stages.
pub fn survival_pipeline_validate(
  run : SurvivalPipelineRun,
  data : SurvivalData,
) -> SurvivalPipelineRun {
  let issues = survival_quality_check(data, survival_default_quality_rules())
  let status = if issues.passed { Succeeded } else { Failed }
  run.add_stage({
    name: "quality",
    status,
    rows_in: data.length(),
    rows_out: data.length(),
    events_out: data.event_count(),
    message: if issues.passed {
      "quality checks passed"
    } else {
      "quality checks failed"
    },
    pipeline_stage: true,
  })
}

///|
/// Record a Kaplan-Meier estimation stage.
pub fn survival_pipeline_fit_km(
  run : SurvivalPipelineRun,
  data : SurvivalData,
) -> SurvivalPipelineRun {
  let fit = kaplan_meier(data)
  run.add_stage({
    name: "kaplan_meier",
    status: Succeeded,
    rows_in: data.length(),
    rows_out: fit.sample_size,
    events_out: fit.event_count,
    message: "KM curve estimated",
    pipeline_stage: true,
  })
}

///|
/// Record a bootstrap uncertainty stage.
pub fn survival_pipeline_bootstrap(
  run : SurvivalPipelineRun,
  data : SurvivalData,
  seed : UInt,
) -> SurvivalPipelineRun {
  if !run.config.run_bootstrap {
    return run
  }
  let estimate = survival_bootstrap_km(
    data,
    run.config.horizon,
    run.config.bootstrap_replicates,
    seed,
  )
  run.add_stage({
    name: "bootstrap",
    status: if estimate.replicates > 0 {
      Succeeded
    } else {
      Failed
    },
    rows_in: data.length(),
    rows_out: estimate.replicates,
    events_out: data.event_count(),
    message: "bootstrap uncertainty estimated",
    pipeline_stage: true,
  })
}

///|
/// Finish a run after all stages.
pub fn survival_pipeline_finish(
  run : SurvivalPipelineRun,
) -> SurvivalPipelineRun {
  let mut failed = false
  for stage in run.stages {
    if stage.status == Failed {
      failed = true
    }
  }
  {
    config: run.config,
    stages: run.stages,
    status: if failed {
      Failed
    } else {
      Succeeded
    },
    manifest: run.manifest,
    pipeline_run: true,
  }
}

///|
/// Return the number of successful stages.
pub fn SurvivalPipelineRun::successful_stages(
  self : SurvivalPipelineRun,
) -> Int {
  let mut count = 0
  for stage in self.stages {
    if stage.status == Succeeded {
      count += 1
    }
  }
  count
}

///|
/// Return a pipeline status label.
pub fn survival_pipeline_status_name(status : SurvivalPipelineStatus) -> String {
  match status {
    Pending => "pending"
    Running => "running"
    Succeeded => "succeeded"
    Failed => "failed"
  }
}

///|
/// Export pipeline stages.
pub fn survival_pipeline_to_csv(run : SurvivalPipelineRun) -> String {
  let builder = StringBuilder::new()
  builder.write_string("stage,status,rows_in,rows_out,events_out,message\n")
  for stage in run.stages {
    builder.write_string(stage.name)
    builder.write_char(',')
    builder.write_string(survival_pipeline_status_name(stage.status))
    builder.write_char(',')
    builder.write_string(stage.rows_in.to_string())
    builder.write_char(',')
    builder.write_string(stage.rows_out.to_string())
    builder.write_char(',')
    builder.write_string(stage.events_out.to_string())
    builder.write_char(',')
    builder.write_string(stage.message)
    builder.write_char('\n')
  }
  builder.to_string()
}

///|
/// Validate pipeline stage accounting.
pub fn survival_pipeline_is_valid(run : SurvivalPipelineRun) -> Bool {
  if run.config.name == "" || run.stages.length() == 0 {
    return false
  }
  for stage in run.stages {
    if stage.name == "" ||
      stage.rows_in < 0 ||
      stage.rows_out < 0 ||
      stage.events_out < 0 {
      return false
    }
  }
  true
}
