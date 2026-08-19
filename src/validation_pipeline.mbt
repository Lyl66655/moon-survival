///|
/// A deterministic fold assignment for reproducible model validation.
pub(all) struct ValidationFold {
  fold : Int
  train_indices : Array[Int]
  test_indices : Array[Int]
  validation_fold : Bool
} derive(Eq, Debug)

///|
/// Metrics recorded for one validation fold.
pub(all) struct ValidationFoldResult {
  fold : Int
  train_size : Int
  test_size : Int
  c_index : Double
  auc : Double
  brier_score : Double
  integrated_error : Double
  validation_result : Bool
} derive(Eq, Debug)

///|
/// Cross-validation summary across folds.
pub(all) struct SurvivalValidationReport {
  folds : Array[ValidationFoldResult]
  mean_c_index : Double
  standard_error_c_index : Double
  mean_auc : Double
  mean_brier_score : Double
  model_name : String
  validation_report : Bool
} derive(Eq, Debug)

///|
/// Bootstrap optimism-correction summary.
pub(all) struct BootstrapValidationResult {
  apparent_c_index : Double
  bootstrap_c_index : Double
  optimism : Double
  corrected_c_index : Double
  replicates : Int
  seed : Int
} derive(Eq, Debug)

///|
/// A prediction vector with its evaluation horizon.
pub(all) struct SurvivalPredictionSet {
  horizon : Double
  probabilities : Array[Double]
  prediction_set : Bool
} derive(Eq, Debug)

///|
/// Create stratified, deterministic folds by cycling through event strata.
pub fn validation_folds(
  data : SurvivalData,
  fold_count : Int,
) -> Array[ValidationFold] {
  if data.length() == 0 || fold_count <= 1 || fold_count > data.length() {
    return []
  }
  let bounded = if fold_count > data.length() {
    data.length()
  } else {
    fold_count
  }
  let assignments = Array::make(data.length(), 0)
  let mut event_cursor = 0
  let mut censor_cursor = 0
  for i in 0..<data.length() {
    if data.get(i).event == 1 {
      assignments[i] = event_cursor % bounded
      event_cursor += 1
    } else {
      assignments[i] = censor_cursor % bounded
      censor_cursor += 1
    }
  }
  let result : Array[ValidationFold] = []
  for fold in 0..<bounded {
    let train : Array[Int] = []
    let holdout_indices : Array[Int] = []
    for i in 0..<data.length() {
      if assignments[i] == fold {
        holdout_indices.push(i)
      } else {
        train.push(i)
      }
    }
    result.push({
      fold,
      train_indices: train,
      test_indices: holdout_indices,
      validation_fold: true,
    })
  }
  result
}

///|
/// Materialize a subset of SurvivalData by sorted row indices.
pub fn survival_subset(
  data : SurvivalData,
  indices : Array[Int],
) -> SurvivalData {
  let times : Array[Double] = []
  let events : Array[Int] = []
  for index in indices {
    if index >= 0 && index < data.length() {
      times.push(data.get(index).time)
      events.push(data.get(index).event)
    }
  }
  SurvivalData::new(times, events)
}

///|
/// Compute the mean of a finite metric vector.
fn validation_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  mean_doubles(values)
}

///|
/// Compute the standard error of a metric vector.
fn validation_standard_error(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let mean = validation_mean(values)
  let mut sum = 0.0
  for value in values {
    let delta = value - mean
    sum += delta * delta
  }
  double_sqrt(sum / (values.length() * (values.length() - 1)).to_double())
}

///|
/// Evaluate binary horizon predictions with C-index, AUC and Brier score.
pub fn validate_prediction_set(
  data : SurvivalData,
  prediction_set : SurvivalPredictionSet,
) -> ValidationFoldResult {
  let predictions = prediction_set.probabilities
  let risk_scores : Array[Double] = []
  for probability in predictions {
    risk_scores.push(1.0 - clamp_probability(probability))
  }
  let concordance = concordance_index(data, risk_scores)
  let auc = time_dependent_auc(data, risk_scores, prediction_set.horizon)
  let brier = ipcw_brier_score(data, predictions, prediction_set.horizon)
  {
    fold: 0,
    train_size: data.length(),
    test_size: data.length(),
    c_index: concordance.c_index,
    auc,
    brier_score: brier,
    integrated_error: brier,
    validation_result: true,
  }
}

///|
/// Fit a Cox model on each training split and score it on the held-out split.
pub fn cox_cross_validate(
  data : CoxData,
  survival_data : SurvivalData,
  fold_count : Int,
  model_name? : String = "cox",
) -> SurvivalValidationReport {
  if data.length() != survival_data.length() {
    return {
      folds: [],
      mean_c_index: 0.0,
      standard_error_c_index: 0.0,
      mean_auc: 0.0,
      mean_brier_score: 0.0,
      model_name,
      validation_report: false,
    }
  }
  let folds = validation_folds(survival_data, fold_count)
  let results : Array[ValidationFoldResult] = []
  for fold in folds {
    let train_times : Array[Double] = []
    let train_events : Array[Int] = []
    let train_covs : Array[Array[Double]] = []
    for index in fold.train_indices {
      let observation = data.get(index)
      train_times.push(observation.time)
      train_events.push(observation.event)
      train_covs.push(observation.covariates)
    }
    let fit = cox_fit(CoxData::new(train_times, train_events, train_covs))
    let test_data = survival_subset(survival_data, fold.test_indices)
    let test_covs : Array[Array[Double]] = []
    for index in fold.test_indices {
      test_covs.push(data.get(index).covariates)
    }
    let risks : Array[Double] = []
    let predictions : Array[Double] = []
    for covariates in test_covs {
      risks.push(fit.risk(covariates))
      predictions.push(
        fit.survival_at(
          covariates,
          if test_data.length() == 0 {
            0.0
          } else {
            test_data.get(test_data.length() - 1).time
          },
        ),
      )
    }
    let concordance = concordance_index(test_data, risks)
    let horizon = if test_data.length() == 0 {
      0.0
    } else {
      test_data.get(test_data.length() - 1).time
    }
    let auc = time_dependent_auc(test_data, risks, horizon)
    let brier = ipcw_brier_score(test_data, predictions, horizon)
    results.push({
      fold: fold.fold,
      train_size: fold.train_indices.length(),
      test_size: fold.test_indices.length(),
      c_index: concordance.c_index,
      auc,
      brier_score: brier,
      integrated_error: brier,
      validation_result: true,
    })
  }
  let c_indices : Array[Double] = []
  let aucs : Array[Double] = []
  let briers : Array[Double] = []
  for result in results {
    c_indices.push(result.c_index)
    aucs.push(result.auc)
    briers.push(result.brier_score)
  }
  {
    folds: results,
    mean_c_index: validation_mean(c_indices),
    standard_error_c_index: validation_standard_error(c_indices),
    mean_auc: validation_mean(aucs),
    mean_brier_score: validation_mean(briers),
    model_name,
    validation_report: results.length() > 0,
  }
}

///|
/// Fit a Kaplan-Meier training curve and validate its marginal prediction.
pub fn km_cross_validate(
  data : SurvivalData,
  fold_count : Int,
  horizon : Double,
) -> SurvivalValidationReport {
  let folds = validation_folds(data, fold_count)
  let results : Array[ValidationFoldResult] = []
  for fold in folds {
    let train = survival_subset(data, fold.train_indices)
    let holdout = survival_subset(data, fold.test_indices)
    let fit = kaplan_meier(train)
    let prediction = fit.survival_at(horizon)
    let predictions : Array[Double] = Array::make(holdout.length(), prediction)
    let risk : Array[Double] = Array::make(holdout.length(), 1.0 - prediction)
    let concordance = concordance_index(holdout, risk)
    let auc = time_dependent_auc(holdout, risk, horizon)
    let brier = ipcw_brier_score(holdout, predictions, horizon)
    results.push({
      fold: fold.fold,
      train_size: train.length(),
      test_size: holdout.length(),
      c_index: concordance.c_index,
      auc,
      brier_score: brier,
      integrated_error: brier,
      validation_result: true,
    })
  }
  let c_indices : Array[Double] = []
  let aucs : Array[Double] = []
  let briers : Array[Double] = []
  for result in results {
    c_indices.push(result.c_index)
    aucs.push(result.auc)
    briers.push(result.brier_score)
  }
  {
    folds: results,
    mean_c_index: validation_mean(c_indices),
    standard_error_c_index: validation_standard_error(c_indices),
    mean_auc: validation_mean(aucs),
    mean_brier_score: validation_mean(briers),
    model_name: "kaplan-meier",
    validation_report: results.length() > 0,
  }
}

///|
/// Deterministic pseudo-random bootstrap index generation.
fn validation_next_seed(seed : Int) -> Int {
  let value = (seed * 1103515245 + 12345) % 2147483647
  if value <= 0 {
    value + 2147483646
  } else {
    value
  }
}

///|
/// Draw a bootstrap sample with replacement.
pub fn validation_bootstrap_indices(
  size : Int,
  seed : Int,
) -> (Array[Int], Int) {
  if size <= 0 {
    return ([], seed)
  }
  let result : Array[Int] = []
  let mut current = if seed <= 0 { 1 } else { seed }
  for _ in 0..<size {
    current = validation_next_seed(current)
    result.push(current % size)
  }
  (result, current)
}

///|
/// Estimate optimism by bootstrap resampling a risk-score validation set.
pub fn bootstrap_optimism_correction(
  data : SurvivalData,
  risk_scores : Array[Double],
  replicates : Int,
  seed : Int,
) -> BootstrapValidationResult {
  if risk_scores.length() != data.length() ||
    data.length() == 0 ||
    replicates <= 0 {
    return {
      apparent_c_index: 0.5,
      bootstrap_c_index: 0.5,
      optimism: 0.0,
      corrected_c_index: 0.5,
      replicates: 0,
      seed,
    }
  }
  let apparent = concordance_index(data, risk_scores).c_index
  let mut current_seed = seed
  let mut optimism_total = 0.0
  for _ in 0..<replicates {
    let (indices, next_seed) = validation_bootstrap_indices(
      data.length(),
      current_seed,
    )
    current_seed = next_seed
    let bootstrap_data = survival_subset(data, indices)
    let bootstrap_scores : Array[Double] = []
    for index in indices {
      bootstrap_scores.push(risk_scores[index])
    }
    let bootstrap_fit = concordance_index(bootstrap_data, bootstrap_scores).c_index
    let original_fit = concordance_index(data, bootstrap_scores).c_index
    optimism_total += bootstrap_fit - original_fit
  }
  let optimism = optimism_total / replicates.to_double()
  let corrected = apparent - optimism
  {
    apparent_c_index: apparent,
    bootstrap_c_index: apparent + optimism,
    optimism,
    corrected_c_index: corrected,
    replicates,
    seed,
  }
}

///|
/// Compute a calibration table for a prediction set.
pub fn validation_calibration(
  data : SurvivalData,
  prediction_set : SurvivalPredictionSet,
  bins : Int,
) -> Array[IpcwCalibrationBin] {
  ipcw_calibration(
    data,
    prediction_set.probabilities,
    prediction_set.horizon,
    bins~,
  )
}

///|
/// Export a cross-validation report.
pub fn SurvivalValidationReport::to_csv(
  self : SurvivalValidationReport,
) -> String {
  let builder = StringBuilder::new()
  builder.write_string("fold,train_size,test_size,c_index,auc,brier_score\n")
  for result in self.folds {
    builder.write_string(result.fold.to_string())
    builder.write_char(',')
    builder.write_string(result.train_size.to_string())
    builder.write_char(',')
    builder.write_string(result.test_size.to_string())
    builder.write_char(',')
    builder.write_string(result.c_index.to_string())
    builder.write_char(',')
    builder.write_string(result.auc.to_string())
    builder.write_char(',')
    builder.write_string(result.brier_score.to_string())
    builder.write_char('\n')
  }
  builder.write_string("mean_c_index,")
  builder.write_string(self.mean_c_index.to_string())
  builder.write_char('\n')
  builder.write_string("standard_error_c_index,")
  builder.write_string(self.standard_error_c_index.to_string())
  builder.write_char('\n')
  builder.write_string("mean_auc,")
  builder.write_string(self.mean_auc.to_string())
  builder.write_char('\n')
  builder.write_string("mean_brier_score,")
  builder.write_string(self.mean_brier_score.to_string())
  builder.write_char('\n')
  builder.to_string()
}

///|
/// Check a validation report's finite metric invariants.
pub fn validation_report_is_valid(report : SurvivalValidationReport) -> Bool {
  if report.folds.length() == 0 ||
    report.mean_c_index.is_nan() ||
    report.mean_auc.is_nan() ||
    report.mean_brier_score.is_nan() {
    return false
  }
  for fold in report.folds {
    if fold.train_size <= 0 ||
      fold.test_size <= 0 ||
      fold.c_index < 0.0 ||
      fold.c_index > 1.0 ||
      fold.auc < 0.0 ||
      fold.auc > 1.0 ||
      fold.brier_score < 0.0 {
      return false
    }
  }
  true
}
