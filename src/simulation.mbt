///|
/// Configuration for reproducible survival-data simulations.
pub(all) struct SimulationConfig {
  sample_size : Int
  baseline_rate : Double
  censor_rate : Double
  shape : Double
  seed : UInt
} derive(Eq, Debug)

///|
/// A simulated right-censored dataset with the latent failure times retained
/// for validation experiments.
pub(all) struct SimulatedSurvival {
  data : SurvivalData
  latent_failure_times : Array[Double]
  censoring_times : Array[Double]
} derive(Eq, Debug)

///|
/// Create a sensible default simulation configuration.
pub fn default_simulation_config(sample_size : Int) -> SimulationConfig {
  {
    sample_size,
    baseline_rate: 0.08,
    censor_rate: 0.025,
    shape: 1.0,
    seed: 20260818U,
  }
}

///|
/// Simulate exponential or Weibull failure times with independent censoring.
pub fn simulate_survival(config : SimulationConfig) -> SimulatedSurvival {
  let rng = DeterministicRng::new(config.seed)
  let times = []
  let events = []
  let latent = []
  let censoring = []
  for _ in 0..<config.sample_size {
    let failure_uniform = clamp_positive(rng.uniform(), floor=0.0000000001)
    let failure = if config.shape <= 0.0 {
      -double_log(failure_uniform) / config.baseline_rate
    } else {
      @math.pow(
        -double_log(failure_uniform) / config.baseline_rate,
        1.0 / config.shape,
      )
    }
    let censor_uniform = clamp_positive(rng.uniform(), floor=0.0000000001)
    let censor = if config.censor_rate <= 0.0 {
      failure + 1000000.0
    } else {
      -double_log(censor_uniform) / config.censor_rate
    }
    let observed = if failure <= censor { failure } else { censor }
    let event = if failure <= censor { 1 } else { 0 }
    times.push(observed)
    events.push(event)
    latent.push(failure)
    censoring.push(censor)
  }
  {
    data: SurvivalData::new(times, events),
    latent_failure_times: latent,
    censoring_times: censoring,
  }
}

///|
/// Simulated Cox data with known coefficients for recovery tests.
pub(all) struct SimulatedCox {
  data : CoxData
  true_coefficients : Array[Double]
  latent_failure_times : Array[Double]
} derive(Eq, Debug)

///|
/// Generate a proportional-hazards sample with normally distributed
/// covariates and independent exponential censoring.
pub fn simulate_cox(
  sample_size : Int,
  true_coefficients : Array[Double],
  baseline_rate : Double,
  censor_rate : Double,
  seed : UInt,
) -> SimulatedCox {
  let rng = DeterministicRng::new(seed)
  let times = []
  let events = []
  let covariates = []
  let latent = []
  for _ in 0..<sample_size {
    let row = []
    for _ in true_coefficients {
      row.push(rng.normal())
    }
    let linear_predictor = dot_product(row, true_coefficients)
    let hazard = baseline_rate * double_exp(linear_predictor)
    let failure_uniform = clamp_positive(rng.uniform(), floor=0.0000000001)
    let failure = -double_log(failure_uniform) / clamp_positive(hazard)
    let censor_uniform = clamp_positive(rng.uniform(), floor=0.0000000001)
    let censor = if censor_rate <= 0.0 {
      failure + 1000000.0
    } else {
      -double_log(censor_uniform) / censor_rate
    }
    times.push(if failure < censor { failure } else { censor })
    events.push(if failure < censor { 1 } else { 0 })
    covariates.push(row)
    latent.push(failure)
  }
  {
    data: CoxData::new(times, events, covariates),
    true_coefficients,
    latent_failure_times: latent,
  }
}

///|
/// A reproducible sample used by benchmark documentation.
pub fn benchmark_fixture(sample_size : Int) -> SimulatedSurvival {
  simulate_survival({
    sample_size,
    baseline_rate: 0.06,
    censor_rate: 0.018,
    shape: 1.35,
    seed: 424242U,
  })
}

///|
/// Bootstrap confidence band for survival at fixed horizons.
pub(all) struct BootstrapKMResult {
  horizons : Array[Double]
  mean : Array[Double]
  lower : Array[Double]
  upper : Array[Double]
  replicates : Int
} derive(Eq, Debug)

///|
/// Deterministic nonparametric bootstrap for a vector of survival queries.
pub fn bootstrap_kaplan_meier(
  data : SurvivalData,
  horizons : Array[Double],
  replicates : Int,
  seed : UInt,
) -> BootstrapKMResult {
  if data.length() == 0 || replicates <= 0 {
    return {
      horizons,
      mean: zeros(horizons.length()),
      lower: zeros(horizons.length()),
      upper: zeros(horizons.length()),
      replicates: 0,
    }
  }
  let all_values : Array[Array[Double]] = []
  let rng = DeterministicRng::new(seed)
  for _ in 0..<replicates {
    let times = []
    let events = []
    for _ in 0..<data.length() {
      let index = (rng.next_u32() % data.length().reinterpret_as_uint()).reinterpret_as_int()
      let observation = data.get(index)
      times.push(observation.time)
      events.push(observation.event)
    }
    let estimate = kaplan_meier(SurvivalData::new(times, events))
    all_values.push(estimate.predict_many(horizons))
  }
  let mean = zeros(horizons.length())
  let lower = zeros(horizons.length())
  let upper = zeros(horizons.length())
  for j in 0..<horizons.length() {
    let values = []
    for replicate in all_values {
      values.push(replicate[j])
      mean[j] += replicate[j]
    }
    mean[j] /= replicates.to_double()
    lower[j] = quantile_doubles(values, 0.025).unwrap_or(0.0)
    upper[j] = quantile_doubles(values, 0.975).unwrap_or(0.0)
  }
  { horizons, mean, lower, upper, replicates }
}

///|
/// Return the empirical censoring fraction of a simulated sample.
pub fn SimulatedSurvival::censoring_fraction(
  self : SimulatedSurvival,
) -> Double {
  if self.data.length() == 0 {
    0.0
  } else {
    self.data.censored_count().to_double() / self.data.length().to_double()
  }
}

///|
/// Compare recovered latent failure mean and observed RMST.
pub(all) struct SimulationDiagnostics {
  sample_size : Int
  observed_events : Int
  censoring_fraction : Double
  latent_mean : Double
  observed_rmst : Double
  median_survival : Double?
} derive(Eq, Debug)

///|
/// Summarize a simulation for a reproducible validation report.
pub fn SimulatedSurvival::diagnostics(
  self : SimulatedSurvival,
  horizon : Double,
) -> SimulationDiagnostics {
  let latent_mean = mean_doubles(self.latent_failure_times)
  let km = kaplan_meier(self.data)
  {
    sample_size: self.data.length(),
    observed_events: self.data.event_count(),
    censoring_fraction: self.censoring_fraction(),
    latent_mean,
    observed_rmst: km.restricted_mean(horizon),
    median_survival: km.median(),
  }
}
