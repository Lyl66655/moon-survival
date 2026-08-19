///|
/// Observation types supported by interval-censored survival analysis.
pub(all) enum IntervalObservationKind {
  ExactEvent
  RightCensored
  LeftCensored
  IntervalEvent
} derive(Eq, Debug)

///|
/// One interval-censored event-time observation.
pub(all) struct IntervalObservation {
  lower : Double
  upper : Double
  kind : IntervalObservationKind
} derive(Eq, Debug)

///|
/// Validation failures for interval-censored records.
pub(all) enum IntervalDataError {
  NegativeTime
  NonFiniteTime
  ReversedInterval
  InvalidRightCensoring
  InvalidLeftCensoring
  EmptyData
} derive(Eq, Debug)

///|
/// Validated interval-censored data.
pub(all) struct IntervalSurvivalData {
  observations : Array[IntervalObservation]
} derive(Eq, Debug)

///|
/// A Turnbull non-parametric maximum-likelihood fit.
pub(all) struct TurnbullFit {
  grid : Array[Double]
  masses : Array[Double]
  iterations : Int
  converged : Bool
  log_likelihood : Double
  sample_size : Int
} derive(Eq, Debug)

///|
/// Create one exact event observation.
pub fn interval_exact(time : Double) -> IntervalObservation {
  { lower: time, upper: time, kind: ExactEvent }
}

///|
/// Create one right-censored observation.
pub fn interval_right_censored(last_seen : Double) -> IntervalObservation {
  { lower: last_seen, upper: last_seen, kind: RightCensored }
}

///|
/// Create one left-censored event observation.
pub fn interval_left_censored(first_seen : Double) -> IntervalObservation {
  { lower: 0.0, upper: first_seen, kind: LeftCensored }
}

///|
/// Create one interval event observation.
pub fn interval_event(lower : Double, upper : Double) -> IntervalObservation {
  { lower, upper, kind: IntervalEvent }
}

///|
/// Validate an interval observation.
pub fn interval_observation_is_valid(
  observation : IntervalObservation,
) -> Result[IntervalObservation, IntervalDataError] {
  if observation.lower.is_nan() ||
    observation.upper.is_nan() ||
    observation.lower.is_inf() ||
    observation.upper.is_inf() {
    return Err(NonFiniteTime)
  }
  if observation.lower < 0.0 || observation.upper < 0.0 {
    return Err(NegativeTime)
  }
  match observation.kind {
    ExactEvent =>
      if (observation.lower - observation.upper).abs() > survival_epsilon {
        Err(ReversedInterval)
      } else {
        Ok(observation)
      }
    RightCensored => Ok(observation)
    LeftCensored =>
      if observation.lower > survival_epsilon {
        Err(InvalidLeftCensoring)
      } else {
        Ok(observation)
      }
    IntervalEvent =>
      if observation.upper < observation.lower {
        Err(ReversedInterval)
      } else {
        Ok(observation)
      }
  }
}

///|
/// Validate and sort interval-censored data.
pub fn IntervalSurvivalData::try_new(
  observations : Array[IntervalObservation],
) -> Result[IntervalSurvivalData, IntervalDataError] {
  if observations.length() == 0 {
    return Err(EmptyData)
  }
  let copied : Array[IntervalObservation] = []
  for observation in observations {
    match interval_observation_is_valid(observation) {
      Ok(value) => copied.push(value)
      Err(error) => return Err(error)
    }
  }
  copied.sort_by(fn(a, b) {
    if a.lower < b.lower {
      -1
    } else if a.lower > b.lower {
      1
    } else {
      0
    }
  })
  Ok({ observations: copied })
}

///|
/// Construct interval data and abort on invalid input.
pub fn IntervalSurvivalData::new(
  observations : Array[IntervalObservation],
) -> IntervalSurvivalData {
  match IntervalSurvivalData::try_new(observations) {
    Ok(value) => value
    Err(_) => panic()
  }
}

///|
/// Number of interval observations.
pub fn IntervalSurvivalData::length(self : IntervalSurvivalData) -> Int {
  self.observations.length()
}

///|
/// Access one interval observation.
pub fn IntervalSurvivalData::get(
  self : IntervalSurvivalData,
  index : Int,
) -> IntervalObservation {
  self.observations[index]
}

///|
/// Copy observations out of the validated data container.
pub fn IntervalSurvivalData::to_array(
  self : IntervalSurvivalData,
) -> Array[IntervalObservation] {
  let result : Array[IntervalObservation] = []
  for observation in self.observations {
    result.push(observation)
  }
  result
}

///|
/// Count exact, censored, and interval event records.
pub fn interval_kind_counts(data : IntervalSurvivalData) -> Array[Int] {
  let result = [0, 0, 0, 0]
  for observation in data.observations {
    match observation.kind {
      ExactEvent => result[0] += 1
      RightCensored => result[1] += 1
      LeftCensored => result[2] += 1
      IntervalEvent => result[3] += 1
    }
  }
  result
}

///|
/// Return all finite event-time bounds used to construct an NPMLE grid.
pub fn interval_time_grid(data : IntervalSurvivalData) -> Array[Double] {
  let values : Array[Double] = [0.0]
  for observation in data.observations {
    values.push(observation.lower)
    if observation.kind != RightCensored {
      values.push(observation.upper)
    }
  }
  values.sort()
  let result : Array[Double] = []
  for value in values {
    if result.length() == 0 ||
      (result[result.length() - 1] - value).abs() > survival_epsilon {
      result.push(value)
    }
  }
  result
}

///|
/// Fit the Turnbull interval-censored NPMLE using EM.
pub fn turnbull_fit(
  data : IntervalSurvivalData,
  max_iterations? : Int = 500,
  tolerance? : Double = 0.000001,
) -> TurnbullFit {
  let grid = interval_time_grid(data)
  if grid.length() == 0 {
    return {
      grid: [],
      masses: [],
      iterations: 0,
      converged: false,
      log_likelihood: 0.0,
      sample_size: 0,
    }
  }
  let masses = []
  for _ in 0..<grid.length() {
    masses.push(1.0 / grid.length().to_double())
  }
  let current = masses
  let limit = if max_iterations <= 0 { 1 } else { max_iterations }
  let target = if tolerance <= 0.0 { 0.000001 } else { tolerance }
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<limit {
    let expected = zeros(grid.length())
    let mut valid_rows = 0
    for observation in data.observations {
      let compatible = interval_compatible_indices(observation, grid)
      let probability = interval_mass(current, compatible)
      if compatible.length() == 0 || probability <= survival_epsilon {
        continue
      }
      valid_rows += 1
      for index in compatible {
        expected[index] += current[index] / probability
      }
    }
    if valid_rows == 0 {
      break
    }
    let mut difference = 0.0
    let denominator = valid_rows.to_double()
    for i in 0..<expected.length() {
      let next = expected[i] / denominator
      difference += (next - current[i]).abs()
      current[i] = next
    }
    iterations = iteration + 1
    if difference <= target {
      converged = true
      break
    }
  }
  {
    grid,
    masses: current,
    iterations,
    converged,
    log_likelihood: interval_turnbull_log_likelihood(data, grid, current),
    sample_size: data.length(),
  }
}

///|
/// Return the fitted CDF at a time.
pub fn TurnbullFit::cdf_at(self : TurnbullFit, time : Double) -> Double {
  let mut total = 0.0
  for i in 0..<self.grid.length() {
    if self.grid[i] <= time {
      total += self.masses[i]
    }
  }
  clamp_probability(total)
}

///|
/// Return the fitted survival probability at a time.
pub fn TurnbullFit::survival_at(self : TurnbullFit, time : Double) -> Double {
  1.0 - self.cdf_at(time)
}

///|
/// Return the discrete mass associated with a grid point.
pub fn TurnbullFit::mass_at(self : TurnbullFit, time : Double) -> Double {
  for i in 0..<self.grid.length() {
    if (self.grid[i] - time).abs() <= survival_epsilon {
      return self.masses[i]
    }
  }
  0.0
}

///|
/// Return the first grid time with CDF at least the requested probability.
pub fn TurnbullFit::quantile(
  self : TurnbullFit,
  probability : Double,
) -> Double? {
  let target = clamp_probability(probability)
  let mut cumulative = 0.0
  for i in 0..<self.grid.length() {
    cumulative += self.masses[i]
    if cumulative >= target {
      return Some(self.grid[i])
    }
  }
  None
}

///|
/// Numerically integrate the Turnbull survival curve.
pub fn TurnbullFit::restricted_mean(
  self : TurnbullFit,
  horizon : Double,
) -> Double {
  if horizon <= 0.0 {
    return 0.0
  }
  let mut area = 0.0
  let mut previous = 0.0
  let mut survival = 1.0
  for i in 0..<self.grid.length() {
    let time = self.grid[i]
    if time >= horizon {
      break
    }
    area += survival * (time - previous)
    survival -= self.masses[i]
    previous = time
  }
  area + survival * (horizon - previous)
}

///|
/// Return a discrete hazard at a grid point.
pub fn TurnbullFit::hazard_at(self : TurnbullFit, time : Double) -> Double {
  let mut index = -1
  for i in 0..<self.grid.length() {
    if (self.grid[i] - time).abs() <= survival_epsilon {
      index = i
      break
    }
  }
  if index < 0 {
    return 0.0
  }
  let mut prior = 0.0
  for i in 0..<index {
    prior += self.masses[i]
  }
  let at_risk = 1.0 - prior
  if at_risk <= survival_epsilon {
    0.0
  } else {
    self.masses[index] / at_risk
  }
}

///|
/// Produce a compact CSV curve for external plotting.
pub fn TurnbullFit::to_csv(self : TurnbullFit) -> String {
  let builder = StringBuilder::new()
  builder.write_string("time,mass,cdf,survival,hazard\n")
  for time in self.grid {
    builder.write_string(
      "\{time},\{self.mass_at(time)},\{self.cdf_at(time)},\{self.survival_at(time)},\{self.hazard_at(time)}\n",
    )
  }
  builder.to_string()
}

///|
/// Convert ordinary right-censored data to interval observations.
pub fn interval_data_from_survival(data : SurvivalData) -> IntervalSurvivalData {
  let observations : Array[IntervalObservation] = []
  for observation in data.to_array() {
    if observation.event == 1 {
      observations.push(interval_exact(observation.time))
    } else {
      observations.push(interval_right_censored(observation.time))
    }
  }
  IntervalSurvivalData::new(observations)
}

///|
/// Estimate the probability an interval contains the event.
pub fn interval_event_probability(
  fit : TurnbullFit,
  lower : Double,
  upper : Double,
) -> Double {
  if upper < lower {
    return 0.0
  }
  clamp_probability(fit.cdf_at(upper) - fit.cdf_at(lower))
}

///|
/// Estimate the probability of surviving beyond a visit schedule.
pub fn interval_visit_survival(
  fit : TurnbullFit,
  visits : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  for visit in visits {
    result.push(fit.survival_at(visit))
  }
  result
}

///|
/// Return whether the fitted masses form a valid probability distribution.
pub fn turnbull_is_valid(fit : TurnbullFit) -> Bool {
  let mut total = 0.0
  for mass in fit.masses {
    if mass < -survival_epsilon || mass.is_nan() || mass.is_inf() {
      return false
    }
    total += mass
  }
  (total - 1.0).abs() <= 0.00001
}

///|
fn interval_compatible_indices(
  observation : IntervalObservation,
  grid : Array[Double],
) -> Array[Int] {
  let result : Array[Int] = []
  for i in 0..<grid.length() {
    let time = grid[i]
    let compatible = match observation.kind {
      ExactEvent => (time - observation.lower).abs() <= survival_epsilon
      RightCensored => time >= observation.lower
      LeftCensored => time <= observation.upper
      IntervalEvent => time >= observation.lower && time <= observation.upper
    }
    if compatible {
      result.push(i)
    }
  }
  result
}

///|
fn interval_mass(masses : Array[Double], indices : Array[Int]) -> Double {
  let mut total = 0.0
  for index in indices {
    total += masses[index]
  }
  total
}

///|
fn interval_turnbull_log_likelihood(
  data : IntervalSurvivalData,
  grid : Array[Double],
  masses : Array[Double],
) -> Double {
  let mut result = 0.0
  for observation in data.observations {
    let indices = interval_compatible_indices(observation, grid)
    result += double_log(interval_mass(masses, indices))
  }
  result
}
