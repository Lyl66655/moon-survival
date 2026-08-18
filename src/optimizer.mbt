///|
/// Optimization controls shared by reusable model-fitting routines.
pub(all) struct OptimizerOptions {
  max_iterations : Int
  tolerance : Double
  learning_rate : Double
  memory : Int
} derive(Eq, Debug)

///|
/// Default robust optimization controls.
pub fn default_optimizer_options() -> OptimizerOptions {
  { max_iterations: 100, tolerance: 0.000001, learning_rate: 1.0, memory: 5 }
}

///|
/// Generic optimization result.
pub(all) struct OptimizationResult {
  parameters : Array[Double]
  objective : Double
  gradient : Array[Double]
  iterations : Int
  evaluations : Int
  converged : Bool
} derive(Eq, Debug)

///|
/// A numerically stable finite-difference gradient for scalar objectives.
pub fn finite_difference_gradient(
  objective : (Array[Double]) -> Double,
  parameters : Array[Double],
  step? : Double = 0.000001,
) -> Array[Double] {
  let gradient = zeros(parameters.length())
  for i in 0..<parameters.length() {
    let plus = []
    let minus = []
    for j in 0..<parameters.length() {
      plus.push(if i == j { parameters[j] + step } else { parameters[j] })
      minus.push(if i == j { parameters[j] - step } else { parameters[j] })
    }
    gradient[i] = (objective(plus) - objective(minus)) / (2.0 * step)
  }
  gradient
}

///|
/// Gradient ascent with Armijo-style backtracking.
pub fn gradient_ascent(
  initial : Array[Double],
  objective : (Array[Double]) -> (Double, Array[Double]),
  options : OptimizerOptions,
) -> OptimizationResult {
  let mut parameters = vector_scale(initial, 1.0)
  let mut evaluation = objective(parameters)
  let mut evaluations = 1
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<options.max_iterations {
    let gradient = evaluation.1
    if vector_norm(gradient) <= options.tolerance {
      converged = true
      iterations = iteration
      break
    }
    let mut step_size = options.learning_rate
    let mut candidate = vector_add(
      parameters,
      vector_scale(gradient, step_size),
    )
    let mut candidate_evaluation = objective(candidate)
    evaluations += 1
    while candidate_evaluation.0 < evaluation.0 && step_size > 0.0000001 {
      step_size *= 0.5
      candidate = vector_add(parameters, vector_scale(gradient, step_size))
      candidate_evaluation = objective(candidate)
      evaluations += 1
    }
    parameters = candidate
    evaluation = candidate_evaluation
    iterations = iteration + 1
  }
  {
    parameters,
    objective: evaluation.0,
    gradient: evaluation.1,
    iterations,
    evaluations,
    converged,
  }
}

///|
/// Limited-memory BFGS maximization for smooth log likelihoods. The two-loop
/// recursion stores only recent parameter and gradient differences, making it
/// suitable for moderate MoonBit data-science models without a dense Hessian.
pub fn lbfgs_maximize(
  initial : Array[Double],
  objective : (Array[Double]) -> (Double, Array[Double]),
  options : OptimizerOptions,
) -> OptimizationResult {
  let mut parameters = vector_scale(initial, 1.0)
  let mut evaluation = objective(parameters)
  let mut evaluations = 1
  let s_history : Array[Array[Double]] = []
  let y_history : Array[Array[Double]] = []
  let rho_history : Array[Double] = []
  let mut iterations = 0
  let mut converged = false
  for iteration in 0..<options.max_iterations {
    let gradient = evaluation.1
    if vector_norm(gradient) <= options.tolerance {
      converged = true
      iterations = iteration
      break
    }
    let direction = lbfgs_direction(gradient, s_history, y_history, rho_history)
    let mut step_size = options.learning_rate
    let mut candidate = vector_add(
      parameters,
      vector_scale(direction, step_size),
    )
    let mut candidate_evaluation = objective(candidate)
    evaluations += 1
    while candidate_evaluation.0 < evaluation.0 && step_size > 0.0000001 {
      step_size *= 0.5
      candidate = vector_add(parameters, vector_scale(direction, step_size))
      candidate_evaluation = objective(candidate)
      evaluations += 1
    }
    let s = vector_subtract(candidate, parameters)
    let y = vector_subtract(candidate_evaluation.1, gradient)
    let curvature = dot_product(s, y)
    if curvature > survival_epsilon {
      if s_history.length() >= options.memory && s_history.length() > 0 {
        ignore(s_history.remove(0))
        ignore(y_history.remove(0))
        ignore(rho_history.remove(0))
      }
      s_history.push(s)
      y_history.push(y)
      rho_history.push(1.0 / curvature)
    }
    parameters = candidate
    evaluation = candidate_evaluation
    iterations = iteration + 1
  }
  {
    parameters,
    objective: evaluation.0,
    gradient: evaluation.1,
    iterations,
    evaluations,
    converged,
  }
}

///|
/// Build an approximate inverse-Hessian direction using the L-BFGS two-loop
/// recursion. The final sign is positive because the public routine maximizes.
fn lbfgs_direction(
  gradient : Array[Double],
  s_history : Array[Array[Double]],
  y_history : Array[Array[Double]],
  rho_history : Array[Double],
) -> Array[Double] {
  let mut direction = vector_scale(gradient, 1.0)
  let alpha = zeros(s_history.length())
  for index in 0..<s_history.length() {
    let reverse = s_history.length() - 1 - index
    alpha[reverse] = rho_history[reverse] *
      dot_product(s_history[reverse], direction)
    for j in 0..<direction.length() {
      direction[j] -= alpha[reverse] * y_history[reverse][j]
    }
  }
  if y_history.length() > 0 {
    let last = y_history.length() - 1
    let yy = dot_product(y_history[last], y_history[last])
    let sy = dot_product(s_history[last], y_history[last])
    if yy > survival_epsilon {
      direction = vector_scale(direction, sy / yy)
    }
  }
  for index in 0..<s_history.length() {
    let beta = rho_history[index] * dot_product(y_history[index], direction)
    for j in 0..<direction.length() {
      direction[j] += s_history[index][j] * (alpha[index] - beta)
    }
  }
  direction
}

///|
/// Numerically maximize a scalar objective by wrapping a finite-difference
/// gradient around the L-BFGS driver.
pub fn lbfgs_maximize_scalar(
  initial : Array[Double],
  objective : (Array[Double]) -> Double,
  options : OptimizerOptions,
) -> OptimizationResult {
  lbfgs_maximize(
    initial,
    parameters => {
      (objective(parameters), finite_difference_gradient(objective, parameters))
    },
    options,
  )
}

///|
/// Minimize a smooth scalar objective by maximizing its negation.
pub fn lbfgs_minimize_scalar(
  initial : Array[Double],
  objective : (Array[Double]) -> Double,
  options : OptimizerOptions,
) -> OptimizationResult {
  let result = lbfgs_maximize_scalar(
    initial,
    parameters => -objective(parameters),
    options,
  )
  {
    parameters: result.parameters,
    objective: -result.objective,
    gradient: vector_scale(result.gradient, -1.0),
    iterations: result.iterations,
    evaluations: result.evaluations,
    converged: result.converged,
  }
}
