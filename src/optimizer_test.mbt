///|
test "finite difference and gradient ascent optimize a concave objective" {
  let objective = (values : Array[Double]) => {
    let difference = values[0] - 3.0
    -(difference * difference)
  }
  let gradient = finite_difference_gradient(objective, [1.0])
  inspect(gradient[0] > 0.0, content="true")
  let result = gradient_ascent(
    [0.0],
    (values : Array[Double]) => {
      let difference = values[0] - 3.0
      (-difference * difference, [-2.0 * difference])
    },
    { max_iterations: 100, tolerance: 0.00001, learning_rate: 0.2, memory: 4 },
  )
  inspect(
    result.parameters[0] > 2.9 && result.parameters[0] < 3.1,
    content="true",
  )
  inspect(result.objective > -0.02, content="true")
}

///|
test "lbfgs scalar minimization finds a two dimensional optimum" {
  let objective = (values : Array[Double]) => {
    let x = values[0] - 2.0
    let y = values[1] + 1.0
    x * x + 2.0 * y * y
  }
  let result = lbfgs_minimize_scalar([0.0, 0.0], objective, {
    max_iterations: 80,
    tolerance: 0.00001,
    learning_rate: 0.5,
    memory: 5,
  })
  inspect(
    result.parameters[0] > 1.8 && result.parameters[0] < 2.2,
    content="true",
  )
  inspect(
    result.parameters[1] < -0.8 && result.parameters[1] > -1.2,
    content="true",
  )
  inspect(result.objective < 0.1, content="true")
}

///|
test "optimizer respects zero iteration boundary" {
  let result = gradient_ascent(
    [1.0],
    (values : Array[Double]) => (values[0], [1.0]),
    { max_iterations: 0, tolerance: 0.00001, learning_rate: 1.0, memory: 2 },
  )
  inspect(result.iterations, content="0")
  inspect(result.evaluations, content="1")
  inspect(result.parameters[0], content="1")
}
