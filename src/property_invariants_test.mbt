///|
test "kaplan meier is monotone for deterministic grids" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0, 5.0], [1, 0, 1, 0, 1])
  let fit = kaplan_meier(data)
  let grid = linspace(0.0, 6.0, 31)
  inspect(curve_is_monotone(fit, grid), content="true")
  for value in fit.predict_many(grid) {
    if value < 0.0 || value > 1.0 {
      fail("survival left probability bounds")
    }
  }
}

///|
test "cumulative hazard is monotone and nonnegative" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 1, 0, 1])
  let fit = nelson_aalen(data)
  let grid = linspace(0.0, 5.0, 20)
  let mut previous = 0.0
  for value in fit.predict_many(grid) {
    if value < previous {
      fail("hazard decreased")
    }
    previous = value
  }
  inspect(previous > 0.0, content="true")
}

///|
test "competing incidence identity holds over all steps" {
  let data = CompetingRiskData::new([1.0, 2.0, 3.0, 4.0, 5.0], [1, 2, 0, 1, 2])
  let fit = cumulative_incidence(data)
  for step in fit.steps {
    let total = step.survival +
      step.cumulative_incidence[0] +
      step.cumulative_incidence[1]
    if !double_is_close(total, 1.0) {
      fail("incidence identity failed")
    }
  }
  inspect(fit.predict_all(5.0).length(), content="2")
}

///|
test "weibull and exponential survival are bounded" {
  let data = benchmark_fixture(50).data
  let exponential = exponential_fit(data)
  let weibull = weibull_fit(data)
  for time in linspace(0.0, 50.0, 20) {
    let e = exponential.survival_at(time)
    let w = weibull.survival_at(time)
    if e < 0.0 || e > 1.0 || w < 0.0 || w > 1.0 {
      fail("parametric survival out of bounds")
    }
  }
  inspect(exponential.mean_time() >= 0.0, content="true")
}

///|
test "risk scores and calibration remain bounded" {
  let data = SurvivalData::new([1.0, 2.0, 3.0, 4.0], [1, 0, 1, 0])
  let scores = [0.9, 0.7, 0.5, 0.1]
  let concordance = concordance_index(data, scores)
  inspect(
    concordance.c_index >= 0.0 && concordance.c_index <= 1.0,
    content="true",
  )
  for bin in calibration_curve(data, [0.1, 0.3, 0.6, 0.9], 3.0, 4) {
    if bin.observed_event_fraction < 0.0 || bin.observed_event_fraction > 1.0 {
      fail("calibration fraction out of bounds")
    }
  }
}
