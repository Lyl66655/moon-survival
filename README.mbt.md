# MoonBit Survival Analysis Toolkit (moon-survival)

A comprehensive survival analysis toolkit for MoonBit, designed for medical follow-ups, equipment lifetime estimation, software failure rates, and user churn analysis.

## Features

- **Kaplan-Meier Estimator**: Non-parametric survival curves with confidence intervals (Greenwood's formula).
- **Nelson-Aalen Estimator**: Cumulative hazard rate estimation.
- **Cox Proportional Hazards**: Basic foundation for partial likelihood evaluation with covariates.
- **Right-Censored Data**: Full support for right-censored observation data.

## Usage

```mbt check
test "Kaplan-Meier Example" {
  // 1 = Event occurred, 0 = Censored
  let times = [10.0, 20.0, 30.0, 40.0, 50.0]
  let events = [1, 0, 1, 1, 0]
  
  // Initialize survival data
  let data = @moon_survival.SurvivalData::new(times, events)
  
  // Estimate Kaplan-Meier survival curve
  let km = @moon_survival.kaplan_meier(data)
  
  // Predict survival probability at time 25.0
  let p = km.predict(25.0)
  inspect(p > 0.0, content="true")
}

test "Nelson-Aalen Example" {
  let times = [5.0, 10.0, 15.0]
  let events = [1, 1, 1]
  
  let data = @moon_survival.SurvivalData::new(times, events)
  let na = @moon_survival.nelson_aalen(data)
  
  // Predict cumulative hazard at time 12.0
  let h = na.predict(12.0)
  inspect(h > 0.0, content="true")
}
```

## Contributing

We welcome contributions to expand the toolkit (e.g., competing risks, time-varying covariates, etc.).
