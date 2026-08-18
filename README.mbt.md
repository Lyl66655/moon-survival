# MoonBit Survival Analysis Toolkit

moon-survival is a practical survival-analysis toolkit for MoonBit. It
supports medical follow-up, reliability engineering, software incident
analysis, and user churn studies with right-censored observations.

## Implemented scope

- Kaplan-Meier curves with Greenwood standard errors, confidence intervals,
  median survival, quantiles, RMST, weighted observations, and CSV export.
- Nelson-Aalen cumulative hazards, variance, confidence intervals, and
  survival transformations.
- Cox proportional hazards with Newton fitting, ridge stabilization, Breslow
  baseline hazards, predictions, residuals, and model diagnostics.
- Competing risks (Aalen-Johansen), log-rank tests, RMST differences, power
  planning, and bootstrap/jackknife resampling.
- CSV ingestion, time-varying covariates, incremental accumulators, grouped
  analysis, reliability models, simulation fixtures, and Markdown reports.
- A validated finite-state Markov model and clinical trial reporting helpers
  for downstream decision analysis.

## Quick start

~~~bash
moon add Lyl66655/moon-survival
moon run --target native cmd/example
~~~

The runnable example covers Kaplan-Meier, Nelson-Aalen, Cox fitting, and
competing risks. The benchmark fixture is reproducible:

~~~bash
moon run --target native --release cmd/benchmark
~~~

On the local Windows benchmark run recorded on 2026-08-18, 5,000 simulated
observations produced 4,403 events and 11.94% censoring. Kaplan-Meier took
1 ms, Nelson-Aalen took 1 ms, and a 500-row, three-covariate Cox fit converged
in four iterations and 54 ms. Timing is machine-dependent; the fixture and
reported sample statistics are deterministic.

## Quality checks

~~~bash
moon fmt --check
moon check --target all --deny-warn
moon test --target wasm-gc --deny-warn
moon test --target native --deny-warn
moon info
~~~

The repository CI runs the same checks on Ubuntu, macOS, and Windows with the
latest stable MoonBit toolchain. The current local regression suite contains
109 passing tests, including malformed input, empty data, all-censored data,
ties, singular matrices, quoting, numerical extremes, and deterministic
performance fixtures.

## Package and license

- GitHub: https://github.com/Lyl66655/moon-survival
- GitLink mirror: https://www.gitlink.org.cn/ylyl/moon-survival
- Mooncakes package: Lyl66655/moon-survival
- License: Apache-2.0

The implementation depends only on MoonBit core packages. No third-party
source was copied into this repository.
