# moon-survival

Survival analysis toolkit for MoonBit with native support for right-censored
data. It targets medical follow-up, reliability engineering, software failure
analysis, and user churn studies.

## Project positioning

The library provides reusable time-to-event estimators, regression models,
diagnostics, reports, and deterministic simulation fixtures for WebAssembly,
Native, and JavaScript targets.

## Core capabilities

- Validated right-censored and weighted observations.
- Kaplan-Meier with Greenwood confidence intervals, quantiles, median survival,
  restricted mean survival time, and export.
- Nelson-Aalen cumulative hazards with variance and confidence intervals.
- Cox proportional hazards with Newton fitting, ridge stabilization, Breslow
  baseline hazards, predictions, residuals, and diagnostics.
- Competing risks, Aalen-Johansen cumulative incidence, log-rank, RMST,
  reliability models, power planning, resampling, and grouped analysis.
- CSV and time-varying data, incremental processing, model selection,
  calibration, Markdown/CSV reports, ASCII visualization, and Markov cohorts.

## Quick start

~~~bash
moon add Lyl66655/moon-survival
moon run --target native cmd/example
~~~

For a reproducible performance run:

~~~bash
moon run --target native --release cmd/benchmark
~~~

## CLI

The repository executables are `cmd/example` for an end-to-end API example and
`cmd/benchmark` for the deterministic performance fixture. The library is
consumed through MoonBit package APIs:

~~~mbt
import {
  "Lyl66655/moon-survival/src",
}
~~~

## Architecture

~~~text
src/           public survival-analysis package
cmd/example/   runnable API example
cmd/benchmark/ deterministic benchmark executable
~~~

The source is split by responsibility: data and numerical foundations,
non-parametric estimators, regression and reliability models, competing risks,
diagnostics, reporting, and visualization.

## Benchmark

The deterministic fixture contains 5,000 observations, 4,403 events, and
11.94% censoring. A reference run on 2026-08-18 with MoonBit 0.1.20260814 and
moonc 0.10.8 measured 1 ms for Kaplan-Meier, 1 ms for Nelson-Aalen, and four
iterations / 54 ms for a 500-row, three-covariate Cox fit. Timings depend on
the host machine.

## Tests and CI

~~~bash
moon fmt --check
moon check --target all --deny-warn
moon test --target wasm-gc --deny-warn
moon test --target native --deny-warn
moon info
~~~

The regression suite contains 109 passing tests, including malformed input,
empty and all-censored data, ties, singular matrices, CSV quoting, numerical
extremes, and deterministic performance fixtures. GitHub Actions runs the
same quality gates on Ubuntu, macOS, and Windows with the latest stable
MoonBit toolchain.

## License

Apache-2.0. See [LICENSE](LICENSE).

## Resources

- [GitHub repository](https://github.com/Lyl66655/moon-survival)
- [GitLink mirror](https://www.gitlink.org.cn/ylyl/moon-survival)
- Mooncakes package: Lyl66655/moon-survival
