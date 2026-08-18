# moon-survival

Survival analysis toolkit for MoonBit, with native support for right-censored
data. It is designed for medical follow-up, reliability engineering, software
failure analysis, and user churn studies.

![MoonBit quality checks](https://github.com/Lyl66655/moon-survival/actions/workflows/check.yml/badge.svg)

## Project positioning

moon-survival provides reusable statistical building blocks for estimating
time-to-event behavior across WebAssembly, Native, and JavaScript targets.
The library is dependency-light and keeps data validation, estimators, model
diagnostics, reporting, and reproducible simulation fixtures in one MoonBit
module.

## Core capabilities

- Right-censored and weighted survival observations with validation and
  summary statistics.
- Kaplan-Meier with Greenwood standard errors, confidence intervals, quantiles,
  median survival, restricted mean survival time, and curve export.
- Nelson-Aalen cumulative hazards with variance, confidence intervals, and
  survival transformations.
- Cox proportional hazards with Newton fitting, ridge stabilization, Breslow
  baseline hazards, predictions, residuals, and proportional-hazards diagnostics.
- Competing risks with Aalen-Johansen cumulative incidence, log-rank tests,
  pairwise comparisons, and RMST differences.
- Reliability models, power planning, bootstrap/jackknife resampling, grouped
  analysis, model selection, and calibration metrics.
- CSV ingestion/export, time-varying covariates, incremental accumulation,
  deterministic simulation, Markdown/CSV reports, ASCII visualization, and
  finite-state Markov cohort analysis.

## Quick start

~~~bash
moon add Lyl66655/moon-survival
moon run --target native cmd/example
~~~

The example executable demonstrates Kaplan-Meier, Nelson-Aalen, Cox, and
competing-risk workflows. A package can import the library package as:

~~~mbt
import {
  "Lyl66655/moon-survival/src",
}
~~~

## CLI and runnable commands

The repository provides two executable MoonBit entrypoints:

~~~bash
# Run the end-to-end example
moon run --target native cmd/example

# Run the deterministic performance fixture
moon run --target native --release cmd/benchmark
~~~

The benchmark command is the project CLI for repeatable local performance
checks; the library itself is consumed through MoonBit package APIs.

## Architecture

~~~text
moon.mod                         module metadata and package version
src/moon.pkg                     public survival-analysis package
src/*.mbt                        estimators, models, diagnostics, reports
src/pkg.generated.mbti           generated public API summary
cmd/example/                     runnable API example
cmd/benchmark/                   deterministic benchmark executable
.github/workflows/check.yml      cross-platform quality checks
LICENSE                           Apache-2.0 license text
~~~

The source package is organized by responsibility: data and numerical
utilities are the foundation; non-parametric estimators, regression models,
reliability and competing-risk models build on them; diagnostics, reporting,
visualization, and command entrypoints consume those public APIs.

## Benchmark

The deterministic fixture contains 5,000 observations, 4,403 events, and
11.94% censoring. A reference run on 2026-08-18 with MoonBit 0.1.20260814 and
moonc 0.10.8 measured:

| Workload | Result |
|---|---:|
| Kaplan-Meier, 5,000 observations | 1 ms |
| Nelson-Aalen, 5,000 observations | 1 ms |
| Cox, 500 rows and 3 covariates | 4 iterations, 54 ms |

The fixture, event counts, censoring fraction, convergence status, and RMST
are deterministic. Timings depend on the host machine.

## Tests

~~~bash
moon fmt --check
moon check --target all --deny-warn
moon test --target wasm-gc --deny-warn
moon test --target native --deny-warn
moon info
~~~

The regression suite contains 109 passing tests covering empty data,
all-censored data, ties, invalid events, malformed matrices, singular
systems, CSV quoting, numerical extremes, model edge cases, and deterministic
performance fixtures.

## CI

GitHub Actions runs on Ubuntu, macOS, and Windows. Each job installs the latest
stable MoonBit toolchain from the official installer and runs formatting,
all-target compilation checks, metadata generation, generated-source drift
detection, WebAssembly GC tests, Native tests, JavaScript tests, and all-target
builds.

## License

Apache-2.0. See [LICENSE](LICENSE).

## Resources

- [GitHub repository](https://github.com/Lyl66655/moon-survival)
- [GitLink mirror](https://www.gitlink.org.cn/ylyl/moon-survival)
- Mooncakes package: Lyl66655/moon-survival
