# MoonBit Survival Analysis Toolkit

moon-survival is a practical MoonBit library for right-censored survival
data, designed for medical follow-up, reliability engineering, software
failure analysis, and user churn.

It includes Kaplan-Meier and Nelson-Aalen estimators, Cox proportional hazards,
competing risks, log-rank tests, RMST, power planning, resampling, reliability
models, CSV/time-varying data, incremental processing, reports, diagnostics,
and a validated Markov cohort model.

## Run it

~~~bash
moon add Lyl66655/moon-survival
moon run --target native cmd/example
moon run --target native --release cmd/benchmark
~~~

The deterministic benchmark fixture has 5,000 observations, 4,403 events, and
11.94% censoring. In the local Windows run on 2026-08-18, KM took 1 ms,
Nelson-Aalen 1 ms, and a 500-row three-covariate Cox fit converged in four
iterations and 54 ms. Timing varies by machine.

## Verification

~~~bash
moon fmt --check
moon check --target all --deny-warn
moon test --target wasm-gc --deny-warn
moon test --target native --deny-warn
moon info
~~~

The local regression suite has 109 passing tests, including boundary cases for
empty/all-censored data, ties, invalid shapes, singular matrices, CSV quoting,
numerical extremes, and deterministic benchmark fixtures. CI runs on all three
major operating systems using the latest stable MoonBit toolchain.

## Links and license

- GitHub: https://github.com/Lyl66655/moon-survival
- GitLink: https://www.gitlink.org.cn/ylyl/moon-survival
- Mooncakes: Lyl66655/moon-survival
- License: Apache-2.0

The project uses only MoonBit core packages and contains no copied
third-party source.
