# Non-Probability Survey Inference

A simulation study evaluating bias-adjustment methods for non-probability survey samples, with a focus on when adjustment works and where its assumptions break down.

## Motivation

Non-probability samples are increasingly common in online surveys and platform-based data collection. They can be timely and cost-efficient, but raise challenges for population inference because inclusion probabilities are unknown and selection may depend on the outcome of interest.

This project asks:

> Under what conditions can adjustment methods reduce bias in non-probability survey samples, and where do their assumptions break down?

## Simulation Design

The simulation varies two dimensions:

**Selection mechanism** — ignorable (depends only on observed X) vs. non-ignorable (also depends directly on Y).

**Auxiliary information quality** — strong (X strongly predicts Y) vs. weak (X weakly predicts Y).

This yields four scenarios: `strong_ignorable`, `strong_nonignorable`, `weak_ignorable`, `weak_nonignorable`.

Each scenario uses a finite population of N = 100,000, a probability sample of n = 2,000, and a non-probability sample of n ≈ 3,000.

## Methods Evaluated

| Family | Role |
|---|---|
| Naive estimator | Baseline |
| Calibration / balancing weights | Main comparison |
| Propensity-score / IPW (logistic, RF, XGB) | Main comparison |
| Mass imputation | Main comparison |
| Doubly robust (logistic, RF, XGB) | Key highlight |
| Non-ignorable selection models | Planned |
| Sensitivity analysis / bounds | Planned |

## Current Status

- [x] Finite population simulation and sampling (`01_simulate_population.R`)
- [x] Baseline estimators under full-population auxiliary information (`02_estimate_methods.R`)
- [x] Reference-sample integration with probability + non-probability samples (`03_reference_sample_integration.R`)
- [x] Monte Carlo evaluation across 500 replicates (`04_monte_carlo_simulations.R`)
- [x] Benchmark against `{nonprobsvy}` package — single replicate (`05_nonprobsvy.R`)
- [x] ML propensity models: random forest and gradient boosting (`06_ml_propensity.R`)
- [ ] Non-ignorable selection extensions and sensitivity analysis

## Results Summary

**Scripts 02–03** established the baseline: adjustment works well under ignorable selection with strong auxiliary variables, partially reduces bias when selection is non-ignorable, and is unreliable when auxiliary variables are weak. Replacing exact population totals with a probability reference sample increases estimation variance but preserves the pattern.

**Script 04** (MC, 500 reps) confirmed these patterns hold across repeated samples. Calibration and mass imputation are numerically equivalent in the current linear setup.

**Script 05** (nonprobsvy benchmark) validated the hand-coded estimators: mass imputation and DR estimates match the `{nonprobsvy}` package to floating-point precision. GEE-calibrated IPW is close to hand-coded Hájek IPW; MLE logit IPW differs due to normalization conventions.

**Script 06** (ML propensity, 500 reps) compared logistic, random forest, and gradient-boosted propensity models for IPW and DR. Key findings:
- DR estimates are nearly identical across all three propensity models in every scenario, demonstrating the double-robustness property empirically.
- IPW is sensitive to propensity model choice. In `strong_ignorable` (where logistic is correctly specified), RF IPW shows higher bias than logistic IPW despite higher ESS, suggesting overly uniform weights. In `weak_ignorable`, RF IPW slightly outperforms logistic IPW.
- ML propensity models do not help under non-ignorable selection — residual bias remains at the same level regardless of propensity model.

## Planned Next Steps

1. Non-ignorable selection extensions
   - sensitivity parameters for outcome-dependent selection;
   - pseudo-likelihood estimators using reference probability samples;
   - compare IPW, DR, and prediction estimators as selection becomes increasingly non-ignorable.

2. Later-stage extensions
   - partial identification and bounds;
   - ML outcome models;
   - high-dimensional auxiliary information.

## Repository Structure

```text
nonprob-survey-inference/
├── README.md
├── nonprob-survey-inference.Rproj
├── scripts/
│   ├── 01_simulate_population.R
│   ├── 02_estimate_methods.R
│   ├── 03_reference_sample_integration.R
│   ├── 04_monte_carlo_simulations.R
│   ├── 05_nonprobsvy.R
│   └── 06_ml_propensity.R
├── notes/
│   ├── Simulation_results_round1.md
│   ├── simulation_results_round2_reference_integration.md
│   ├── simulation_results_round3_nonprobsvy_benchmark.md
│   └── simulation_results_round3_ml_propensity.md
├── R/
├── data/
│   └── processed/
├── outputs/
│   └── tables/
└── references/
```
