# Non-Probability Survey Inference

This repository contains an early-stage methodological project on inference from non-probability survey samples. 
The project uses simulation studies to evaluate how different bias-adjustment methods perform under varying assumptions about sample selection and auxiliary information.

The long-term goal is to develop a transparent and reproducible workflow for assessing when non-probability samples can support reliable social and official statistics.

## Motivation

Non-probability samples are increasingly used in online surveys, opt-in panels, platform-based data collection, and mixed-source statistical systems. 
These data sources can be timely and cost-efficient, but they raise important challenges for population inference because sample inclusion probabilities are unknown and selection may be systematically related to the outcome of interest.

This project asks:

> Under what conditions can adjustment methods reduce bias in non-probability survey samples, and where do their assumptions break down?

The project is motivated by debates in survey methodology and official statistics about how to responsibly integrate new forms of data collection while maintaining transparency, validity, and reproducibility.

## Core Simulation Design

The initial simulation framework varies two key dimensions.

### 1. Selection Mechanism

- **Ignorable selection**: sample inclusion depends only on observed auxiliary variables \(X\).
- **Non-ignorable selection**: sample inclusion also depends directly on the outcome \(Y\).

The non-ignorable scenarios are included as diagnostic cases in which the ignorability assumption is violated. 
Standard adjustment methods based only on observed auxiliary variables are not expected to fully remove bias in these settings.

### 2. Quality of Auxiliary Information

- **Strong auxiliary variables**: observed auxiliary variables are strongly related to the outcome.
- **Weak auxiliary variables**: observed auxiliary variables only weakly explain the outcome, while unobserved factors play a larger role.

This design allows the project to examine not only which methods perform better, but also under what informational and assumption conditions adjustment becomes credible.

## Current Simulation Scenarios

The first script generates four scenario-specific finite populations and samples:

| Auxiliary information | Selection mechanism | Scenario              |
| --------------------- | ------------------- | --------------------- |
| Strong                | Ignorable           | `strong_ignorable`    |
| Strong                | Non-ignorable       | `strong_nonignorable` |
| Weak                  | Ignorable           | `weak_ignorable`      |
| Weak                  | Non-ignorable       | `weak_nonignorable`   |

For each scenario, the script creates:

- a finite population;
- a simple random probability sample;
- a biased non-probability sample;
- the true population mean of the outcome;
- naive sample estimates and initial bias summaries.

## Methodological Framework

The project will compare several families of adjustment methods for non-probability survey data.

| Method family                   | Main idea                                                                         | Core assumption                                                     | Data requirement                                                    | Role in this project         |
| ------------------------------- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------- | ---------------------------- |
| Naive estimator                 | Use the non-probability sample directly                                           | Sample is representative                                            | Non-probability sample with Y                                       | Baseline                     |
| Calibration / balancing weights | Directly adjust weights to match known auxiliary totals or margins                | Selection bias is captured by observed X                            | Population benchmarks or reference sample                           | Main comparison              |
| Propensity-score / IPW methods  | Estimate P(S=1 ∣ X), then weight by inverse propensity                            | Ignorability given X; adequate propensity model                     | Non-probability sample plus population or reference data on X       | Main comparison              |
| Mass imputation / prediction    | Estimate E(Y ∣ X), then predict for the target population or reference sample    | Ignorability given X; adequate outcome model                        | Non-probability sample with (Y, X), plus reference or population X  | Main comparison              |
| MRP                             | Multilevel outcome model plus poststratification over population cells            | Ignorability within modeled poststratification cells                | Population cell counts                                              | Later extension              |
| Doubly robust methods           | Combine propensity and outcome models                                             | Ignorability given X; one of the two models correctly specified     | Information for both selection and outcome modeling                 | Key methodological highlight |
| Non-ignorable selection models  | Explicitly model selection depending on Y or unobserved factors                   | Additional structural assumptions, proxies, or instruments          | Stronger assumptions or additional information                      | Later extension              |
| Bounds / partial identification | Estimate a plausible range rather than a single point                             | Weaker assumptions defining bounds                                  | Assumptions for bounding                                            | Later extension              |
| Sensitivity analysis            | Vary assumptions about unobserved selection and inspect robustness                | Diagnostic rather than full correction                              | Sensitivity parameters                                              | Later extension              |

## Current Status

- [x] Finite population simulation and sampling (`scripts/01_simulate_population.R`)
- [x] Baseline estimator comparison under full-population auxiliary information (`scripts/02_estimate_methods.R`)
  - probability sample mean;
  - naive non-probability sample mean;
  - calibration / balancing weights;
  - estimated IPW;
  - oracle IPW as a simulation benchmark (full-population setting only);
  - mass imputation;
  - doubly robust estimation.
- [x] Reference-sample integration using probability + non-probability samples (`scripts/03_reference_sample_integration.R`)
  - probability sample mean as design-based benchmark;
  - naive non-probability mean;
  - calibrated integration (targets estimated from p-sample design weights);
  - sample-membership IPW (inverse odds of np-sample membership);
  - mass imputation to p-sample;
  - doubly robust integration.
  - *Note:* oracle IPW is excluded here because the true selection probabilities are not available in the reference-sample setting.
- [x] Simulation results documented (`notes/Simulation_results_round1.md`, `notes/simulation_results_round2_reference_integration.md`)
- [ ] Hájek-normalized IPW and weight diagnostics
- [ ] Repeated simulation runs and performance evaluation
- [ ] Non-ignorable selection extensions and sensitivity analysis

Generated outputs (`data/processed/`, `outputs/`) are excluded from version control.

## Results Summary

**Round 1** (`02_estimate_methods.R`) used full-population auxiliary information as an ideal benchmark.
Observed-\(X\)-based adjustment methods work well when selection is ignorable and auxiliary variables are informative.
Under non-ignorable selection, these methods may reduce bias but do not fully remove it.

**Round 2** (`03_reference_sample_integration.R`) replaced full-population auxiliary information with a probability sample as reference.
The broad pattern is consistent with Round 1: adjustment works best when selection is ignorable and auxiliary variables are strong, partially reduces bias when strong auxiliary variables absorb some non-ignorable selection, and performs poorly when auxiliary variables are weak.
In the `strong_ignorable` scenario, adjusted bias dropped from 1.73 to 0.33–0.40 in Round 2, compared to 0.07 in Round 1 — the gap reflects sampling variability in the reference sample replacing exact population totals.
In the `weak_ignorable` scenario, adjustment performed slightly worse than the naive estimate in a single run, suggesting that adjustment is not automatically beneficial when auxiliary variables are weakly related to the outcome.

Both rounds confirm the core logic of the project: the success of adjustment depends jointly on the selection mechanism, the quality of auxiliary information, and the reference information available for adjustment.

## Planned Next Steps

1. Add stability and diagnostic measures for the implemented estimators:
   - Hájek-normalized IPW;
   - minimum and maximum weights;
   - weight distributions;
   - effective sample size.

2. Run repeated simulations to evaluate estimator performance:
   - mean bias;
   - absolute bias;
   - variance;
   - RMSE;
   - coverage where applicable.

3. Compare results across information settings:
   - full population auxiliary information — completed in Round 1;
   - probability reference sample only — completed in Round 2;
   - partial auxiliary information;
   - marginal population benchmarks only.

4. Examine model and method extensions:
   - bounded calibration weights;
   - raking / marginal calibration;
   - nonlinear outcome models;
   - richer auxiliary variables;
   - alternative propensity-score specifications.

5. Explore non-ignorable selection extensions:
   - sensitivity analysis;
   - selection models;
   - partial identification or bounds.

## Repository Structure

```text
nonprob-survey-inference/
├── README.md
├── nonprob-survey-inference.Rproj
├── scripts/
│   ├── 01_simulate_population.R
│   ├── 02_estimate_methods.R
│   └── 03_reference_sample_integration.R
├── notebooks/
├── notes/
│   ├── simulation_results_round1.md
│   └── simulation_results_round2_reference_integration.md
├── R/
├── data/
│   ├── raw/
│   └── processed/
├── outputs/
│   ├── figures/
│   └── tables/
└── references/
```

## Project Positioning

This project builds on survey methodology, statistical modeling, and data quality research. 
It treats non-probability sample adjustment not simply as a technical weighting problem, but as a question of assumptions, 
auxiliary information, and responsible statistical inference.

A central theme is that adjustment methods are only as credible as the data and assumptions that support them, 
and the project evaluates both where adjustment works and where it breaks down.
