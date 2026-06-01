# Simulation Results — Round 1

**Scripts:** `01_simulate_population.R` → `02_estimate_methods.R`  
**Output:** `outputs/tables/estimator_comparison.csv`

---

## Setup

- Population size: N = 100,000
- Probability sample size: n ≈ 2,000
- Non-probability sample size: n ≈ 3,000–3,100
- Auxiliary variables: age, female, higher_education, urban
- Outcome: y (continuous)

Four scenarios crossing two dimensions:

| Scenario | Aux quality | Selection mechanism |
|---|---|---|
| strong_ignorable | Strong | Ignorable |
| strong_nonignorable | Strong | Non-ignorable |
| weak_ignorable | Weak | Ignorable |
| weak_nonignorable | Weak | Non-ignorable |

---

## Results: Absolute Bias by Method and Scenario

| Method | strong_ignorable | strong_nonignorable | weak_ignorable | weak_nonignorable |
|---|---|---|---|---|
| probability_sample_mean | 0.54 | 0.02 | 0.03 | 0.02 |
| nonprob_naive_mean | 1.73 | 7.43 | 0.22 | 5.16 |
| calibration_mean | **0.07** | 2.16 | 0.25 | 4.87 |
| estimated_ipw_mean | **0.09** | 2.24 | **0.07** | 4.60 |
| oracle_ipw_mean | 0.74 | **0.51** | 0.58 | **0.31** |
| mass_imputation_mean | **0.07** | 2.16 | 0.25 | 4.87 |
| doubly_robust_mean | 0.14 | 2.08 | 0.29 | 4.85 |

*Note: results from a single simulation run; variance not yet assessed.*

---

## Key Findings

### 1. Strong auxiliary variables + ignorable selection
In the `strong_ignorable` scenario, the naive non-probability sample estimate is biased, 
but all observed-X adjustment methods substantially reduce the bias.

This is the ideal setting for standard adjustment methods:

 - selection depends only on observed \(X\);
 - \(X\) is strongly related to \(Y\);
 - therefore, adjusting for \(X\) removes most of the selection-induced bias.

This scenario supports the basic expectation that calibration, IPW, mass imputation, 
and doubly robust estimation can work well when the ignorability assumption is plausible 
and auxiliary variables are informative.

### 2. Strong auxiliary variables + Non-ignorable selection
In the `strong_nonignorable` scenario, the naive non-probability sample estimate has a much larger bias.

Observed-X adjustment methods reduce the bias substantially, but do not fully remove it.
This is expected because the true selection mechanism depends directly on \(Y\), while the adjustment methods only use observed \(X\).

The result illustrates an important point:

> Observed-X adjustment methods can partially reduce bias under non-ignorable selection when \(X\)
is strongly related to \(Y\), but residual bias remains because the ignorability assumption is violated.

The oracle IPW estimator performs better because it uses the true simulated selection probabilities,
including the part of selection that depends on \(Y\). This is useful as a simulation benchmark, 
but it is not generally available in real applications. 

### 3. Weak auxiliary variables + ignorable selection

In the `weak_ignorable` scenario, the naive non-probability bias is already small. 
This is because selection depends on \(X\), but \(X\) is only weakly related to \(Y\).
Therefore, sample imbalance in \(X\) does not strongly translate into bias in \(Y\).

Adjustment methods do not show large improvements in this scenario. 
Some methods may even introduce small additional variation or bias in a single run.

This highlights that adjustment is most useful when auxiliary variables are related both to selection and to the outcome.

### 4. Weak auxiliary variables + non-ignorable selection

The `weak_nonignorable` scenario is the most difficult case.

The naive non-probability sample estimate is strongly biased, and observed-X adjustment
methods only reduce the bias slightly. This is expected because:

 - selection depends directly on \(Y\);
 - observed \(X\) is weakly related to \(Y\);
 - therefore, the available auxiliary variables do not capture the main source of selection bias.
 
The scenario shows the limits of standard observed-X adjustment methods.
 

### Interpretation

The first results support the central design logic of the project:

> The success of non-probability sample adjustment depends jointly on the selection mechanism and the quality of auxiliary information.

The results suggest that the key question is not simply which estimator performs best, 
but under what assumptions and information conditions each estimator is credible.

Observed-X adjustment methods perform well when selection is ignorable and auxiliary variables
are informative. When selection is non-ignorable, these methods may reduce bias but are not expected
to fully correct it. When auxiliary variables are weak, even sophisticated methods have limited information
with which to adjust the sample.

### Notes on individual methods

#### Calibration and mass imputation

In the current implementation, calibration and mass imputation give identical or nearly identical results. 
This is not necessarily an error.

Under linear calibration and linear outcome regression using the same auxiliary variables, 
calibration and regression prediction can be closely related or algebraically equivalent. 
This is an important methodological point to document.

Future versions may examine whether this equivalence changes when using:

 - bounded calibration weights;
 - raking;
 - nonlinear outcome models;
 - richer auxiliary variables;
 - package-based calibration using survey::calibrate().

#### Doubly robust estimator

The doubly robust estimator combines an outcome model and a propensity model. 
In this implementation, it can be interpreted as:

\[
\text{DR} = \text{outcome prediction} + \text{IPW-weighted residual correction}.
\]

In the current single-run results, the doubly robust estimator performs similarly to calibration and mass imputation.
It does not fully correct bias under non-ignorable selection, because its double-robust property still relies on ignorability given observed \(X\). 
When selection depends directly on \(Y\), observed-\(X\)-based propensity and outcome models are not sufficient to remove all bias.

---

## Working Interpretation from Round 1

| Condition | Tentative interpretation  |
|---|---|
| Ignorable selection, strong aux | Calibration, estimated IPW, mass imputation, and DR all work well in this run |
| Ignorable selection, weak aux | Gains are small because naive bias is already limited; estimated IPW performs well in this run |
| Non-ignorable selection | Observed-\(X\) methods are insufficient; additional information, sensitivity analysis, or non-ignorable models are needed |
| Benchmark / simulation only | Oracle IPW helps diagnose what is lost when the true selection mechanism is unknown |

---

## Limitations and Next Steps

- **Single run only** — results may be sensitive to simulation randomness; Monte Carlo (≥100 replications) needed to assess variance and RMSE
- **Linear models throughout** — calibration, outcome model, and propensity model are all linear/logistic; non-linear outcome relationships not yet explored
- **Binary non-ignorability** — the current setup treats ignorable vs non-ignorable as a binary switch; gradations of non-ignorability would be more realistic

---
