# Simulation Results — Round 2: Reference-Sample Integration

**Date:** 2026-06-03  
**Script:** `03_reference_sample_integration.R`  
**Output:** `outputs/tables/integration_comparison.csv`

---

## Purpose

Round 2 extends the analysis from `02_estimate_methods.R`.

In `02`, adjustment methods used full population auxiliary information. This is an idealized setting because the true population distribution of \(X\) is known.

In `03`, full population \(X\) is no longer used directly. Instead, the probability sample is used as a reference sample. This makes the setup more realistic:

> The probability sample provides representativeness; the non-probability sample provides additional outcome information.

The goal is to evaluate whether probability + non-probability sample integration can reduce bias when only reference-sample information is available.

---

## Main Difference from Round 1

| Round | Reference information | Main question |
|---|---|---|
| Round 1 / `02` | Full population auxiliary totals | Can adjustment work under ideal population-X information? |
| Round 2 / `03` | Probability sample as reference | Can adjustment still work when population-X information is estimated from a probability sample? |

---

## Results: Absolute Bias

| Method | strong_ignorable | strong_nonignorable | weak_ignorable | weak_nonignorable |
|---|---:|---:|---:|---:|
| probability_sample_mean | 0.54 | 0.02 | 0.03 | 0.02 |
| np_naive_mean | 1.73 | 7.43 | 0.22 | 5.16 |
| calibrated_integration | 0.33 | 2.21 | 0.32 | 4.83 |
| membership_ipw | 0.40 | 2.12 | 0.38 | 4.83 |
| mass_imputation_integration | 0.33 | 2.21 | 0.32 | 4.83 |
| doubly_robust_integration | 0.39 | 2.13 | 0.36 | 4.81 |

*Single simulation run; results are diagnostic rather than final.*

---

## Main Findings

Compared with Round 1, the same broad pattern remains:

- Adjustment works best when selection is ignorable and auxiliary variables are strong.
- Non-ignorable selection leaves residual bias even after adjustment.
- Weak auxiliary variables limit the value of adjustment.
- Calibration and mass imputation remain numerically identical in the current linear setup.

The main new finding is that reference-sample integration is less ideal than full-population-X adjustment. In the `strong_ignorable` scenario, adjusted bias is reduced from 1.73 to about 0.33–0.40, but not as close to zero as in Round 1. This is expected because the probability sample is only an estimate of the population structure and contains sampling variability.

In the `strong_nonignorable` scenario, adjustment still reduces bias substantially, from 7.43 to about 2.1–2.2. This suggests that strong auxiliary variables can partially absorb selection bias, even when selection depends directly on \(Y\). However, observed-\(X\) methods cannot fully remove the bias.

In the `weak_ignorable` scenario, the naive bias is already small, and adjustment performs slightly worse in this single run. This suggests that adjustment is not automatically beneficial when auxiliary variables are weakly related to the outcome.

In the `weak_nonignorable` scenario, all adjusted methods perform poorly. This remains the hardest case: selection depends directly on \(Y\), while observed \(X\) provides little information about \(Y\).

---

## Takeaway

Round 2 confirms the core logic of the project:

> Integration methods can improve non-probability sample estimates when the reference sample is informative and auxiliary variables are strong, but they cannot overcome weak auxiliary information or non-ignorable selection.

The next step is Monte Carlo simulation to assess whether these single-run patterns hold across repeated samples.
