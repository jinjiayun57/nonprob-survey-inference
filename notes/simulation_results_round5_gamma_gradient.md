# Simulation Results — Round 5: Non-Ignorability Gradient

**Script:** `07_nonignorability_gradient.R`
**Plot:** `08_plot_gamma_results.R`
**Output:** `outputs/tables/gamma_performance.csv`, `outputs/tables/gamma_diagnostics.csv`, `outputs/figures/gamma_bias_plot.png`

---

## Purpose

Rounds 1–4 treated the selection mechanism as a binary switch: ignorable (selection depends only on observed X) vs. non-ignorable (selection additionally depends on Y, via a single fixed coefficient of 0.60 on standardized Y).

This round replaces that switch with a continuous non-ignorability parameter, gamma, in the selection model:

```
linear_predictor = base_predictor(X) + gamma * z_y
```

`gamma = 0` reproduces the ignorable scenarios; `gamma = 0.60` reproduces the non-ignorable scenarios from rounds 1–4. The grid extends to `gamma = 1.2` to trace out how performance changes as selection becomes progressively more outcome-dependent, rather than comparing only two points.

---

## Design

- Gamma grid: 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2
- Aux quality: strong, weak
- `N_SIM = 200` replicates per (aux_quality x gamma) cell
- Same fixed-population seed (2026) and replicate seed (42) as previous rounds
- Estimators: naive non-probability mean, calibration, mass imputation, membership IPW, doubly robust (logistic propensity + linear outcome model, as in rounds 1–4)
- Probability sample mean included as a reference (SRS, unaffected by gamma)

---

## Diagnostic: realized selection–Y association

`gamma` is a model coefficient, not a directly interpretable quantity. `gamma_diagnostics.csv` reports the realized correlation between Y and the non-probability selection probability, holding the overall selection rate fixed at ~3%.

| gamma | strong | weak |
|---:|---:|---:|
| 0.0 | 0.267 | 0.064 |
| 0.2 | 0.482 | 0.317 |
| 0.4 | 0.618 | 0.496 |
| 0.6 | 0.686 | 0.595 |
| 0.8 | 0.711 | 0.639 |
| 1.0 | 0.709 | 0.651 |
| 1.2 | 0.695 | 0.648 |

The correlation rises with gamma but plateaus (and slightly declines) beyond gamma ≈ 0.8 in both regimes. This is a logistic-saturation effect: with the selection rate pinned at ~3%, larger gamma sharpens the selection probability into more of a step function over a thin high-Y tail, rather than spreading the Y-dependence further. Beyond gamma ≈ 0.8, larger gamma is not meaningfully "more non-ignorable" in realized correlation terms — but as the bias results below show, it still translates into growing estimation bias.

---

## Results: Mean Bias

**Strong auxiliary information**

| Method | γ=0 | γ=0.2 | γ=0.4 | γ=0.6 | γ=0.8 | γ=1.0 | γ=1.2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| prob_sample_mean | 0.004 | -0.008 | -0.051 | 0.013 | 0.003 | 0.006 | 0.012 |
| np_naive_mean | 1.874 | 3.758 | 5.601 | 7.363 | 9.039 | 10.617 | 12.027 |
| calibration | -0.008 | 0.710 | 1.369 | 2.134 | 2.862 | 3.602 | 4.328 |
| membership_ipw | 0.010 | 0.727 | 1.395 | 2.183 | 2.907 | 3.653 | 4.433 |
| doubly_robust | 0.003 | 0.715 | 1.370 | 2.130 | 2.835 | 3.535 | 4.208 |

**Weak auxiliary information**

| Method | γ=0 | γ=0.2 | γ=0.4 | γ=0.6 | γ=0.8 | γ=1.0 | γ=1.2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| prob_sample_mean | 0.001 | 0.008 | 0.006 | -0.013 | 0.016 | -0.034 | -0.027 |
| np_naive_mean | 0.367 | 2.035 | 3.621 | 5.218 | 6.719 | 8.136 | 9.432 |
| calibration | -0.008 | 1.642 | 3.207 | 4.813 | 6.321 | 7.784 | 9.163 |
| membership_ipw | -0.009 | 1.639 | 3.200 | 4.807 | 6.307 | 7.768 | 9.141 |
| doubly_robust | -0.010 | 1.638 | 3.199 | 4.804 | 6.306 | 7.767 | 9.139 |

(Mass imputation is omitted — numerically identical to calibration in every cell, as in rounds 1–4.) Variance is essentially flat across gamma for every method (0.02–0.08), so RMSE tracks |mean bias| almost exactly; the gradient is a bias story, not a variance story.

---

## Main Findings

**Differences among the X-based adjustment methods are small once gamma > 0.**
At gamma = 0, all adjustment methods are approximately unbiased, as expected. 
For gamma > 0, calibration, membership IPW, and DR track each other closely and grow together. 
In the strong-auxiliary regime, they remain far below the naive non-probability mean, but their differences from one another are small relative to the bias introduced by outcome-dependent selection. 
This suggests that, in this design, estimator choice matters less than whether the available auxiliary variables can account for the outcome-related component of selection.

**Doubly robust estimation does not protect against non-ignorable selection.**
DR performs very similarly to calibration and IPW across the gamma gradient. 
This is consistent with the previous ML-propensity comparison: DR can reduce sensitivity to propensity-model specification, 
but it does not address violations of the ignorability assumption itself. 
Once selection depends on Y through a component not captured by X, outcome and propensity models built only on X cannot remove that bias.

**Residual bias grows roughly linearly in gamma, even after the selection–Y correlation plateaus.**
Per-step bias increments for calibration are roughly constant (~0.66–0.77 per 0.2 gamma) for strong, 
and decline only mildly (1.65 → 1.38) for weak — both continue increasing well past gamma = 0.8, where the correlation diagnostic above has flattened. 
Selection becomes more "threshold-like" (extreme-Y units selected with near-certainty) rather than more correlated in the linear sense, 
but this still shifts the selected mean further from the population mean.

**Auxiliary quality determines how much of the Y-driven bias adjustment can absorb.**
At the same gamma, adjusted bias is much smaller under strong auxiliary information than under weak auxiliary information. 
For example, at gamma = 0.6, calibration bias is 2.13 in the strong regime but 4.81 in the weak regime. Adjustment methods work by proxying selection through observed X; 
when X is strongly related to Y, part of the outcome-driven selection bias can be absorbed. When X is a poor proxy for Y, 
that channel is largely closed, and the adjusted curves move close to the naive non-probability curve.

**Probability sample mean remains unbiased throughout**, confirming that gamma only affects non-probability selection, 
not the SRS reference sample — a useful internal consistency check.

---

## Takeaway

The gamma gradient shows that the main boundary condition for non-probability adjustment is not which X-based estimator is used, 
but whether the available auxiliary variables can account for the outcome-related part of selection.
Calibration, IPW, and doubly robust estimation behave very similarly once selection becomes outcome-dependent. 
Strong auxiliary information absorbs part of the bias, while weak auxiliary information offers little protection. 
This motivates the next extension: methods or sensitivity analyses that explicitly address the unobserved outcome-dependent component of selection, 
rather than relying only on adjustment by observed X.

