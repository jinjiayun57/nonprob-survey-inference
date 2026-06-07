# Simulation Results — Round 3: ML Propensity Models

**Script:** `06_ml_propensity.R`  
**Output:** `outputs/tables/ml_propensity_performance.csv`

---

## Purpose

This round extends the Monte Carlo framework by replacing the logistic propensity model with random forest (`ranger`) and gradient boosting (`xgboost`). The outcome model is held fixed as a linear model to isolate the effect of propensity model choice. Six estimators are compared: logistic/RF/XGB × IPW/DR.

The simulation uses the same population seed (2026), replicate seed (42), sample sizes, and 2×2 scenario design as the previous Monte Carlo evaluation (`N_SIM = 500`).

---

## Design

The simulation compares three propensity models:

- logistic regression;
- random forest (`ranger`);
- gradient boosting (`xgboost`).

For each propensity model, two estimators are evaluated:

- IPW;
- doubly robust estimation.

The same four simulation scenarios are used as in previous rounds:

- `strong_ignorable`;
- `strong_nonignorable`;
- `weak_ignorable`;
- `weak_nonignorable`.

Performance is summarized using mean bias, MAE, RMSE, and effective sample size.

---

## Results: Mean Bias

| Method | strong_ignorable | strong_nonignorable | weak_ignorable | weak_nonignorable |
|---|---:|---:|---:|---:|
| prob_sample_mean | 0.010 | 0.000 | 0.002 | 0.009 |
| np_naive_mean | 1.871 | 7.443 | 0.394 | 5.237 |
| logistic_ipw | 0.002 | 2.163 | 0.019 | 4.837 |
| rf_ipw | -0.127 | 2.803 | 0.003 | 4.832 |
| xgb_ipw | 0.035 | 2.299 | 0.025 | 4.845 |
| logistic_dr | 0.003 | 2.123 | 0.018 | 4.836 |
| rf_dr | 0.003 | 2.121 | 0.016 | 4.835 |
| xgb_dr | 0.001 | 2.121 | 0.017 | 4.837 |

---

## Main Findings

**DR is stable across propensity models.**  
Across all scenarios, `logistic_dr`, `rf_dr`, and `xgb_dr` produce nearly identical mean bias. Once the linear outcome model is included, switching the propensity model has little impact on the final DR estimate.

**IPW is sensitive to propensity model choice.**  
In `strong_ignorable`, logistic IPW performs best. This is expected because the current ignorable selection mechanism is generated from a relatively simple function of the observed auxiliary variables. RF IPW has higher bias in this scenario, despite a higher effective sample size.

**ML propensity models may help under model uncertainty, but not uniformly.**  
In `weak_ignorable`, RF IPW performs slightly better than logistic IPW. This suggests that flexible propensity modeling can sometimes improve the weighting estimator, but the gain is scenario-dependent and should be evaluated together with weight stability and RMSE.

**Non-ignorable selection remains the binding constraint.**  
In both non-ignorable scenarios, all propensity-based estimators remain substantially biased. Flexible ML propensity models cannot remove bias caused by direct outcome-dependent selection when the reference sample does not observe the outcome.

---

## Takeaway

The results suggest a nuanced role for machine-learning propensity models. When the parametric propensity model is close to correctly specified, logistic weighting remains highly competitive. When propensity specification is uncertain, ML-based weighting can sometimes improve IPW, but the improvement is not automatic.

For DR estimators, propensity model choice matters much less in the current design. The main practical implication is that DR estimation provides a more stable strategy when propensity model specification is uncertain, while ML propensity models are most useful for improving the weighting component itself.