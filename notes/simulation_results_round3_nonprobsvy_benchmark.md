# Simulation Results — Round 3: nonprobsvy Benchmark

**Script:** `05_nonprobsvy.R`  
**Output:** `outputs/tables/nonprobsvy_benchmark_single_rep.csv`

---

## Purpose

This round benchmarks the hand-coded reference-sample integration estimators against the `{nonprobsvy}` R package.
The goal is not to introduce a new simulation design, but to check whether the project's manually implemented estimators align with a dedicated external package for non-probability survey inference.

The benchmark uses the same four scenarios as previous rounds:

- `strong_ignorable`;
- `strong_nonignorable`;
- `weak_ignorable`;
- `weak_nonignorable`.

This is a single-replicate benchmark rather than a full Monte Carlo evaluation.

---

## Methods Compared

The `{nonprobsvy}` package is used to estimate:

- IPW via MLE logit;
- calibrated IPW via GEE logit;
- mass imputation GLM;
- doubly robust GLM.

These estimates are compared with the hand-coded estimators from `03_reference_sample_integration.R`.

---

## Main Findings

The package benchmark confirms that the hand-coded mass imputation and doubly robust estimators are closely aligned with `{nonprobsvy}` output.
Mass imputation and DR estimates match to numerical precision or near numerical precision, suggesting that the project implementation is consistent with the package-based reference-sample prediction framework.

The calibrated IPW estimator from `{nonprobsvy}` is close to the hand-coded calibrated / membership-weighting results, although the estimators are not identical because they use different calibration and estimating-equation constructions.

The broad substantive pattern remains unchanged: adjustment works well under ignorable selection, partially reduces bias under strong non-ignorable selection, and performs poorly when auxiliary variables are weak and selection is non-ignorable.

---

## Takeaway

This round validates the hand-coded estimators against an external package implementation. It supports the reliability of the existing simulation workflow before extending the project to machine-learning propensity models and non-ignorable selection methods.
