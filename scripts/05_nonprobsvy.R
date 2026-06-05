############################################################
# 05_nonprobsvy_benchmark.R
#
# Purpose:
#   Benchmark the hand-coded reference-sample integration
#   estimators against the {nonprobsvy} package.
#
#   Uses the same population, samples, auxiliary variables, and
#   outcome as scripts/01--04.
#
# Output:
#   - outputs/tables/nonprobsvy_benchmark_single_rep.csv
#   - outputs/tables/nonprobsvy_package_results_single_rep.csv
#   - outputs/tables/nonprobsvy_weight_diagnostics_single_rep.csv
#   - outputs/tables/nonprobsvy_raw_fits_single_rep.rds
############################################################

library(survey)
library(nonprobsvy)

# -----------------------------
# 1. Reuse existing project code
# -----------------------------

source("scripts/01_simulate_population.R")

# Clear the single-run outputs that 01 created.
# This script manages its own fixed populations and samples.
rm(simulated_scenarios, scenario_summary)

source("scripts/03_reference_sample_integration.R")

# Clear the single-run outputs that 03 created.
# We only reuse its helper functions and integrate_one_scenario().
rm(integration_results)

# -----------------------------
# 2. Simulation settings
# -----------------------------

N_POP <- 100000
N_PROB <- 2000
N_NONPROB <- 3000  # expected; actual varies each draw

SCENARIOS <- list(
  strong_ignorable    = list(aux = "strong", sel = "ignorable"),
  strong_nonignorable = list(aux = "strong", sel = "nonignorable"),
  weak_ignorable      = list(aux = "weak",   sel = "ignorable"),
  weak_nonignorable   = list(aux = "weak",   sel = "nonignorable")
)

X_VARS <- c("age", "female", "higher_education", "urban")
SELECTION_FORMULA <- ~ age + female + higher_education + urban
TARGET_FORMULA    <- ~ y
OUTCOME_FORMULA   <- y ~ age + female + higher_education + urban

# -----------------------------
# 3. Generate fixed populations
# -----------------------------

set.seed(2026)

fixed_population <- list(
  strong = simulate_population(n_population = N_POP, aux_quality = "strong"),
  weak   = simulate_population(n_population = N_POP, aux_quality = "weak")
)


# -----------------------------
# 4. Helper functions
# -----------------------------

make_probability_design <- function(p_sample, population_size){
  p_sample$design_weight <- population_size / nrow(p_sample)
  svydesign(ids = ~1, weights = ~design_weight, data = p_sample)
}

extract_nonprobsvy_estimate <- function(fit, method_name, scenario_name, aux_quality,
                                        selection_mechanism, true_mean){
  extracted <- extract(fit)
  estimate <- as.numeric(extracted$mean[1])
  
  data.frame(
    scenario            = scenario_name,
    aux_quality         = aux_quality,
    selection_mechanism = selection_mechanism,
    method              = method_name,
    source              = "nonprobsvy",
    estimate            = estimate,
    true_mean           = true_mean,
    bias                = estimate - true_mean,
    absolute_bias       = abs(estimate - true_mean),
    SE = if("SE" %in% names(extracted)) as.numeric(extracted$SE[1]) else NA_real_,
    lower_bound = if ("lower_bound" %in% names(extracted)) as.numeric(extracted$lower_bound[1]) else NA_real_,
    upper_bound = if ("upper_bound" %in% names(extracted)) as.numeric(extracted$upper_bound[1]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

extract_weight_diagnostics <- function(fit, method_name, scenario_name) {
  w   <- weights(fit)
  q   <- as.numeric(quantile(w, probs = c(0.25, 0.50, 0.75)))
  ess <- sum(w)^2 / sum(w^2)
  
  data.frame(
    scenario      = scenario_name,
    method        = method_name,
    n_weights     = length(w),
    min_weight    = min(w),
    q25_weight    = q[1],
    median_weight = q[2],
    mean_weight   = mean(w),
    q75_weight    = q[3],
    max_weight    = max(w),
    sum_weight    = sum(w),
    ess           = ess,
    stringsAsFactors = FALSE
  )
}

format_hand_coded_results <- function(hand_results) {
  out <- hand_results[, c("scenario", "aux_quality", "selection_mechanism",
                          "method", "estimate", "true_mean",
                          "bias", "absolute_bias")]
  out$source      <- "hand_coded"
  out$SE          <- NA_real_
  out$lower_bound <- NA_real_
  out$upper_bound <- NA_real_
  out
}

# -----------------------------
# 5. Single-scenario benchmark
# -----------------------------

run_nonprobsvy_one_scenario <- function(scenario_name, base_population,
                                       selection_mechanism){
  population <- add_sampling_indicators(
    population                = base_population,
    selection_mechanism       = selection_mechanism,
    n_probability_sample      = N_PROB,
    expected_n_nonprob_sample = N_NONPROB
  )
  
  true_mean <- mean(population$y)
  aux_quality <- unique(population$aux_quality)
  selection_mechanism <- unique(population$selection_mechanism)
  
  p_sample <- population[population$probability_sample == 1, ]
  np_sample <- population[population$nonprob_sample == 1, ]
  
  p_design <- make_probability_design(p_sample, nrow(population))
  
  # Hand-coded benchmark
  hand_results <- format_hand_coded_results(
    integrate_one_scenario(scenario_name = scenario_name, population = population)
  )
  
  # nonprbsvy fit
  fit_ipw_mle <- nonprob(
    data = np_sample, selection = SELECTION_FORMULA,
    target = TARGET_FORMULA, svydesign = p_design,
    method_selection = "logit", se = FALSE
  )
  
  fit_ipw_gee <- nonprob(
    data = np_sample, selection = SELECTION_FORMULA,
    target = TARGET_FORMULA, svydesign = p_design,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
    se = FALSE
  )
  
  fit_mi_glm <- nonprob(
    data = np_sample, outcome = OUTCOME_FORMULA,
    svydesign = p_design,
    method_outcome = "glm", family_outcome = "gaussian", se = FALSE
  )
  
  fit_dr_glm <- nonprob(
    data = np_sample, selection = SELECTION_FORMULA,
    outcome = OUTCOME_FORMULA, svydesign = p_design,
    method_selection = "logit",
    method_outcome = "glm", family_outcome = "gaussian", se = FALSE
  )
  
  fits <- list(
    ipw_mle             = fit_ipw_mle,
    ipw_gee             = fit_ipw_gee,
    mass_imputation_glm = fit_mi_glm,
    doubly_robust_glm   = fit_dr_glm
  )
  
  method_names <- c(
    "nonprobsvy_ipw_mle_logit",
    "nonprobsvy_ipw_gee_logit_calibrated",
    "nonprobsvy_mass_imputation_glm",
    "nonprobsvy_doubly_robust_glm"
  )
  
  package_results <- do.call(rbind, Map(function(fit, name)
      extract_nonprobsvy_estimate(fit, name, scenario_name,
                                  aux_quality, selection_mechanism, true_mean),
    fits, method_names
  ))
  
  weight_diagnostics <- do.call(rbind, Map(
    function(fit, name)
      extract_weight_diagnostics(fit, name, scenario_name),
    fits[c("ipw_mle", "ipw_gee", "doubly_robust_glm")],
    method_names[c(1, 2, 4)]
  ))
  
  list(
    comparison         = rbind(hand_results, package_results),
    package_results    = package_results,
    weight_diagnostics = weight_diagnostics,
    raw_fits           = fits
  )
  
}


# -----------------------------
# 6. Run all scenarios
# -----------------------------

set.seed(42)

benchmark_results <- lapply(names(SCENARIOS), function(scenario_name) {
  cfg <- SCENARIOS[[scenario_name]]
  cat("Running:", scenario_name, "\n")
  run_nonprobsvy_one_scenario(
    scenario_name       = scenario_name,
    base_population     = fixed_population[[cfg$aux]],
    selection_mechanism = cfg$sel
  )
})
names(benchmark_results) <- names(SCENARIOS)

comparison_results <- do.call(rbind, lapply(benchmark_results, `[[`, "comparison"))
package_results    <- do.call(rbind, lapply(benchmark_results, `[[`, "package_results"))
weight_diagnostics <- do.call(rbind, lapply(benchmark_results, `[[`, "weight_diagnostics"))

# -----------------------------
# 7. Save outputs
# -----------------------------

write.csv(comparison_results,
          "outputs/tables/nonprobsvy_benchmark_single_rep.csv",
          row.names = FALSE)

write.csv(package_results,
          "outputs/tables/nonprobsvy_package_results_single_rep.csv",
          row.names = FALSE)

write.csv(weight_diagnostics,
          "outputs/tables/nonprobsvy_weight_diagnostics_single_rep.csv",
          row.names = FALSE)

saveRDS(benchmark_results,
        "outputs/tables/nonprobsvy_raw_fits_single_rep.rds")

cat("\nDone. Package benchmark summary:\n")
print(package_results[, c("scenario", "method", "estimate",
                          "true_mean", "bias", "absolute_bias")],
      row.names = FALSE)

cat("\nWeight diagnostics:\n")
print(weight_diagnostics, row.names = FALSE)























