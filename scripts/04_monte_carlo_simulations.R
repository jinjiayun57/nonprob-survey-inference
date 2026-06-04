############################################################
# 04_monte_carlo_simulation.R
#
# Purpose:
#   Monte Carlo evaluation of reference-sample integration
#   estimators across the 2x2 simulation scenarios.
#
#   Population is fixed (generated once per scenario).
#   Each replicate draws new probability and non-probability
#   samples from the fixed population.
#
# Output:
#   - outputs/tables/mc_raw_results.csv    (all replicates)
#   - outputs/tables/mc_performance.csv    (summary metrics)
############################################################

source("scripts/01_simulate_population.R")

# Clear the single-run outputs that 01 created
# this script manages its own populations
rm(simulated_scenarios, scenario_summary)

# -----------------------------
# 0. Reuse helper functions from 03
# -----------------------------

source("scripts/03_reference_sample_integration.R")

# Clear the single-run outputs that 03 created.
rm(integration_results)

# -----------------------------
# 1. Simulation settings
# -----------------------------

N_SIM <- 500
N_POP <- 100000
N_PROB <- 2000
N_NONPROB <- 3000 #expected; actual varies each draw

SCENARIOS <- list(
  strong_ignorable = list(aux = "strong", sel = "ignorable"),
  strong_nonignorable = list(aux = "strong", sel = "nonignorable"),
  weak_ignorable      = list(aux = "weak",   sel = "ignorable"),
  weak_nonignorable   = list(aux = "weak",   sel = "nonignorable")
)

# -----------------------------
# 2. Generate fixed populations (one per aux_quality)
# -----------------------------

set.seed(2026)

fixed_population <- list(
  strong = simulate_population(n_population = N_POP, aux_quality = "strong"),
  weak = simulate_population(n_population = N_POP, aux_quality = "weak")
)

# -----------------------------
# 3. Single-replicate function
# -----------------------------
#
# Draws new p-sample and np-sample from the fixed population,
# runs all estimators, returns a one-row-per-method data frame.

run_one_replicate <- function(rep_id, base_population, selection_mechanism, scenario_name){
  
  population <- add_sampling_indicators(
    population = base_population,
    selection_mechanism = selection_mechanism,
    n_probability_sample = N_PROB,
    expected_n_nonprob_sample = N_NONPROB
  )
  
  result <- integrate_one_scenario(
    scenario_name = scenario_name,
    population = population
  )
  
  result$rep_id <- rep_id
  return(result)
}

# -----------------------------
# 4. Run Monte Carlo loop
# -----------------------------

set.seed(42)

mc_raw <- do.call(rbind, lapply(names(SCENARIOS), function(scenario_name){
  cfg <- SCENARIOS[[scenario_name]]
  base_pop <- fixed_population[[cfg$aux]]
  
  cat("Running scenario:", scenario_name, "\n")
  
  do.call(rbind, lapply(seq_len(N_SIM), function(rep_id){
    run_one_replicate(
      rep_id = rep_id,
      base_population = base_pop,
      selection_mechanism = cfg$sel,
      scenario_name = scenario_name 
    )
  }))
}))

# -----------------------------
# 5. Compute performance metrics
# -----------------------------
#
# For each scenario x method:
#   - mean_bias    : mean of (estimate - true_mean)
#   - abs_bias     : mean of |estimate - true_mean|
#   - variance     : var of estimates across replicates
#   - rmse         : sqrt(mean((estimate - true_mean)^2))

# Exclude the "population_truth" row -- it is constant and not an estimator.
mc_estimators <- mc_raw[mc_raw$method != "population_truth", ]

mc_performance <- do.call(rbind, lapply(
  split(mc_estimators, list(mc_estimators$scenario, mc_estimators$method)),
  function(df){
    errors <- df$estimate - df$true_mean
    
    data.frame(
      scenario = df$scenario[1],
      aux_quality         = df$aux_quality[1],
      selection_mechanism = df$selection_mechanism[1],
      method              = df$method[1],
      n_replicates        = nrow(df),
      true_mean           = df$true_mean[1],
      mean_estimate       = mean(df$estimate, na.rm = TRUE),
      mean_bias           = mean(errors,        na.rm = TRUE),
      mae                 = mean(abs(errors),   na.rm = TRUE),
      variance            = var(df$estimate,    na.rm = TRUE),
      rmse                = sqrt(mean(errors^2, na.rm = TRUE)),
      stringsAsFactors    = FALSE
    )
  }
))

# sort for readability

mc_performance <- mc_performance[
  order(mc_performance$scenario, mc_performance$method), ]


# -----------------------------
# 6. Save outputs
# -----------------------------

write.csv(
  mc_raw,
  file      = "outputs/tables/mc_raw_results.csv",
  row.names = FALSE
)

write.csv(
  mc_performance,
  file      = "outputs/tables/mc_performance.csv",
  row.names = FALSE
)

cat("\nDone. Performance summary:\n")
print(mc_performance[, c("scenario", "method", "mean_bias", "mae",
                         "variance", "rmse")],
      row.names = FALSE)












