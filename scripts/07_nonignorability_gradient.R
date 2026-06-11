############################################################
# 07_nonignorability_gradient.R
#
# Purpose:
#   Previous rounds (01-06) treated the selection mechanism as a
#   binary switch: "ignorable" (selection depends only on observed X)
#   vs. "nonignorable" (selection additionally depends on Y, with a
#   single fixed coefficient).
#
#   This script replaces that switch with a continuous
#   non-ignorability parameter, gamma, controlling how strongly
#   selection depends on Y:
#
#     linear_predictor = base_predictor(X) + gamma * z_y
#
#     gamma = 0    -> ignorable selection (as in *_ignorable scenarios)
#     gamma = 0.60 -> matches the *_nonignorable scenarios used in 01-06
#     gamma > 0.60 -> selection more strongly driven by Y than before
#
#   For each auxiliary-quality regime (strong / weak) and each gamma
#   in a grid, we run N_SIM replicates of the same reference-sample
#   integration estimators used in 03/04, and summarize bias and
#   RMSE as a function of gamma.
#
#   Goal: trace out *how* (gradually vs. abruptly) each estimator's
#   performance degrades as selection becomes less ignorable, rather
#   than comparing only two discrete points.
#
# Output:
#   - outputs/tables/gamma_raw_results.csv     (all replicates)
#   - outputs/tables/gamma_performance.csv     (summary metrics)
#   - outputs/tables/gamma_diagnostics.csv     (realized selection-Y
#                                                association per gamma)
############################################################

# -----------------------------
# 0. Reuse helper / estimator functions from 01 and 03
# -----------------------------

source("scripts/01_simulate_population.R")

# Clear single-run outputs that 01 created; this script manages its
# own populations and sampling indicators.
rm(simulated_scenarios, scenario_summary)

source("scripts/03_reference_sample_integration.R")

# Clear single-run outputs that 03 created (built on the binary design).
rm(simulated_scenarios, integration_results)


# -----------------------------
# 1. Continuous-gamma sampling indicator function
# -----------------------------
#
# Generalizes add_sampling_indicators() from 01_simulate_population.R.
# The auxiliary part of the selection model (the dependence on
# observed X) is held fixed at the same coefficients used throughout
# the project. Only the Y-dependence is varied, via gamma.

add_sampling_indicators_gamma <- function(
    population,
    gamma = 0,
    n_probability_sample = 2000,
    expected_n_nonprob_sample = 3000) {

  n_population <- nrow(population)

  # Probability sample: simple random sample, as in 01.
  probability_sample <- rep(0, n_population)
  probability_sample[sample(seq_len(n_population), size = n_probability_sample)] <- 1

  z_age   <- standardize(population$age)
  z_edu   <- standardize(population$higher_education)
  z_urban <- standardize(population$urban)
  z_female <- standardize(population$female)
  z_y     <- standardize(population$y)

  # Same X-dependent base predictor as in 01_simulate_population.R.
  base_predictor <- 0.50 * z_edu + 0.35 * z_urban - 0.25 * z_age + 0.15 * z_female

  linear_predictor <- base_predictor + gamma * z_y

  target_rate <- expected_n_nonprob_sample / n_population
  intercept <- calibrate_intercept(linear_predictor, target_rate)

  nonprob_selection_probability <- plogis(intercept + linear_predictor)
  nonprob_sample <- draw_bernoulli_sample(nonprob_selection_probability)

  population$selection_mechanism <- if (gamma == 0) "ignorable" else paste0("gamma_", gamma)
  population$gamma <- gamma
  population$probability_sample <- probability_sample
  population$nonprob_sample <- nonprob_sample
  population$nonprob_selection_probability <- nonprob_selection_probability

  return(population)
}


# -----------------------------
# 2. Simulation settings
# -----------------------------

N_SIM <- 200          # replicates per (aux_quality x gamma) cell
N_POP <- 100000
N_PROB <- 2000
N_NONPROB <- 3000     # expected; actual varies each draw

# Gamma grid: 0 = ignorable (as in *_ignorable scenarios),
# 0.60 = matches *_nonignorable scenarios from 01-06.
# Extending beyond 0.60 shows what happens under selection that is
# even more strongly outcome-driven than previously considered.
GAMMA_GRID <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2)

AUX_QUALITIES <- c("strong", "weak")


# -----------------------------
# 3. Fixed populations (one per aux_quality)
# -----------------------------
#
# Same population-generation seed as 02-06, so results are
# comparable to earlier rounds.

set.seed(2026)

fixed_population <- list(
  strong = simulate_population(n_population = N_POP, aux_quality = "strong"),
  weak   = simulate_population(n_population = N_POP, aux_quality = "weak")
)


# -----------------------------
# 4. Diagnostics: realized selection-Y association per gamma
# -----------------------------
#
# nonprob_selection_probability is a deterministic function of
# (population, gamma); only the Bernoulli draw is random. We use this
# to characterize how strongly selection actually depends on Y at
# each gamma, separately for the strong- and weak-auxiliary
# populations (the relationship between z_y and the outcome differs
# between the two).

gamma_diagnostics <- do.call(rbind, lapply(AUX_QUALITIES, function(aux) {
  base_pop <- fixed_population[[aux]]

  do.call(rbind, lapply(GAMMA_GRID, function(g) {
    pop_g <- add_sampling_indicators_gamma(
      population = base_pop,
      gamma = g,
      n_probability_sample = N_PROB,
      expected_n_nonprob_sample = N_NONPROB
    )

    data.frame(
      aux_quality = aux,
      gamma = g,
      cor_y_selection_prob = cor(pop_g$y, pop_g$nonprob_selection_probability),
      mean_selection_prob  = mean(pop_g$nonprob_selection_probability)
    )
  }))
}))

write.csv(
  gamma_diagnostics,
  file = "outputs/tables/gamma_diagnostics.csv",
  row.names = FALSE
)


# -----------------------------
# 5. Single-replicate function
# -----------------------------

run_one_replicate_gamma <- function(rep_id, base_population, gamma, scenario_name) {

  population <- add_sampling_indicators_gamma(
    population = base_population,
    gamma = gamma,
    n_probability_sample = N_PROB,
    expected_n_nonprob_sample = N_NONPROB
  )

  result <- integrate_one_scenario(
    scenario_name = scenario_name,
    population = population
  )

  result$gamma <- gamma
  result$rep_id <- rep_id
  return(result)
}


# -----------------------------
# 6. Run Monte Carlo loop over aux_quality x gamma
# -----------------------------

set.seed(42)

gamma_raw <- do.call(rbind, lapply(AUX_QUALITIES, function(aux) {
  base_pop <- fixed_population[[aux]]

  do.call(rbind, lapply(GAMMA_GRID, function(g) {

    scenario_name <- paste0(aux, "_gamma_", g)
    cat("Running scenario:", scenario_name, "\n")

    do.call(rbind, lapply(seq_len(N_SIM), function(rep_id) {
      run_one_replicate_gamma(
        rep_id = rep_id,
        base_population = base_pop,
        gamma = g,
        scenario_name = scenario_name
      )
    }))
  }))
}))


# -----------------------------
# 7. Performance metrics by aux_quality x gamma x method
# -----------------------------

gamma_estimators <- gamma_raw[gamma_raw$method != "population_truth", ]

gamma_performance <- do.call(rbind, lapply(
  split(gamma_estimators, list(gamma_estimators$aux_quality,
                                gamma_estimators$gamma,
                                gamma_estimators$method)),
  function(df) {
    if (nrow(df) == 0) return(NULL)

    errors <- df$estimate - df$true_mean

    data.frame(
      aux_quality   = df$aux_quality[1],
      gamma         = df$gamma[1],
      method        = df$method[1],
      n_replicates  = nrow(df),
      true_mean     = df$true_mean[1],
      mean_estimate = mean(df$estimate, na.rm = TRUE),
      mean_bias     = mean(errors,      na.rm = TRUE),
      mae           = mean(abs(errors), na.rm = TRUE),
      variance      = var(df$estimate,  na.rm = TRUE),
      rmse          = sqrt(mean(errors^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
))

gamma_performance <- gamma_performance[
  order(gamma_performance$aux_quality,
        gamma_performance$method,
        gamma_performance$gamma), ]


# -----------------------------
# 8. Save outputs
# -----------------------------

write.csv(
  gamma_raw,
  file = "outputs/tables/gamma_raw_results.csv",
  row.names = FALSE
)

write.csv(
  gamma_performance,
  file = "outputs/tables/gamma_performance.csv",
  row.names = FALSE
)

cat("\nDone. Bias by method across the gamma gradient:\n")
print(
  gamma_performance[, c("aux_quality", "method", "gamma", "mean_bias", "rmse")],
  row.names = FALSE
)
