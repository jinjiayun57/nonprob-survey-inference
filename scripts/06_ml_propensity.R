############################################################
# 06_ml_propensity.R
#
# Purpose:
#   Compare logistic, random-forest (ranger), and gradient-
#   boosted (xgboost) propensity models for reference-sample
#   IPW and doubly robust (DR) estimation across the 2x2
#   simulation scenarios.
#
#   Population is fixed (generated once per scenario).
#   Each replicate draws new probability and non-probability
#   samples. The outcome model (linear) is held fixed across
#   all methods to isolate the effect of propensity model
#   choice.
#
#   Extends 04_monte_carlo_simulations.R.
#   Same population seed (2026), replicate seed (42), and
#   sample sizes as 04.
#
# Requirements:
#   install.packages(c("ranger", "xgboost"))
#
# Output:
#   - outputs/tables/ml_propensity_raw.csv
#   - outputs/tables/ml_propensity_performance.csv
############################################################

library(ranger)
library(xgboost)

# -------------------------------------------------------
# 0. Reuse helpers from 01 and 03
# -------------------------------------------------------

source("scripts/01_simulate_population.R")
rm(simulated_scenarios, scenario_summary)

source("scripts/03_reference_sample_integration.R")
rm(integration_results)

# Functions now in scope:
#   simulate_population(), add_sampling_indicators()
#   clip_probabilities(), weighted_mean()
#   mass_imputation_integration_mean()

# -------------------------------------------------------
# 1. Settings
# -------------------------------------------------------

N_SIM     <- 500
N_POP     <- 100000
N_PROB    <- 2000
N_NONPROB <- 3000

X_VARS <- c("age", "female", "higher_education", "urban")

RF_NUM_TREES   <- 500
RF_MIN_NODE    <- 10

XGB_NROUNDS <- 100

SCENARIOS <- list(
  strong_ignorable = list(aux = "strong", sel = "ignorable"),
  strong_nonignorable = list(aux = "strong", sel = "nonignorable"),
  weak_ignorable  = list(aux = "weak",   sel = "ignorable"),
  weak_nonignorable = list(aux = "weak",   sel = "nonignorable")
)

# -------------------------------------------------------
# 2. Fixed populations
# -------------------------------------------------------

set.seed(2026)

fixed_population <- list(
  strong = simulate_population(n_population = N_POP, aux_quality = "strong"),
  weak   = simulate_population(n_population = N_POP, aux_quality = "weak")
)

# -------------------------------------------------------
# 3. Shared weight utilities
# -------------------------------------------------------

# Hájek inverse-odds weights from P(np | X)

ipw_weights <- function(pi_hat_np){
  pi_hat_np <- clip_probabilities(pi_hat_np)
  (1 - pi_hat_np) / pi_hat_np
}

weight_diagnostics <- function(w){
  list(
    min_weight = min(w),
    max_weight = max(w),
    ess = sum(w)^2 / sum(w^2)
  )
}

# -------------------------------------------------------
# 4. Propensity model functions
#
# Each takes p_sample and np_sample; returns a list with:
#   weights    : Hájek inverse-odds weights for np units
#   min_weight, max_weight, ess : diagnostics
# -------------------------------------------------------

propensity_logistic <- function(p_sample, np_sample){
  pooled <- rbind(
    data.frame(p_sample[, X_VARS], R = 0),
    data.frame(np_sample[, X_VARS], R = 1)
  )
  fit <- glm(R ~ age + female + higher_education + urban,
             data = pooled, family = binomial())
  pi_hat_all <- predict(fit, type = "response")
  pi_hat_np <- pi_hat_all[(nrow(p_sample) + 1) : nrow(pooled)]
  w <- ipw_weights(pi_hat_np)
  c(list(weights = w), weight_diagnostics(w))
}

propensity_rf <- function(p_sample, np_sample){
  pooled <- rbind(
    data.frame(p_sample[, X_VARS], R = as.factor(0)),
    data.frame(np_sample[, X_VARS], R = as.factor(1))
  )
  fit <- ranger(
    R ~ age + female + higher_education + urban,
    data = pooled,
    probability = TRUE,
    num.trees = RF_NUM_TREES,
    min.node.size = RF_MIN_NODE
  )
  
  # predict() returns a matrix; column "1" = P(np | X)
  pi_hat_np <- predict(fit, data = np_sample[, X_VARS])$predictions[, "1"]
  w <- ipw_weights(pi_hat_np)
  c(list(weights = w), weight_diagnostics(w))
}

propensity_xgb <- function(p_sample, np_sample) {
  X_pooled <- rbind(
    p_sample[, X_VARS],
    np_sample[, X_VARS]
  )
  
  y_pooled <- factor(c(rep(0, nrow(p_sample)), rep(1, nrow(np_sample))),
                     levels = c(0, 1),
                     labels = c("p", "np"))
  
  fit <- xgboost(
    x = X_pooled,
    y = y_pooled,
    objective = "binary:logistic",
    nrounds = XGB_NROUNDS,
    max_depth = 4,
    learning_rate = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    nthreads = 1,
    verbosity = 0
  )
  
  pi_hat_np <- predict(fit, 
                       newdata = np_sample[, X_VARS],
                       type = "response")
  
  w <- ipw_weights(pi_hat_np)
  c(list(weights = w), weight_diagnostics(w))
}

PROPENSITY_MODELS <- list(
  logistic = propensity_logistic,
  rf       = propensity_rf,
  xgb      = propensity_xgb
)

# -------------------------------------------------------
# 5. IPW and DR estimators
#
# Both use Hájek-normalized weights.
# DR outcome model: linear, fitted in np-sample.
# The outcome model is held fixed across all propensity
# model choices to isolate the propensity model effect.
# -------------------------------------------------------

ipw_estimate <- function(y_np, w){
  weighted_mean(y_np, w)
}

dr_estimate <- function(p_sample, np_sample, w){
  outcome_model <- lm(
    y ~ age + female + higher_education + urban,
    data = np_sample
  )
  m_hat_p <- predict(outcome_model, newdata = p_sample)
  m_hat_np <- predict(outcome_model, newdata = np_sample)
  
  mean(m_hat_p) + weighted_mean(np_sample$y - m_hat_np, w)
}

# -------------------------------------------------------
# 6. Single-replicate function
# -------------------------------------------------------

run_one_replicate_ml <- function(rep_id, base_population, 
                                 selection_mechanism, scenario_name){
  population <- add_sampling_indicators(
    population = base_population,
    selection_mechanism = selection_mechanism,
    n_probability_sample = N_PROB,
    expected_n_nonprob_sample = N_NONPROB
  )

true_mean <- mean(population$y)
aux_quality <- unique(population$aux_quality)
p_sample <- population[population$probability_sample == 1, ]
np_sample <- population[population$nonprob_sample == 1, ]

# --- Benchmarks (no propensity model needed) ---
benchmarks <- data.frame(
  rep_id = rep_id,
  scenario = scenario_name,
  aux_quality = aux_quality,
  selection_mechanism = selection_mechanism,
  method = c("prob_sample_mean", "np_naive_mean"),
  estimate = c(mean(p_sample$y), mean(np_sample$y)),
  true_mean = true_mean,
  min_weight          = NA_real_,
  max_weight          = NA_real_,
  ess                 = NA_real_,
  stringsAsFactors    = FALSE
)

# --- IPW and DR for each propensity model ---

ml_rows <- do.call(rbind, lapply(names(PROPENSITY_MODELS), function(model_name){
  
  prop <- tryCatch(
    PROPENSITY_MODELS[[model_name]](p_sample, np_sample),
    error = function(e) {
      warning("Propensity model '", model_name, "' failed in rep ", rep_id,
              ": ", conditionMessage(e))
      NULL
    }
  )
methods_out <- paste0(model_name, c("_ipw", "_dr"))

if (is.null(prop)) {
  return(data.frame(
    rep_id = rep_id, scenario = scenario_name,
    aux_quality = aux_quality, selection_mechanism = selection_mechanism,
    method = methods_out, estimate = NA_real_, true_mean = true_mean,
    min_weight = NA_real_, max_weight = NA_real_, ess = NA_real_,
    stringsAsFactors = FALSE
  ))
}

estimates <- c(
  ipw_estimate(np_sample$y, prop$weights),
  dr_estimate(p_sample, np_sample, prop$weights)
)

data.frame(
  rep_id              = rep_id,
  scenario            = scenario_name,
  aux_quality         = aux_quality,
  selection_mechanism = selection_mechanism,
  method              = methods_out,
  estimate            = estimates,
  true_mean           = true_mean,
  min_weight          = prop$min_weight,
  max_weight          = prop$max_weight,
  ess                 = prop$ess,     # same weights for both rows
  stringsAsFactors    = FALSE
)
}))

rbind(benchmarks, ml_rows)
}

# -------------------------------------------------------
# 7. Monte Carlo loop
# -------------------------------------------------------

set.seed(42)

ml_raw <- do.call(rbind, lapply(names(SCENARIOS), function(scenario_name){
  cfg <- SCENARIOS[[scenario_name]]
  base_pop <- fixed_population[[cfg$aux]]
  
  cat("Running:", scenario_name, "\n")
  
  do.call(rbind, lapply(seq_len(N_SIM), function(rep_id){
    run_one_replicate_ml(
      rep_id = rep_id,
      base_population = base_pop,
      selection_mechanism = cfg$sel,
      scenario_name = scenario_name
    )
  }))
}))

# -------------------------------------------------------
# 8. Performance metrics
# -------------------------------------------------------

ml_performance <- do.call(rbind, lapply(
  split(ml_raw, list(ml_raw$scenario, ml_raw$method)),
  function(df){
    errors <- df$estimate - df$true_mean
    data.frame(
      scenario            = df$scenario[1],
      aux_quality         = df$aux_quality[1],
      selection_mechanism = df$selection_mechanism[1],
      method              = df$method[1],
      n_replicates        = sum(!is.na(df$estimate)),
      true_mean           = df$true_mean[1],
      mean_estimate       = mean(df$estimate,  na.rm = TRUE),
      mean_bias           = mean(errors,        na.rm = TRUE),
      mae                 = mean(abs(errors),   na.rm = TRUE),
      variance            = var(df$estimate,    na.rm = TRUE),
      rmse                = sqrt(mean(errors^2, na.rm = TRUE)),
      mean_ess            = mean(df$ess,        na.rm = TRUE),
      mean_max_weight     = mean(df$max_weight, na.rm = TRUE),
      stringsAsFactors    = FALSE
    )
  }
))

ml_performance <- ml_performance[
  order(ml_performance$scenario, ml_performance$method), ]

# -------------------------------------------------------
# 9. Save outputs
# -------------------------------------------------------

write.csv(ml_raw,
          "outputs/tables/ml_propensity_raw.csv",
          row.names = FALSE)

write.csv(ml_performance,
          "outputs/tables/ml_propensity_performance.csv",
          row.names = FALSE)

cat("\nDone. Performance summary:\n")
print(
  ml_performance[, c("scenario", "method", "mean_bias",
                     "mae", "rmse", "mean_ess")],
  row.names = FALSE
)
