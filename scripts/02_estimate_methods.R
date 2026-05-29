############################################################
# 02_estimate_methods.R
#
# Purpose:
#   Estimate the population mean of Y using several methods
#   for each simulated scenario created by:
#
#   scripts/01_simulate_population.R
#
# Methods implemented:
#   1. Probability sample mean
#   2. Naive non-probability sample mean
#   3. Calibration / balancing weights
#   4. Estimated inverse probability weighting (IPW)
#   5. Oracle IPW using true simulated selection probabilities
#   6. Mass imputation / outcome prediction
#   7. Doubly robust estimator
#
# Output:
#   - outputs/tables/estimator_comparison.csv
############################################################

# -----------------------------
# 1. Load simulated scenarios
# -----------------------------

simulated_scenarios <- readRDS("data/processed/simulated_samples.rds")

# -----------------------------
# 2. Helper functions
# -----------------------------

clip_probabilities <- function(p, lower = 0.001, upper = 0.999){
  pmin(pmax(p, lower), upper)
}

# Calibration estimator
#
# This function directly adjusts non-probability sample weights so that
# weighted auxiliary-variable totals match the finite population totals.
#
# This is a simple linear calibration implementation.
# It is useful for simulation and learning purposes.
calibration_mean <- function(population, sample_indicator) {
  sample_data <- population[sample_indicator == 1, ]
 
  # Calibration variables
  # Include intercept automatically
 x_population <- model.matrix(
  ~ age + female + higher_education + urban, 
  data = population
 )
 
 x_sample <- model.matrix(
   ~ age + female + higher_education + urban,
   data = sample_data
 )
 
 # Population totals of auxiliary variables
 population_totals <- colSums(x_population)
 
 # Base weights for non-probability sample
 base_weights <- rep(1, nrow(sample_data))
 
 # Sample totals under base weights
 sample_totals <- colSums(x_sample * base_weights)
 
 # Calibration adjustment
 #
 # We solve for lambda in:
 #  sum_s w_i x_i = X_N
 # with
 #  w_i = d_i * (1 + x_i' lambda)
 #
 # where d_i are base weights
 xtx <- t(x_sample) %*% x_sample
 difference <- population_totals - sample_totals
 
 lambda <- tryCatch(
   solve(xtx, difference),
   error = function(e) { return(rep(NA_real_, ncol(x_sample)))}
 )
 if (any(is.na(lambda))) {
   return(
     list(
       estimate = NA_real_,
       min_weight = NA_real_,
       max_weight = NA_real_
     )
   )
 }
 
 calibration_factors <- as.vector(1 + x_sample %*% lambda)
 calibrated_weights <- base_weights * calibration_factors
 
 estimate <- sum(calibrated_weights * sample_data$y) / nrow(population)
 
 return(
   list(
     estimate = estimate,
     min_weight = min(calibrated_weights),
     max_weight = max(calibrated_weights)
   )
 )
 
 }

# Estimated IPW estimator
#
# This estimates P(S = 1 | X) using a logistic regression model
# In this simulation, we use the full finite population X information
# to estimate selection probabilities. This reflects a setting where
# the population frame or rich benchmark data are available
# 
# The model uses only observed X variables, not Y
# Therefore, it is expected to be limited under non-ignorable selection


estimated_ipw_mean <- function(population, sample_indicator){
  propensity_model <- glm(
    sample_indicator ~ age + female + higher_education + urban,
    data = population,
    family = binomial()
  )
  
  p_hat <- predict(propensity_model, type = "response")
  p_hat <- clip_probabilities(p_hat)
  
  sample_data <- population[population$sample_indicator == 1, ]
  sample_p_hat <- p_hat[sample_indicator == 1]
  
  estimate <- sum(sample_data$y / sample_p_hat) / nrow(population)
  
  return(estimate)
}


# Oracle IPW estimator
# 
# This uses the true simulated non-probability selection probabilities
# It is not available in real applications
# It is included only as a simulation benchmark
oracle_ipw_mean <- function(population, sample_indicator){
  sample_data <- population[sample_indicator == 1,]
  sample_true_p <- population$nonprob_selection_probability[sample_indicator == 1]
  sample_true_p <- clip_probabilities(sample_true_p)
  
  estimate <- sum(sample_data$y / sample_true_p) / nrow(population)
  
  return(estimate)
}























