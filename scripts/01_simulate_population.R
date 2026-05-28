##########################
# Purpose:
# Create a simulated finite population and draw probability and non-probability samples
# under a 2*2 design:
#
# 1. Auxiliary variable quality:
#   - strong: observed X variables strongly predict Y
#   - weak: observed X variables weakly predict Y
#
# 2. Selection mechanism:
#   - ignorable: non-probability sample selection depends only on X
#   - nonignorable: selection depends on X and Y
#
# Output:
#   - data/processed/simulated_samples.rds
#   - outputs/tables/scenario_summary.csv
###############################

# -----------------------------
# 0. Helper functions
# -----------------------------

standardize <- function(x){as.numeric(scale(x))}

calibrate_intercept <- function(linear_predictor, target_rate) {
  # Find an intercept so that mean(plogis(intercept + linear_predictor))
  # is approximately equal to target_rate.
  f <- function(intercept) {
    mean(plogis(intercept + linear_predictor)) - target_rate
  }
  
  uniroot(f, interval = c(-50, 50))$root
}

draw_bernoulli_sample <- function(probabilities) {
  rbinom(length(probabilities), size = 1, prob = probabilities)
}


# -----------------------------
# 1. Simulate finite population
# -----------------------------
set.seed(2026)

simulate_population <- function(n_population = 100000,
                                aux_quality = c("strong", "weak")){
  
  aux_quality <- match.arg(aux_quality)
  
 # Observed auxiliary variables X
  age <- round(runif(n_population, min = 18, max = 80))
  female <- rbinom(n_population, size = 1, prob = 0.52)
  higher_education <- rbinom(n_population, size = 1, prob = 0.35)
  urban <- rbinom(n_population, size = 1, prob = 0.65)
  
  # Unobserved variable U
  u_unobserved <- rnorm(n_population)
  
  # Standardized versions for outcome generation
  z_age <- standardize(age)
  z_female <- standardize(female)
  z_edu <- standardize(higher_education)
  z_urban <- standardize(urban)
  z_u <- standardize(u_unobserved)
  
 # Outcome Y
 # In the strong auxiliary scenario, observed X explains Y well.
 # In the weak auxiliary scenario, Y is driven more by unobserved U.
  
  if(aux_quality == "strong"){
    y <- 50 + 6.0 * z_age + 4.5 * z_edu + 2.5 * z_urban + 1.5 * z_female + 1.0 * z_u + 
      rnorm(n_population, mean = 0, sd = 6)
  }
  
  if(aux_quality == "weak") {
    y <- 50 + 1.0 * z_age + 0.8 * z_edu + 0.5 * z_urban + 0.3 * z_female + 6.0 * z_u + 
      rnorm(n_population, mean = 0, sd = 6)
  }
  
  population <- data.frame(
    id = seq_len(n_population),
    aux_quality = aux_quality,
    age = age, 
    female = female,
    higher_education = higher_education,
    urban = urban,
    u_unobserved = u_unobserved,
    y = y
  )
  
  return(population)
}

# -----------------------------
# 2. Add sampling indicators
# -----------------------------

add_sampling_indicators <- function(
    population,
    selection_mechanism = c("ignorable", "nonignorable"),
    n_probability_sample = 2000,
    expected_n_nonprob_sample = 3000){
  
  selection_mechanism <- match.arg(selection_mechanism)
  n_population <- nrow(population)
  
  # Probability sample:
  # simple random sample from the finite population.
  probability_sample <- rep(0, n_population)
  probability_sample[sample(seq_len(n_population), size = n_probability_sample)] <- 1
  
  # Non-probability sample:
  # selection probability is biased.
  z_age <- standardize(population$age)
  z_edu <- standardize(population$higher_education)
  z_urban <- standardize(population$urban)
  z_female <- standardize(population$female)
  z_y <- standardize(population$y)
  
  # Base selection tendency depending on observed X.
  # This creates overrepresentation of some groups.
  linear_predictor <- 0.50 * z_edu + 0.35 * z_urban - 0.25 * z_age + 0.15 * z_female
  
  # Non-ignorable selection:
  # selection also depends directly on the outcome Y.
  if(selection_mechanism == "nonignorable"){
    linear_predictor <- linear_predictor + 0.60 * z_y
  }
  
  target_rate <- expected_n_nonprob_sample / n_population
  intercept <- calibrate_intercept(linear_predictor, target_rate)
  
  nonprob_selection_probability <- plogis(intercept + linear_predictor)
  nonprob_sample <- draw_bernoulli_sample(nonprob_selection_probability)
  
  population$selection_mechanism <- selection_mechanism
  population$probability_sample <- probability_sample
  population$nonprob_sample <- nonprob_sample
  population$nonprob_selection_probability <- nonprob_selection_probability
  
  return(population)
  
}


# -----------------------------
# 3. Generate 2 x 2 scenarios
# -----------------------------

auxiliary_scenarios <- c("strong", "weak")
selection_scenarios <- c("ignorable", "nonignorable")

simulated_scenarios <- list()

for (aux in auxiliary_scenarios){
  base_population <- simulate_population(aux_quality = aux)
  
  for (selection in selection_scenarios){
    scenario_name <- paste(aux, selection, sep = "_")
    
    scenario_data <- add_sampling_indicators(
      population = base_population,
      selection_mechanism = selection
    )
    
    simulated_scenarios[[scenario_name]] <- scenario_data
  }
}

# -----------------------------
# 4. Create scenario summary
# -----------------------------

scenario_summary <- do.call(
  rbind,
  lapply(names(simulated_scenarios), function(name) {
    dat <- simulated_scenarios[[name]]
    
    data.frame(
      scenario = name,
      aux_quality = unique(dat$aux_quality),
      selection_mechanism = unique(dat$selection_mechanism),
      population_n = nrow(dat),
      probability_sample_n = sum(dat$probability_sample),
      nonprob_sample_n = sum(dat$nonprob_sample),
      population_mean_y = mean(dat$y),
      probability_sample_mean_y = mean(dat$y[dat$probability_sample == 1]),
      nonprob_naive_mean_y = mean(dat$y[dat$nonprob_sample == 1])
    )
  })
)

scenario_summary$probability_sample_bias <-
  scenario_summary$probability_sample_mean_y -
  scenario_summary$population_mean_y

scenario_summary$nonprob_naive_bias <-
  scenario_summary$nonprob_naive_mean_y -
  scenario_summary$population_mean_y


# -----------------------------
# 5. Save outputs
# -----------------------------

saveRDS(
  simulated_scenarios,
  file = "data/processed/simulated_samples.rds"
)

write.csv(
  scenario_summary,
  file = "outputs/tables/scenario_summary.csv",
  row.names = FALSE
)







