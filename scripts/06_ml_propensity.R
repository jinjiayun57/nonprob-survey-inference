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

                                                                                                                                                                                                                                        DR 估计量的逻辑是两项之和：
                                                                                                                                                                                                                                             μ^DR=1np∑i∈spm^(Xi)⏟prediction component+∑i∈snpwi (Yi−m^(Xi))∑i∈snpwi⏟residual correction\hat{\mu}_{DR} = \underbrace{\frac{1}{n_p}\sum_{i \in s_p} \hat{m}(X_i)}_{\text{prediction component}} + \underbrace{\frac{\sum_{i \in s_{np}} w_i \,(Y_i - \hat{m}(X_i))}{\sum_{i \in s_{np}} w_i}}_{\text{residual correction}}μ^​DR​=prediction componentnp​1​i∈sp​∑​m^(Xi​)​​+residual correction∑i∈snp​​wi​∑i∈snp​​wi​(Yi​−m^(Xi​))​​​
                                                                                                                                                                                                                                             第一项：在 p-sample 上对结果模型 m^(X)\hat{m}(X)
                                                                                                                                                                                                                                             m^(X) 取均值——这是 mass imputation 部分。
                                                                                                                                                                                                                                             第二项：用 同一套 IPW 权重 对 np-sample 的残差做 Hájek 加权均值——这是 propensity 部分。
                                                                                                                                                                                                                                             所以 propensity_rf / propensity_xgb 返回的 prop$weights 同时传给了 ipw_estimate（只用第二项）和 dr_estimate（用第二项做残差修正）。两个估计量共用同一个 w，这正是脚本里 ess 那一列对 _ipw 和 _dr 行数值相同的原因。
                                                                                                                                                                                                                                             DR 的"双重稳健"性质来自：只要 m^\hat{m}
                                                                                                                                                                                                                                             m^ 和 ww
                                                                                                                                                                                                                                             w 中有一个接近正确，估计量就趋于一致——outcome model 误设时 propensity 兜底，propensity 误设时 outcome model 兜底。You said: 所以这个ipw和dr都是基于上面machine learning的方法算出来的propensity求的, 不是logistic model求出来的propensity所以这个ipw和dr都是基于上面machine learning的方法算出来的propensity求的, 不是logistic model求出来的propensityClaude responded: 对。对。
                                                                                                                                                                                                                                             三个 propensity_* 函数各自独立估计 π^i\hat{\pi}_i
                                                                                                                                                                                                                                             π^i​：
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             propensity_logistic → logistic 回归的 π^\hat{\pi}
                                                                                                                                                                                                                                             π^
                                                                                                                                                                                                                                               propensity_rf → ranger 概率森林的 π^\hat{\pi}
                                                                                                                                                                                                                                             π^
                                                                                                                                                                                                                                               propensity_xgb → xgboost 的 π^\hat{\pi}
                                                                                                                                                                                                                                             π^
                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                               然后 ipw_estimate 和 dr_estimate 只是接收 w（已经算好的权重），不知道也不关心这个 w 是哪个模型算出来的。
                                                                                                                                                                                                                                             所以在 run_one_replicate_ml 的循环里：
                                                                                                                                                                                                                                             rlapply(names(PROPENSITY_MODELS), function(model_name) {
                                                                                                                                                                                                                                               prop <- PROPENSITY_MODELS[[model_name]](p_sample, np_sample)
                                                                                                                                                                                                                                               ipw_estimate(np_sample$y, prop$weights)   # 用这个模型的 w
                                                                                                                                                                                                                                               dr_estimate(p_sample, np_sample, prop$weights)  # 用同一个 w
                                                                                                                                                                                                                                             })
                                                                                                                                                                                                                                             每次迭代用的是当前那个模型的 prop$weights，互不干扰。最终输出的 rf_ipw / rf_dr 全程只用了 RF 的 propensity，xgb_ipw / xgb_dr 全程只用了 XGB 的 propensity。You said: 啊好的 我还以为这里的ipw和dr是我们以前传统方法里那种ipw和dr就是logistic的, 看来全部更新成新方法了啊好的 我还以为这里的ipw和dr是我们以前传统方法里那种ipw和dr就是logistic的, 看来全部更新成新方法了Claude responded: 对，可以这样理解脚本的结构：对，可以这样理解脚本的结构：
                                                                                                                                                                                                                                             04 里（传统方法）：
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             membership_ipw = logistic propensity + Hájek weights
                                                                                                                                                                                                                                             doubly_robust_integration = logistic propensity + linear outcome
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             06 里：
                                                                                                                                                                                                                                             把 propensity 估计这一步抽象成了可替换的模块，所以：
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             logistic_ipw / logistic_dr = 和 04 完全一样，作为 baseline
                                                                                                                                                                                                                                             rf_ipw / rf_dr = RF propensity + 同样的 Hájek / DR 框架
                                                                                                                                                                                                                                             xgb_ipw / xgb_dr = XGB propensity + 同样的 Hájek / DR 框架
                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                             outcome model（linear lm）在所有 DR 变体里保持不变，这样比较结果的差异就能归因到 propensity 模型本身，而不是两个东西同时变。
                                                                                                                                                                                                                                             所以跑完之后，如果 rf_dr 比 logistic_dr 在某个场景下 bias 更小，就说明 RF 在那个场景里对 propensity 的刻画更准确。You said: runonereplicateml(1, basepopulation, "ignorable", "weak_ignorable")run_one_replicate_ml(1, base_population, "ignorable", "weak_ignorable")
                                                                                                                                                                                                                                             Warning messages:
                                                                                                                                                                                                                                               1: In throw_err_or_depr_msg("Parameter(s) have been removed from this function: ",  :
                                                                                                                                                                                                                                                                             Parameter(s) have been removed from this function: params. This warning will become an error in a future version.
                                                                                                                                                                                                                                                                           2: In throw_err_or_depr_msg("Passed unrecognized parameters: ", paste(head(names_unrecognized),  :
                                                                                                                                                                                                                                                                                                                                                   Passed unrecognized parameters: verbose. This warning will become an error in a future version.
                                                                                                                                                                                                                                                                                                                                                 3: In throw_err_or_depr_msg("Parameter '", match_old, "' has been renamed to '",  :
                                                                                                                                                                                                                                                                                                                                                                               Parameter 'data' has been renamed to 'x'. This warning will become an error in a future version.
                                                                                                                                                                                                                                                                                                                                                                             4: In value[[3L]](cond) :
                                                                                                                                                                                                                                                                                                                                                                               Propensity model 'xgb' failed in rep 1: argument "y" is missing, with no defaultShow moreClaude responded: xgboost 2.Thought processThought processxgboost 2.0+ 改了 API：xgboost() 高层函数重写了，params/verbose/data 这些参数都变了。修复方法是换用底层的 xgb.train()，它的接口更稳定。
                                                                                                                                                                                                                                                                                                                                                                             把 XGB_PARAMS 里的 nthread 保留，去掉任何 verbose，然后把 propensity_xgb 改成这样：
                                                                                                                                                                                                                                                                                                                                                                             rXGB_PARAMS <- list(
                                                                                                                                                                                                                                                                                                                                                                               objective        = "binary:logistic",
                                                                                                                                                                                                                                                                                                                                                                               eta              = 0.1,
                                                                                                                                                                                                                                                                                                                                                                               max_depth        = 4,
                                                                                                                                                                                                                                                                                                                                                                               subsample        = 0.8,
                                                                                                                                                                                                                                                                                                                                                                               colsample_bytree = 0.8,
                                                                                                                                                                                                                                                                                                                                                                               nthread          = 1
                                                                                                                                                                                                                                                                                                                                                                             )
                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                             propensity_xgb <- function(p_sample, np_sample) {
                                                                                                                                                                                                                                                                                                                                                                               X_pooled <- as.matrix(rbind(
                                                                                                                                                                                                                                                                                                                                                                                 p_sample[, X_VARS],
                                                                                                                                                                                                                                                                                                                                                                                 np_sample[, X_VARS]
                                                                                                                                                                                                                                                                                                                                                                               ))
                                                                                                                                                                                                                                                                                                                                                                               y_pooled <- c(rep(0L, nrow(p_sample)), rep(1L, nrow(np_sample)))
                                                                                                                                                                                                                                                                                                                                                                               dtrain   <- xgboost::xgb.DMatrix(data = X_pooled, label = y_pooled)
                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                               fit <- xgboost::xgb.train(
                                                                                                                                                                                                                                                                                                                                                                                 params  = XGB_PARAMS,
                                                                                                                                                                                                                                                                                                                                                                                 data    = dtrain,
                                                                                                                                                                                                                                                                                                                                                                                 nrounds = XGB_NROUNDS,
                                                                                                                                                                                                                                                                                                                                                                                 verbose = 0          # verbose 是 xgb.train() 自己的参数，不放进 params
                                                                                                                                                                                                                                                                                                                                                                               )
                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                               X_np  <- as.matrix(np_sample[, X_VARS])
                                                                                                                                                                                                                                                                                                                                                                               dtest <- xgboost::xgb.DMatrix(data = X_np)
                                                                                                                                                                                                                                                                                                                                                                               pi_hat_np <- predict(fit, newdata = dtest)
                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                               w <- ipw_weights(pi_hat_np)
                                                                                                                                                                                                                                                                                                                                                                               c(list(weights = w), weight_diagnostics(w))
                                                                                                                                                                                                                                                                                                                                                                             }
                                                                                                                                                                                                                                                                                                                                                                             核心变化两处：xgboost() → xgb.train()，prediction 的输入也包成 xgb.DMatrix。xgb.train 的 params + data(DMatrix) + nrounds 这套接口在各版本里一直稳定。You said: Warning messages:Warning messages:
                                                                                                                                                                                                                                                                                                                                                                               1: In throw_err_or_depr_msg("Parameter(s) have been removed from this function: ",  :
                                                                                                                                                                                                                                                                                                                                                                                                             Parameter(s) have been removed from this function: params. This warning will become an error in a future version.
                                                                                                                                                                                                                                                                                                                                                                                                           2: In throw_err_or_depr_msg("Passed unrecognized parameters: ", paste(head(names_unrecognized),  :
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   Passed unrecognized parameters: verbose. This warning will become an error in a future version.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 3: In throw_err_or_depr_msg("Parameter '", match_old, "' has been renamed to '",  :
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Parameter 'data' has been renamed to 'x'. This warning will become an error in a future version.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             4: In value[[3L]](cond) :
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Propensity model 'xgb' failed in rep 1: argument "y" is missing, with no default
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             >还是这样啊. 你直接来我项目里debugShow moreClaude responded: API 全清楚了。问题有两个来源：Searched the web, ran 2 commands, loaded toolsSearched the web, ran 2 commands, loaded toolsAPI 全清楚了。问题有两个来源：
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             问题1：xgboost() 高层函数在新版重写了，params/data/verbose 全部变了。但 xgb.train() 本身没变，只是 params 现在要用 xgb.params() 创建，不能是普通 list()。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             问题2：xgb.DMatrix 的 data 参数改名为 x，作为命名参数传会报警。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             把 XGB_PARAMS 和 propensity_xgb 改成这样：
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             r# 删掉原来的 XGB_PARAMS <- list(...)，改成下面这个
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             XGB_PARAMS <- xgboost::xgb.params(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               objective        = "binary:logistic",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               eta              = 0.1,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               max_depth        = 4L,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               subsample        = 0.8,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               colsample_bytree = 0.8,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               nthread          = 1L
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             )
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             propensity_xgb <- function(p_sample, np_sample) {
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               X_pooled <- as.matrix(rbind(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 p_sample[, X_VARS],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 np_sample[, X_VARS]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               y_pooled <- c(rep(0, nrow(p_sample)), rep(1, nrow(np_sample)))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               # 第一个参数用位置传入，避免 data= 命名触发 rename 警告
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               dtrain <- xgboost::xgb.DMatrix(X_pooled, label = y_pooled)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               fit <- xgboost::xgb.train(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 params  = XGB_PARAMS,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 data    = dtrain,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 nrounds = XGB_NROUNDS,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 verbose = 0          # verbose 是 xgb.train 自己的参数，不进 params
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               )
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               X_np  <- as.matrix(np_sample[, X_VARS])
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               dtest <- xgboost::xgb.DMatrix(X_np)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               pi_hat_np <- predict(fit, dtest)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               w <- ipw_weights(pi_hat_np)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               c(list(weights = w), weight_diagnostics(w))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             核心三处：list() → xgb.params()；xgb.DMatrix 的 data 改为位置传入；verbose 留在 xgb.train() 外层而不是 params 里。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Sources:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Fit XGBoost Model — xgb.train
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             You said: propensityxgb <- function(psample, np_sample) {propensity_xgb <- function(p_sample, np_sample) {
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               X_pooled <- rbind(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 p_sample[, X_VARS],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 np_sample[, X_VARS]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               )  y_pooled <- factor(c(rep(0, nrow(p_sample)), rep(1, nrow(np_sample))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     levels = c(0, 1),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     labels = c("p", "np"))  fit <- xgboost(
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     )  pi_hat_np <- predict(fit, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             newdata = np_sample[, X_VARS],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             type = "response")  w <- ipw_weights(pi_hat_np)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             c(list(weights = w), weight_diagnostics(w))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             我改成了这样 可以了Show lessClaude responded: 新的高层接口用起来确实更干净。几点值得记一下：Thought processThought process新的高层接口用起来确实更干净。几点值得记一下：
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             eta 改名成了 learning_rate，nthread 变成了 nthreads——这是新版向 Python 接口对齐的结果。y 用 factor 传入，predict(..., type = "response") 返回的就是 P(class = "np")，也就是你要的 π^i\hat{\pi}_i
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             π^i​，逻辑上完全正确。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             你也不需要 as.matrix() 了，新接口直接接受 data frame，省了一步。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             可以把 XGB_PARAMS 这个全局变量也删掉，因为现在参数直接写在函数里了，保持一致就好。You said: mc 循环预计多久可以跑好, 一个scenario 都跑了两三分钟了还没好mc 循环预计多久可以跑好, 一个scenario 都跑了两三分钟了还没好Claude responded: RF 是瓶颈。每次 replicate 要在 ~5000 行（2000 p-sample + 3000 np-sample）上跑 500 棵树，单核大概 0.Thought processThought processRF 是瓶颈。每次 replicate 要在 ~5000 行（2000 p-sample + 3000 np-sample）上跑 500 棵树，单核大概 0.3–0.8 秒。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             粗算：500 reps × ~0.5s × 4 scenarios ≈ 15–30 分钟。如果现在跑了 2–3 分钟还没完一个 scenario，说明每 rep 偏慢，总时长可能到 30 分钟以上。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             两个办法可以加速：
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             1. 减少 RF 树数（最直接）
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             rRF_NUM_TREES <- 200   # 从 500 降到 200，速度提升约 2.5x
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             4 个特征的数据集，200 棵树已经足够稳定。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2. 开多线程
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ranger 和 xgboost 都支持并行。把线程数从 1 改成你机器的核数：
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             rRF_NUM_THREADS <- parallel::detectCores() - 1  # 或直接写 4/6/8
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             # ranger 里加 num.threads = RF_NUM_THREADS
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             fit <- ranger::ranger(..., num.threads = RF_NUM_THREADS)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             # xgboost 里 nthreads 改成同样的值
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             fit <- xgboost(..., nthreads = RF_NUM_THREADS)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             注意 xgboost 之前设 nthreads = 1 是为了 reproducibility，改多线程后同一 seed 结果可能有微小浮点差异，但对 MC 均值影响可以忽略。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             现在正在跑的这次建议等它跑完看结果，之后调参再重跑。AskSonnet 4.6









































