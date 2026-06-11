############################################################
# 08_plot_gamma_results.R
#
# Purpose:
#   Visualize results from 07_nonignorability_gradient.R:
#   mean bias of each integration estimator as a function of the
#   non-ignorability parameter gamma, faceted by aux_quality.
#
# Input:
#   - outputs/tables/gamma_performance.csv
#
# Output:
#   - outputs/figures/gamma_bias_plot.png
############################################################

library(ggplot2)

gamma_performance <- read.csv("outputs/tables/gamma_performance.csv")

# Methods to plot. mass_imputation_integration is omitted because it
# is numerically identical to calibrated_integration under the linear
# model used here (see notes/round1-2).
methods_to_plot <- c(
  "np_naive_mean",
  "calibrated_integration",
  "membership_ipw",
  "doubly_robust_integration",
  "prob_sample_mean"
)

method_labels <- c(
  np_naive_mean              = "Naive (unadjusted)",
  calibrated_integration     = "Calibration",
  membership_ipw             = "IPW",
  doubly_robust_integration  = "Doubly robust",
  prob_sample_mean           = "Probability sample (reference)"
)

plot_data <- subset(gamma_performance, method %in% methods_to_plot)
plot_data$method <- factor(plot_data$method,
                            levels = methods_to_plot,
                            labels = method_labels[methods_to_plot])

plot_data$aux_quality <- factor(
  plot_data$aux_quality,
  levels = c("strong", "weak"),
  labels = c("Strong auxiliary information", "Weak auxiliary information")
)

p <- ggplot(plot_data, aes(x = gamma, y = mean_bias,
                            color = method, linetype = method)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  facet_wrap(~ aux_quality) +
  scale_linetype_manual(values = c(
    "Naive (unadjusted)" = "dashed",
    "Calibration" = "solid",
    "IPW" = "solid",
    "Doubly robust" = "solid",
    "Probability sample (reference)" = "dotted"
  )) +
  labs(
    x = expression("Non-ignorability parameter " * gamma),
    y = "Mean bias",
    color = NULL,
    linetype = NULL,
    title = "Bias of integration estimators as selection becomes increasingly outcome-dependent",
    subtitle = "Mass imputation omitted (overlaps calibration exactly)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  filename = "outputs/figures/gamma_bias_plot.png",
  plot = p,
  width = 10, height = 4.5, dpi = 150
)

print(p)
