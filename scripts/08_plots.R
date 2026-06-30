############################################################
# 08_plots.R
#
# Purpose:
#   All visualisations for the non-probability simulation project.
#   Merged from former scripts 08–10 and 12.
#
# Sections:
#   1. Gamma gradient – bias & RMSE line plots
#   2. Gamma gradient – distribution (violin + boxplot)
#   3. Monte Carlo scenarios – RMSE dot plot
#   4. Monte Carlo scenarios – distribution (violin + boxplot)
#
# Inputs:
#   outputs/tables/gamma_performance.csv
#   outputs/tables/gamma_raw_results.csv
#   outputs/tables/mc_performance.csv
#   outputs/tables/mc_raw_results.csv
#
# Outputs:
#   outputs/figures/gamma_bias_plot.png
#   outputs/figures/gamma_rmse_plot.png
#   outputs/figures/gamma_distribution_plot.png
#   outputs/figures/scenario_rmse_plot.png
#   outputs/figures/mc_violin_plot.png
############################################################

library(ggplot2)
library(patchwork)

# ── Shared palette ────────────────────────────────────────────────────────────
METHODS_CORE <- c(
  "np_naive_mean",
  "calibrated_integration",
  "membership_ipw",
  "doubly_robust_integration",
  "prob_sample_mean"
)

COLORS <- c(
  "Naive (unadjusted)"             = "grey50",
  "Naive"                          = "#999999",
  "Calibration / mass imputation"  = "#3aa655",
  "Calibration"                    = "#3aa655",
  "Calibration / MI"               = "#3aa655",
  "IPW"                            = "#e07b00",
  "Doubly robust"                  = "#6a5acd",
  "Probability sample (reference)" = "black",
  "Prob sample (ref)"              = "#222222"
)


# ══════════════════════════════════════════════════════════════════════════════
# 1. Gamma gradient – bias & RMSE line plots
# ══════════════════════════════════════════════════════════════════════════════

gamma_performance <- read.csv("outputs/tables/gamma_performance.csv")

method_labels_gamma <- c(
  np_naive_mean              = "Naive (unadjusted)",
  calibrated_integration     = "Calibration / mass imputation",
  membership_ipw             = "IPW",
  doubly_robust_integration  = "Doubly robust",
  prob_sample_mean           = "Probability sample (reference)"
)

method_linetypes <- c(
  "Naive (unadjusted)"             = "dashed",
  "Calibration / mass imputation"  = "solid",
  "IPW"                            = "solid",
  "Doubly robust"                  = "solid",
  "Probability sample (reference)" = "dotted"
)

gd <- subset(gamma_performance, method %in% METHODS_CORE)
gd$method <- factor(gd$method,
                    levels = METHODS_CORE,
                    labels = method_labels_gamma[METHODS_CORE])
gd$aux_quality <- factor(gd$aux_quality,
                         levels = c("strong", "weak"),
                         labels = c("Strong auxiliary information",
                                    "Weak auxiliary information"))

# Bias plot
p_gamma_bias <- ggplot(gd, aes(x = gamma, y = mean_bias,
                                color = method, linetype = method,
                                shape = method)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0.6, color = "grey70", linewidth = 0.4,
             linetype = "dotted") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ aux_quality) +
  scale_color_manual(values = COLORS) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = rep(16, 5),
                     breaks = levels(gd$method)) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_y_continuous(limits = c(NA, ceiling(max(gd$mean_bias, na.rm = TRUE)))) +
  labs(x = expression("Non-ignorability parameter " * gamma),
       y = "Mean bias of\nestimated population mean",
       color = NULL, linetype = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position   = "bottom",
        legend.key.width  = unit(1.8, "lines"),
        legend.text       = element_text(size = 11),
        panel.grid        = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_text(size = 11, face = "bold"))

ggsave("outputs/figures/gamma_bias_plot.png",
       p_gamma_bias, width = 9, height = 3.6, dpi = 150)
cat("Saved: outputs/figures/gamma_bias_plot.png\n")

# RMSE plot
p_gamma_rmse <- p_gamma_bias +
  aes(y = rmse) +
  scale_y_continuous(limits = c(0, ceiling(max(gd$rmse, na.rm = TRUE)))) +
  labs(y = "RMSE of\nestimated population mean")

ggsave("outputs/figures/gamma_rmse_plot.png",
       p_gamma_rmse, width = 9, height = 3.6, dpi = 150)
cat("Saved: outputs/figures/gamma_rmse_plot.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Gamma gradient – distribution (violin + boxplot)
# ══════════════════════════════════════════════════════════════════════════════

method_labels_short <- c(
  np_naive_mean             = "Naive (unadjusted)",
  calibrated_integration    = "Calibration",
  membership_ipw            = "IPW",
  doubly_robust_integration = "Doubly robust",
  prob_sample_mean          = "Probability sample (reference)"
)

raw_gamma <- read.csv("outputs/tables/gamma_raw_results.csv")
rd_gamma  <- subset(raw_gamma, method %in% METHODS_CORE)
rd_gamma$method <- factor(rd_gamma$method,
                           levels = METHODS_CORE,
                           labels = method_labels_short[METHODS_CORE])
rd_gamma$aux_quality <- factor(rd_gamma$aux_quality,
                                levels = c("strong", "weak"),
                                labels = c("Strong auxiliary information",
                                           "Weak auxiliary information"))
rd_gamma$gamma <- factor(rd_gamma$gamma)

true_means_gamma <- aggregate(true_mean ~ aux_quality, data = rd_gamma, FUN = mean)

p_gamma_dist <- ggplot(rd_gamma, aes(x = gamma, y = estimate,
                                      color = method, fill = method)) +
  geom_hline(data = true_means_gamma,
             aes(yintercept = true_mean),
             color = "grey30", linewidth = 0.5, linetype = "dashed") +
  geom_violin(alpha = 0.15, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.18, alpha = 0.7, linewidth = 0.4,
               outlier.size = 0.4, outlier.alpha = 0.4) +
  facet_grid(aux_quality ~ method, scales = "free_y") +
  scale_color_manual(values = COLORS, guide = "none") +
  scale_fill_manual(values  = COLORS, guide = "none") +
  labs(x = expression("Non-ignorability parameter " * gamma),
       y = "Estimated population mean") +
  theme_classic(base_size = 10) +
  theme(strip.background   = element_blank(),
        strip.text.x       = element_text(size = 8.5, face = "bold"),
        strip.text.y       = element_text(size = 9,   face = "bold"),
        panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
        axis.text.x        = element_text(size = 7.5))

ggsave("outputs/figures/gamma_distribution_plot.png",
       p_gamma_dist, width = 13, height = 5.5, dpi = 150)
cat("Saved: outputs/figures/gamma_distribution_plot.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Monte Carlo scenarios – RMSE dot (lollipop) plot
# ══════════════════════════════════════════════════════════════════════════════

mc_performance <- read.csv("outputs/tables/mc_performance.csv")

method_labels_mc <- c(
  np_naive_mean             = "Naive (unadjusted)",
  calibrated_integration    = "Calibration",
  membership_ipw            = "IPW",
  doubly_robust_integration = "Doubly robust",
  prob_sample_mean          = "Probability sample (reference)"
)

mc_plot <- subset(mc_performance, method %in% METHODS_CORE)
mc_plot$method <- factor(mc_plot$method,
                          levels = rev(METHODS_CORE),
                          labels = rev(method_labels_mc[METHODS_CORE]))
mc_plot$aux_quality <- factor(mc_plot$aux_quality,
                               levels = c("strong", "weak"),
                               labels = c("Strong auxiliary information",
                                          "Weak auxiliary information"))
mc_plot$selection_mechanism <- factor(mc_plot$selection_mechanism,
                                       levels = c("ignorable", "nonignorable"),
                                       labels = c("Ignorable selection",
                                                  "Non-ignorable selection"))

p_scenario_rmse <- ggplot(mc_plot, aes(x = rmse, y = method, color = method)) +
  geom_segment(aes(x = 0, xend = rmse, yend = method),
               linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 3) +
  facet_grid(selection_mechanism ~ aux_quality) +
  scale_color_manual(values = COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
  labs(x = "RMSE of estimated population mean", y = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
        strip.background    = element_blank(),
        strip.text          = element_text(size = 11, face = "bold"),
        axis.text.y         = element_text(size = 10))

ggsave("outputs/figures/scenario_rmse_plot.png",
       p_scenario_rmse, width = 9, height = 5, dpi = 150)
cat("Saved: outputs/figures/scenario_rmse_plot.png\n")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Monte Carlo scenarios – distribution (violin + boxplot)
# ══════════════════════════════════════════════════════════════════════════════

method_labels_violin <- c(
  np_naive_mean             = "Naive",
  calibrated_integration    = "Calibration / MI",
  membership_ipw            = "IPW",
  doubly_robust_integration = "Doubly robust",
  prob_sample_mean          = "Prob sample (ref)"
)

raw_mc  <- read.csv("outputs/tables/mc_raw_results.csv")
rd_mc   <- subset(raw_mc, method %in% METHODS_CORE)
rd_mc$method <- factor(rd_mc$method,
                        levels  = rev(METHODS_CORE),
                        labels  = rev(method_labels_violin[METHODS_CORE]))

label_colors <- COLORS[rev(method_labels_violin[METHODS_CORE])]
x_rng <- range(rd_mc$estimate, na.rm = TRUE)
x_lim <- x_rng + diff(x_rng) * c(-0.03, 0.03)

make_mc_panel <- function(aux, sel, title, show_y) {
  sub <- rd_mc[rd_mc$aux_quality == aux & rd_mc$selection_mechanism == sel, ]
  tm  <- mean(sub$true_mean, na.rm = TRUE)

  p <- ggplot(sub, aes(x = estimate, y = method,
                        color = method, fill = method)) +
    geom_vline(xintercept = tm, color = "grey30",
               linewidth = 0.6, linetype = "dashed") +
    geom_violin(alpha = 0.20, linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.18, alpha = 0.75, linewidth = 0.4,
                 outlier.size = 0.5, outlier.alpha = 0.4,
                 outlier.color = "grey60") +
    scale_color_manual(values = COLORS, guide = "none") +
    scale_fill_manual(values  = COLORS, guide = "none") +
    scale_x_continuous(limits = x_lim) +
    labs(x = "Estimated population mean", y = NULL, title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title         = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text.x        = element_text(size = 8),
          axis.title.x       = element_text(size = 8.5),
          panel.grid.major.x = element_line(color = "#eeeeee", linewidth = 0.4))

  if (show_y) {
    p + theme(axis.text.y = element_text(size = 8.5, color = label_colors))
  } else {
    p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
}

p_mc_violin <- (
  make_mc_panel("strong", "ignorable",    "Strong aux\nIgnorable",     TRUE)  |
  make_mc_panel("strong", "nonignorable", "Strong aux\nNon-ignorable", FALSE) |
  make_mc_panel("weak",   "ignorable",    "Weak aux\nIgnorable",       FALSE) |
  make_mc_panel("weak",   "nonignorable", "Weak aux\nNon-ignorable",   FALSE)
) +
  plot_annotation(
    title    = "Distribution of population estimates across 500 Monte Carlo simulations",
    subtitle = "N = 100,000   NP sample ≈3%   dashed line = true mean",
    theme    = theme(plot.title    = element_text(size = 11, hjust = 0.5),
                     plot.subtitle = element_text(size = 9,  hjust = 0.5))
  )

ggsave("outputs/figures/mc_violin_plot.png",
       p_mc_violin, width = 14, height = 5, dpi = 150)
cat("Saved: outputs/figures/mc_violin_plot.png\n")
