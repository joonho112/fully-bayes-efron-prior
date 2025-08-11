################################################################################
#
# Replication Code for:
# "Fully Bayesian Inference for Meta-Analytic Deconvolution 
#  Using Efron's Log-Spline Prior"
#
# Part 3: Figures and Tables for Publication
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: January 2025
#
# Description:
# This script generates all figures and tables presented in the manuscript,
# including density comparisons, performance metrics, coverage analysis,
# and diagnostic plots for assessing the validity of uncertainty quantification.
#
################################################################################

# ==============================================================================
# 1. SETUP AND PACKAGE LOADING
# ==============================================================================

# Clear workspace and memory
rm(list = ls(all.names = TRUE))
gc()

# Load required packages
required_packages <- c(
  "tidyverse", "ggplot2", "dplyr", "tidyr", "purrr",
  "scales", "patchwork", "viridis", "see"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Set working directory
# Users should modify these paths according to their system
work_dir <- "."  # Set your working directory here
data_dir <- "./data"  # Set your data directory here
fig_dir <- "./figures"  # Set your figures directory here

# Create directories if they don't exist
for (dir in c(data_dir, fig_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

setwd(work_dir)

# ==============================================================================
# 2. LOAD PROCESSED DATA
# ==============================================================================

cat("Loading processed simulation results...\n")

# Load the collected theta estimates from all methods
# This file should be created by running Part 2
df <- readRDS(file.path(data_dir, "collected_theta_results.rds"))

# Display data structure
cat("\nData structure:\n")
cat("  Reliability levels (I):", unique(df$I), "\n")
cat("  Heterogeneity levels (R):", unique(df$R), "\n")
cat("  Estimators:", unique(df$estimator), "\n")
cat("  Parameters:", unique(df$parameter), "\n")
cat("  Total observations:", nrow(df), "\n\n")

# ==============================================================================
# 3. DEFINE COLOR SCHEMES AND THEMES
# ==============================================================================

# Define consistent color palette for estimators
estimator_colors <- c(
  "MCMC" = "#1f77b4",        # Blue
  "EB" = "#ff7f0e",          # Orange
  "Laplace" = "#2ca02c",     # Green
  "Variational" = "#d62728", # Red
  "Pathfinder" = "#9467bd",  # Purple
  "Penalized MLE" = "#8c564b", # Brown
  "MAP" = "#e377c2"          # Pink
)

# Set default theme
theme_set(theme_bw(base_size = 11))

# ==============================================================================
# 4. FIGURE 3: DENSITY COMPARISONS
# ==============================================================================

cat("Creating Figure 3: Density comparisons...\n")

#' Plot Density Comparison
#'
#' This function creates faceted density plots comparing true, observed,
#' and estimated distributions across simulation conditions.
#'
#' @param data Input data frame with simulation results
#' @param I_subset Vector of reliability values to include
#' @param R_subset Vector of heterogeneity values to include
#' @param estimator_subset Vector of estimators to compare
#' @param parameter_subset Parameter type to visualize
#' @return ggplot object

plot_density_comparison <- function(data, 
                                    I_subset = c(0.5, 0.7, 0.9),
                                    R_subset = c(1, 9),
                                    estimator_subset = c("MCMC", "EB"),
                                    parameter_subset = "theta_mean") {
  
  # Filter data
  df_filtered <- data %>%
    filter(
      I %in% I_subset,
      R %in% R_subset,
      estimator %in% estimator_subset,
      parameter %in% parameter_subset
    )
  
  # Prepare data for plotting
  df_plot <- df_filtered %>%
    select(I, R, estimator, site_id, theta_true, theta_hat, mean) %>%
    pivot_longer(
      cols = c(theta_true, theta_hat, mean),
      names_to = "variable_type",
      values_to = "value"
    ) %>%
    filter(!is.na(value))
  
  # Create facet labels
  df_plot <- df_plot %>%
    mutate(
      facet_label = paste0("I = ", I, ", R = ", R),
      hetero_lab = if_else(R == 1, "Homoscedastic", "Heteroscedastic (R = 9)")
    )
  
  # Build plot
  p <- ggplot(df_plot, aes(x = value)) +
    # Histograms for true and observed distributions
    geom_histogram(
      data = filter(df_plot, variable_type == "theta_true"),
      aes(y = after_stat(density), fill = "True"),
      bins = 30, alpha = 0.3, color = NA
    ) +
    geom_histogram(
      data = filter(df_plot, variable_type == "theta_hat"),
      aes(y = after_stat(density), fill = "Observed"),
      bins = 30, alpha = 0.3, color = NA
    ) +
    # Density curves
    geom_density(
      data = filter(df_plot, variable_type == "theta_true"),
      aes(color = "True"),
      size = 1.2
    ) +
    geom_density(
      data = filter(df_plot, variable_type == "theta_hat"),
      aes(color = "Observed"),
      size = 1.2
    ) +
    geom_density(
      data = filter(df_plot, variable_type == "mean"),
      aes(color = estimator),
      size = 1.2
    ) +
    # Faceting
    facet_grid(
      rows = vars(I),
      cols = vars(hetero_lab),
      labeller = labeller(I = as_labeller(function(x) paste("I =", x)))
    ) +
    # Scales
    scale_fill_manual(
      name = "Distribution",
      values = c("True" = "#d62728", "Observed" = "gray60"),
      labels = c("True" = expression(theta[true]), 
                 "Observed" = expression(hat(theta)))
    ) +
    scale_color_manual(
      name = "Distribution",
      values = c("True" = "#d62728", "Observed" = "gray40",
                 "MCMC" = estimator_colors["MCMC"],
                 "EB" = estimator_colors["EB"])
    ) +
    # Labels
    labs(
      x = expression(theta),
      y = "Density",
      title = "Figure 3. Distribution of True, Observed, and Estimated Site-Specific Effects"
    ) +
    # Theme
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    coord_cartesian(xlim = c(-5, 5))
  
  return(p)
}

# Create Figure 3
fig3 <- plot_density_comparison(
  df,
  estimator_subset = c("MCMC", "EB"),
  parameter_subset = "theta_mean"
)

# Save Figure 3
ggsave(
  filename = file.path(fig_dir, "Figure_03_Density_Comparisons.pdf"),
  plot = fig3,
  width = 10,
  height = 10,
  dpi = 300
)

cat("  Figure 3 saved.\n")

# ==============================================================================
# 5. FIGURE 4: PERFORMANCE METRICS
# ==============================================================================

cat("Creating Figure 4: Performance metrics...\n")

#' Calculate Performance Metrics
#'
#' This function computes RMSE, correlation, and other metrics for each
#' estimator under different simulation conditions.
#'
#' @param data Input data frame
#' @param estimator_subset Estimators to evaluate
#' @param parameter_subset Parameter type to evaluate
#' @return Data frame with performance metrics

calculate_performance_metrics <- function(data,
                                          I_subset = c(0.5, 0.7, 0.9),
                                          R_subset = c(1, 9),
                                          estimator_subset = c("MCMC", "EB"),
                                          parameter_subset = "theta_mean") {
  
  metrics <- data %>%
    filter(
      I %in% I_subset,
      R %in% R_subset,
      estimator %in% estimator_subset,
      parameter %in% parameter_subset,
      !is.na(mean)
    ) %>%
    group_by(I, R, estimator, parameter) %>%
    summarise(
      n = n(),
      RMSE = sqrt(mean((mean - theta_true)^2)),
      MAE = mean(abs(mean - theta_true)),
      Bias = mean(mean - theta_true),
      Correlation = cor(mean, theta_true),
      .groups = "drop"
    )
  
  return(metrics)
}

#' Plot Performance Metrics
#'
#' Creates bar plots of performance metrics across conditions

plot_performance_metrics <- function(metrics, metric = "RMSE") {
  
  p <- ggplot(metrics, aes(x = factor(R), y = .data[[metric]], 
                           fill = estimator, group = estimator)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             alpha = 0.8) +
    geom_text(
      aes(label = sprintf("%.3f", round(.data[[metric]], 3))),
      position = position_dodge(width = 0.8),
      vjust = -0.5, size = 3
    ) +
    facet_grid(
      parameter ~ I,
      scales = "free_y",
      labeller = labeller(I = function(x) paste("I =", x))
    ) +
    scale_fill_manual(values = estimator_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "R (Heterogeneity Ratio)",
      y = metric,
      fill = "Estimator",
      title = paste("Figure 4.", metric, "by Simulation Conditions")
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}

# Calculate metrics
metrics <- calculate_performance_metrics(
  df,
  estimator_subset = c("MCMC", "EB", "MAP"),
  parameter_subset = "theta_mean"
)

# Create performance plots
fig4_rmse <- plot_performance_metrics(metrics, metric = "RMSE")
fig4_corr <- plot_performance_metrics(metrics, metric = "Correlation")

# Combine plots
fig4 <- fig4_rmse / fig4_corr

# Save Figure 4
ggsave(
  filename = file.path(fig_dir, "Figure_04_Performance_Metrics.pdf"),
  plot = fig4,
  width = 9,
  height = 10,
  dpi = 300
)

cat("  Figure 4 saved.\n")

# ==============================================================================
# 6. FIGURE 5: 90% CREDIBLE INTERVALS FOR RANDOM SITES
# ==============================================================================

cat("Creating Figure 5: Credible intervals for random sites...\n")

#' Plot Credible Intervals
#'
#' Shows 90% credible intervals for a random sample of sites

plot_credible_intervals <- function(data, 
                                    condition_I, 
                                    condition_R,
                                    n_sites = 150,
                                    seed = 123) {
  
  set.seed(seed)
  
  # Filter and sample data
  df_filtered <- data %>%
    filter(
      I == condition_I,
      R == condition_R,
      estimator == "MCMC",
      parameter == "theta_rep"
    ) %>%
    sample_n(n_sites) %>%
    arrange(theta_true) %>%
    mutate(
      site_index = row_number(),
      covers_true = theta_true >= q5 & theta_true <= q95,
      color = ifelse(covers_true, "black", "blue")
    )
  
  # Calculate coverage rate
  coverage_rate <- mean(df_filtered$covers_true)
  
  # Create plot
  p <- ggplot(df_filtered, aes(x = site_index)) +
    # Error bars
    geom_errorbar(
      aes(ymin = q5, ymax = q95, color = color),
      width = 0, size = 0.5
    ) +
    # Point estimates
    geom_point(
      aes(y = mean, color = color),
      size = 1.5, shape = 16
    ) +
    # True values
    geom_point(
      aes(y = theta_true),
      color = "red", size = 2, shape = 17
    ) +
    # Reference line
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    # Scales
    scale_color_identity() +
    # Labels
    labs(
      x = "Site Index (sorted by true effect)",
      y = expression(theta),
      title = sprintf("90%% Credible Intervals for %d Random Sites", n_sites),
      subtitle = sprintf("I = %.1f, R = %d | Coverage: %.1f%% | Red triangles: true values",
                         condition_I, condition_R, coverage_rate * 100)
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  return(p)
}

# Create plots for two conditions
fig5_a <- plot_credible_intervals(df, condition_I = 0.5, condition_R = 9)
fig5_b <- plot_credible_intervals(df, condition_I = 0.7, condition_R = 9)

# Combine plots
fig5 <- fig5_a / fig5_b

# Save Figure 5
ggsave(
  filename = file.path(fig_dir, "Figure_05_Credible_Intervals.pdf"),
  plot = fig5,
  width = 10,
  height = 10,
  dpi = 300
)

cat("  Figure 5 saved.\n")

# ==============================================================================
# 7. FIGURE 6: COVERAGE RATES
# ==============================================================================

cat("Creating Figure 6: Coverage rates...\n")

#' Calculate Coverage Rates
#'
#' Computes empirical coverage of credible/confidence intervals

calculate_coverage <- function(data,
                               estimator_subset = c("MCMC", "EB"),
                               parameter_subset = c("theta_rep", "theta_mean"),
                               confidence_level = 90) {
  
  # Determine quantile columns
  lower_q <- paste0("q", (100 - confidence_level) / 2)
  upper_q <- paste0("q", 100 - (100 - confidence_level) / 2)
  
  coverage_rates <- data %>%
    filter(
      estimator %in% estimator_subset,
      parameter %in% parameter_subset
    ) %>%
    group_by(I, R, estimator, parameter) %>%
    summarise(
      n_total = n(),
      n_covered = sum(theta_true >= .data[[lower_q]] & 
                        theta_true <= .data[[upper_q]], na.rm = TRUE),
      coverage_rate = n_covered / n_total,
      se_coverage = sqrt(coverage_rate * (1 - coverage_rate) / n_total),
      .groups = "drop"
    ) %>%
    mutate(nominal_level = confidence_level / 100)
  
  return(coverage_rates)
}

#' Plot Coverage Rates

plot_coverage_rates <- function(coverage_data) {
  
  p <- ggplot(coverage_data, aes(x = factor(R), y = coverage_rate,
                                 fill = estimator)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             alpha = 0.8) +
    geom_hline(aes(yintercept = nominal_level),
               linetype = "dashed", color = "red") +
    geom_text(
      aes(label = sprintf("%.1f%%", coverage_rate * 100)),
      position = position_dodge(width = 0.8),
      vjust = -0.5, size = 3
    ) +
    facet_grid(
      parameter ~ I,
      labeller = labeller(I = function(x) paste("I =", x))
    ) +
    scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.1))
    ) +
    scale_fill_manual(values = estimator_colors) +
    labs(
      x = "R (Heterogeneity Ratio)",
      y = "Coverage Rate",
      fill = "Estimator",
      title = "Figure 6. 90% Credible/Confidence Interval Coverage Rates",
      subtitle = "Red dashed line indicates nominal 90% coverage"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}

# Calculate coverage
coverage_results <- calculate_coverage(
  df,
  estimator_subset = c("MCMC", "EB"),
  parameter_subset = c("theta_rep", "theta_mean")
)

# Create plot
fig6 <- plot_coverage_rates(coverage_results)

# Save Figure 6
ggsave(
  filename = file.path(fig_dir, "Figure_06_Coverage_Rates.pdf"),
  plot = fig6,
  width = 8,
  height = 7,
  dpi = 300
)

cat("  Figure 6 saved.\n")

# ==============================================================================
# 8. FIGURE 7: CALIBRATION PLOT
# ==============================================================================

cat("Creating Figure 7: Calibration plot...\n")

#' Plot Calibration
#'
#' Creates calibration plots comparing nominal vs empirical coverage

plot_calibration <- function(data,
                             estimator_subset = c("MCMC", "EB"),
                             parameter = "theta_rep") {
  
  # Calculate empirical coverage at various nominal levels
  nominal_levels <- seq(0.1, 0.95, by = 0.05)
  
  calibration_data <- map_dfr(nominal_levels, function(level) {
    z_score <- qnorm((1 + level) / 2)
    
    data %>%
      filter(
        estimator %in% estimator_subset,
        parameter == !!parameter,
        !is.na(sd)
      ) %>%
      group_by(I, R, estimator) %>%
      summarise(
        nominal = level,
        empirical = mean(abs((theta_true - mean) / sd) <= z_score),
        .groups = "drop"
      )
  })
  
  # Create plot
  p <- ggplot(calibration_data, aes(x = nominal, y = empirical,
                                    color = estimator, shape = estimator)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "gray50") +
    geom_line(size = 1) +
    geom_point(size = 2.2) +
    facet_grid(
      R ~ I,
      labeller = labeller(
        I = function(x) paste("I =", x),
        R = function(x) paste("R =", x)
      )
    ) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = estimator_colors) +
    coord_equal() +
    labs(
      x = "Nominal Coverage",
      y = "Empirical Coverage",
      color = "Estimator",
      shape = "Estimator",
      title = "Figure 7. Calibration Plot: Nominal vs Empirical Coverage",
      subtitle = "Perfect calibration follows the diagonal line"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Create calibration plot
fig7 <- plot_calibration(df, estimator_subset = c("MCMC", "EB"))

# Save Figure 7
ggsave(
  filename = file.path(fig_dir, "Figure_07_Calibration_Plot.pdf"),
  plot = fig7,
  width = 7,
  height = 7,
  dpi = 300
)

cat("  Figure 7 saved.\n")

# ==============================================================================
# 9. FIGURE 8: SEQUENTIAL COVERAGE
# ==============================================================================

cat("Creating Figure 8: Sequential coverage plot...\n")

#' Plot Sequential Coverage
#'
#' Shows how coverage rates evolve as observations accumulate

plot_sequential_coverage <- function(data,
                                     estimator_subset = c("MCMC"),
                                     parameter = "theta_rep") {
  
  seq_coverage <- data %>%
    filter(
      estimator %in% estimator_subset,
      parameter == !!parameter
    ) %>%
    group_by(I, R, estimator) %>%
    arrange(site_id) %>%
    mutate(
      covers = theta_true >= q5 & theta_true <= q95,
      cumulative_coverage = cummean(covers),
      observation_number = row_number()
    ) %>%
    filter(observation_number %% 10 == 0)  # Thin for visualization
  
  p <- ggplot(seq_coverage, aes(x = observation_number, y = cumulative_coverage,
                                color = estimator)) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    geom_line(size = 1) +
    facet_wrap(~ paste("I =", I, ", R =", R), scales = "free_x") +
    scale_y_continuous(labels = scales::percent, limits = c(0.7, 1)) +
    scale_color_manual(values = estimator_colors) +
    labs(
      x = "Number of Observations",
      y = "Cumulative Coverage Rate",
      color = "Estimator",
      title = "Figure 8. Sequential Coverage Rate",
      subtitle = "Should converge to 90% (red line) as sample size increases"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}

# Create sequential coverage plot
fig8 <- plot_sequential_coverage(df, estimator_subset = c("MCMC"))

# Save Figure 8
ggsave(
  filename = file.path(fig_dir, "Figure_08_Sequential_Coverage.pdf"),
  plot = fig8,
  width = 8,
  height = 7,
  dpi = 300
)

cat("  Figure 8 saved.\n")

# ==============================================================================
# 10. FIGURE 9: ADDITIONAL VALIDITY CHECKS
# ==============================================================================

cat("Creating Figure 9: Additional validity checks...\n")

#' PIT Histogram
#'
#' Probability Integral Transform histogram for calibration assessment

plot_pit_histogram <- function(data,
                               estimator_subset = c("MCMC"),
                               parameter = "theta_rep") {
  
  pit_data <- data %>%
    filter(
      estimator %in% estimator_subset,
      parameter == !!parameter,
      !is.na(mean) & !is.na(sd)
    ) %>%
    mutate(pit = pnorm(theta_true, mean = mean, sd = sd))
  
  p <- ggplot(pit_data, aes(x = pit)) +
    geom_histogram(
      bins = 20, fill = "skyblue", color = "black",
      alpha = 0.7, aes(y = after_stat(density))
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
    facet_grid(estimator ~ paste("I =", I, ", R =", R)) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      x = "PIT Value",
      y = "Density",
      title = "Probability Integral Transform (PIT) Histogram"
    ) +
    theme_bw()
  
  return(p)
}

#' Q-Q Plot of Standardized Residuals

plot_zscore_qq <- function(data,
                           estimator_subset = c("MCMC"),
                           parameter = "theta_rep") {
  
  qq_data <- data %>%
    filter(
      estimator %in% estimator_subset,
      parameter == !!parameter,
      !is.na(sd) & sd > 0
    ) %>%
    mutate(z_score = (theta_true - mean) / sd) %>%
    group_by(I, R, estimator) %>%
    arrange(z_score) %>%
    mutate(
      theoretical_quantiles = qnorm(ppoints(n())),
      empirical_quantiles = z_score
    )
  
  p <- ggplot(qq_data, aes(x = theoretical_quantiles, y = empirical_quantiles)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_point(alpha = 0.5, size = 0.8) +
    facet_grid(estimator ~ paste("I =", I, ", R =", R)) +
    labs(
      x = "Theoretical Quantiles",
      y = "Empirical Quantiles",
      title = "Q-Q Plot of Standardized Residuals"
    ) +
    theme_bw()
  
  return(p)
}

#' Interval Width vs Error Plot

plot_width_error_tradeoff <- function(data,
                                      estimator_subset = c("MCMC"),
                                      parameter = "theta_rep",
                                      sample_size = 500) {
  
  set.seed(123)
  
  tradeoff_data <- data %>%
    filter(
      estimator %in% estimator_subset,
      parameter == !!parameter
    ) %>%
    group_by(I, R, estimator) %>%
    sample_n(min(n(), sample_size)) %>%
    ungroup() %>%
    mutate(
      interval_width = q95 - q5,
      covers = theta_true >= q5 & theta_true <= q95,
      abs_error = abs(theta_true - mean)
    )
  
  p <- ggplot(tradeoff_data, aes(x = interval_width, y = abs_error)) +
    geom_point(aes(color = covers), alpha = 0.6, size = 1.5) +
    geom_smooth(aes(group = covers, color = covers),
                method = "loess", se = FALSE) +
    facet_grid(estimator ~ paste("I =", I, ", R =", R)) +
    scale_color_manual(
      values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
      labels = c("TRUE" = "Covers", "FALSE" = "Misses")
    ) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
      x = "90% Interval Width (log scale)",
      y = "Absolute Error (log scale)",
      color = "Coverage",
      title = "Interval Width vs. Estimation Error"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Create validity check plots
fig9_pit <- plot_pit_histogram(df)
fig9_qq <- plot_zscore_qq(df)
fig9_tradeoff <- plot_width_error_tradeoff(df)

# Combine plots
fig9 <- fig9_pit / fig9_qq / fig9_tradeoff

# Save Figure 9
ggsave(
  filename = file.path(fig_dir, "Figure_09_Validity_Checks.pdf"),
  plot = fig9,
  width = 12,
  height = 10,
  dpi = 300
)

cat("  Figure 9 saved.\n")

# ==============================================================================
# 11. CREATE SUMMARY TABLES
# ==============================================================================

cat("\nCreating summary tables...\n")

# Table 1: Performance metrics summary
table1 <- metrics %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  arrange(I, R, estimator) %>%
  select(I, R, estimator, RMSE, MAE, Bias, Correlation)

cat("\nTable 1: Performance Metrics Summary\n")
print(table1)

# Save as CSV
write.csv(
  table1,
  file = file.path(data_dir, "Table_1_Performance_Metrics.csv"),
  row.names = FALSE
)

# Table 2: Coverage rates summary
table2 <- coverage_results %>%
  mutate(
    Coverage = sprintf("%.1f%%", coverage_rate * 100),
    SE = sprintf("(%.1f%%)", se_coverage * 100)
  ) %>%
  select(I, R, estimator, parameter, Coverage, SE)

cat("\nTable 2: Coverage Rates Summary\n")
print(table2)

# Save as CSV
write.csv(
  table2,
  file = file.path(data_dir, "Table_2_Coverage_Rates.csv"),
  row.names = FALSE
)

# ==============================================================================
# 12. SESSION INFORMATION
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("All figures and tables generated successfully!\n")
cat("\nOutput files:\n")
cat("  Figures:\n")
for (i in 3:9) {
  cat(sprintf("    - Figure_%02d_*.pdf\n", i))
}
cat("  Tables:\n")
cat("    - Table_1_Performance_Metrics.csv\n")
cat("    - Table_2_Coverage_Rates.csv\n")
cat("\nSession Information:\n")
print(sessionInfo())