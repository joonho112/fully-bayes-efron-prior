################################################################################
#
# Replication Code for:
# "Fully Bayesian Inference for Meta-Analytic Deconvolution 
#  Using Efron's Log-Spline Prior"
#
# Appendix C: Grid Resolution and Bounds Sensitivity Analysis
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: January 2025
#
# Description:
# This script evaluates the sensitivity of our fully Bayesian inference
# framework to the discretization scheme (grid resolution and bounds).
# We use Variational Bayes for computational efficiency across multiple
# configurations, demonstrating robustness to grid specification choices.
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
  "tidyverse", "cmdstanr", "splines", "patchwork", 
  "scales", "viridis", "knitr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "cmdstanr") {
      install.packages("cmdstanr", 
                       repos = c("https://mc-stan.org/r-packages/", 
                                 getOption("repos")))
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Set working directory
# Users should modify these paths according to their system
work_dir <- "."  # Set your working directory here
data_dir <- "./data"  # Set your data directory here
stan_dir <- "./stan"  # Set your Stan files directory here
output_dir <- "./appendix_c"  # Output directory for sensitivity analysis

# Create directories if they don't exist
for (dir in c(data_dir, stan_dir, output_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

setwd(work_dir)

# ==============================================================================
# 2. LOAD TEST SCENARIO DATA
# ==============================================================================

cat("Loading simulated data for sensitivity analysis...\n")

# Load the pre-simulated dataset
df_sim <- readRDS(file.path(data_dir, "simulated_data_varied_by_I_and_R.rds"))

# Select a challenging scenario for analysis: I=0.7, R=9
# This represents moderate reliability with heteroscedastic measurement errors
df_test <- df_sim %>%
  filter(I == 0.7, R == 9) %>%
  slice(1) %>%
  pull(data) %>%
  .[[1]]

# Extract data components
K <- nrow(df_test)
theta_true <- df_test$theta_true
theta_hat <- df_test$theta_hat
sigma <- df_test$se

cat(sprintf("Test scenario: K=%d sites, I=0.7, R=9\n", K))
cat(sprintf("True effect range: [%.2f, %.2f]\n", 
            min(theta_true), max(theta_true)))

# ==============================================================================
# 3. DEFINE GRID CONFIGURATIONS
# ==============================================================================

#' Prepare Stan Data with Different Grid Configurations
#'
#' This function creates Stan data with varying grid resolutions and bounds
#' to test sensitivity to discretization choices.
#'
#' @param theta_hat Observed effect estimates
#' @param sigma Standard errors
#' @param theta_true True effects (for bounds calculation)
#' @param L Number of grid points
#' @param M Degrees of freedom for splines (fixed at 6)
#' @param bound_expansion Factor to expand grid bounds beyond data range
#' @return List containing Stan data

prepare_stan_data_grid <- function(theta_hat, sigma, theta_true,
                                   L = 101, M = 6, 
                                   bound_expansion = 0.5) {
  
  # Calculate grid bounds based on true values
  range_true <- range(theta_true)
  range_expansion <- diff(range_true) * bound_expansion
  
  # Create evenly spaced grid
  grid <- seq(
    from = range_true[1] - range_expansion,
    to = range_true[2] + range_expansion,
    length.out = L
  )
  
  # Create natural cubic spline basis matrix
  B <- splines::ns(grid, df = M, intercept = FALSE)
  
  # Return Stan data list
  list(
    K = length(theta_hat),
    theta_hat = theta_hat,
    sigma = sigma,
    L = L,
    grid = grid,
    M = M,
    B = B
  )
}

# Define grid configurations to test (3x3 factorial design)
grid_configs <- expand_grid(
  L = c(51, 101, 201),                    # Grid resolution: coarse, moderate, fine
  bound_expansion = c(0.25, 0.5, 1.0)     # Bounds: conservative, moderate, generous
) %>%
  mutate(
    config_id = row_number(),
    config_label = paste0("L", L, "_bound", bound_expansion)
  )

cat("\nGrid configurations to test:\n")
print(grid_configs)

# ==============================================================================
# 4. COMPILE STAN MODEL
# ==============================================================================

cat("\nCompiling Stan model...\n")

# Check if Stan model file exists
stan_file <- file.path(stan_dir, "efron_log_spline_model.stan")

if (!file.exists(stan_file)) {
  stop("Stan model file not found. Please ensure 'efron_log_spline_model.stan' ",
       "is in the stan directory.")
}

# Compile model without threading (required for VB)
compiled_model <- cmdstan_model(
  stan_file,
  cpp_options = list(stan_threads = FALSE),
  stanc_options = list("O1")
)

cat("Model compiled successfully.\n")

# ==============================================================================
# 5. RUN VARIATIONAL BAYES FOR EACH CONFIGURATION
# ==============================================================================

#' Fit Model Using Variational Bayes and Extract Results
#'
#' This function runs Variational Bayes (ADVI) for a given configuration
#' and extracts point estimates and credible intervals.
#'
#' @param stan_data Prepared Stan data
#' @param compiled_model Compiled Stan model
#' @param config_id Configuration identifier
#' @param seed Random seed for reproducibility
#' @return Tibble with results for all sites

fit_vb_and_extract <- function(stan_data, compiled_model, 
                               config_id = 1, seed = 123) {
  
  cat(sprintf("Fitting configuration %d (L=%d, bound expansion=%.2f)...\n", 
              config_id, stan_data$L, 
              stan_data$B[nrow(stan_data$B), ncol(stan_data$B)]))
  
  tryCatch({
    # Start timing
    start_time <- Sys.time()
    
    # Run Variational Bayes with ADVI algorithm
    fit_vb <- compiled_model$variational(
      data = stan_data,
      seed = seed + config_id,
      algorithm = "meanfield",    # Mean-field approximation
      iter = 10000,               # Number of iterations
      grad_samples = 1,           # Gradient samples per iteration
      elbo_samples = 100,         # ELBO samples for convergence check
      eta = 0.1,                  # Step size
      adapt_iter = 50,            # Adaptation iterations
      output_samples = 4000,      # Posterior samples to draw
      refresh = 0                 # Suppress iteration output
    )
    
    # Calculate elapsed time
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    
    # Extract posterior summaries
    theta_mean_summary <- fit_vb$summary("theta_mean")
    theta_rep_summary <- fit_vb$summary("theta_rep")
    
    # Get posterior samples for credible intervals
    theta_rep_draws <- fit_vb$draws("theta_rep", format = "matrix")
    
    # Calculate 90% credible intervals
    ci_90 <- apply(theta_rep_draws, 2, quantile, probs = c(0.05, 0.95))
    
    # Compile results
    results <- tibble(
      config_id = config_id,
      site_id = 1:K,
      theta_true = theta_true,
      theta_hat = theta_hat,
      sigma = sigma,
      # Point estimates
      theta_mean_est = theta_mean_summary$mean,
      theta_rep_est = theta_rep_summary$mean,
      # Uncertainty quantification
      theta_rep_sd = theta_rep_summary$sd,
      ci_lower = ci_90[1,],
      ci_upper = ci_90[2,],
      # Computational cost
      elapsed_time = elapsed_time
    )
    
    cat(sprintf("  Completed in %.2f seconds\n", elapsed_time))
    
    return(results)
    
  }, error = function(e) {
    warning(sprintf("Failed config %d: %s\n", config_id, e$message))
    return(NULL)
  })
}

cat("\nRunning Variational Bayes for all configurations...\n")

# Run VB for all configurations
results_list <- grid_configs %>%
  mutate(
    stan_data = pmap(
      list(L, bound_expansion),
      ~ prepare_stan_data_grid(
        theta_hat = theta_hat,
        sigma = sigma,
        theta_true = theta_true,
        L = ..1,
        bound_expansion = ..2
      )
    ),
    vb_results = map2(
      stan_data, config_id,
      ~ fit_vb_and_extract(.x, compiled_model, .y)
    )
  )

# Combine all results
results_all <- results_list %>%
  select(L, bound_expansion, config_label, vb_results) %>%
  unnest(vb_results) %>%
  filter(!is.na(theta_mean_est))  # Remove any failed configurations

cat("Variational Bayes fitting completed for all configurations.\n")

# ==============================================================================
# 6. CALCULATE PERFORMANCE METRICS
# ==============================================================================

cat("\nCalculating performance metrics...\n")

# Calculate comprehensive performance metrics for each configuration
performance_metrics <- results_all %>%
  group_by(config_id, L, bound_expansion, config_label) %>%
  summarise(
    # Point estimation accuracy
    rmse_mean = sqrt(mean((theta_mean_est - theta_true)^2)),
    rmse_rep = sqrt(mean((theta_rep_est - theta_true)^2)),
    correlation = cor(theta_mean_est, theta_true),
    
    # Uncertainty calibration
    coverage_90 = mean(theta_true >= ci_lower & theta_true <= ci_upper),
    
    # Interval properties
    mean_interval_width = mean(ci_upper - ci_lower),
    
    # Computational efficiency
    computation_time = mean(elapsed_time),
    
    .groups = "drop"
  ) %>%
  # Add reference configuration indicator (our default: L=101, bound=0.5)
  mutate(
    is_reference = (L == 101 & bound_expansion == 0.5)
  )

# Print summary
cat("\nPerformance Metrics Summary:\n")
print(performance_metrics %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# ==============================================================================
# 7. CREATE VISUALIZATIONS
# ==============================================================================

cat("\nCreating sensitivity analysis plots...\n")

# Set consistent theme
theme_custom <- theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "bottom"
  )

# Panel A: Point Estimation Accuracy (RMSE)
p1_rmse <- performance_metrics %>%
  pivot_longer(
    cols = c(rmse_mean, rmse_rep),
    names_to = "estimator",
    values_to = "rmse",
    names_prefix = "rmse_"
  ) %>%
  mutate(
    estimator = factor(estimator, 
                       levels = c("mean", "rep"),
                       labels = c("Posterior Mean", "Posterior Draw"))
  ) %>%
  ggplot(aes(x = factor(L), y = rmse, color = factor(bound_expansion))) +
  geom_point(aes(shape = is_reference), size = 3.5) +
  geom_line(aes(group = bound_expansion), alpha = 0.6) +
  facet_wrap(~ estimator) +
  scale_color_viridis_d(name = "Bound Expansion Factor") +
  scale_shape_manual(values = c(16, 17), 
                     labels = c("Alternative", "Reference (L=101, β=0.5)"),
                     name = "") +
  labs(
    x = "Number of Grid Points (L)",
    y = "Root Mean Squared Error",
    title = "A. Point Estimation Accuracy"
  ) +
  theme_custom

# Panel B: Coverage Probability
p2_coverage <- performance_metrics %>%
  ggplot(aes(x = factor(L), y = coverage_90, 
             color = factor(bound_expansion))) +
  geom_point(aes(shape = is_reference), size = 3.5) +
  geom_line(aes(group = bound_expansion), alpha = 0.6) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "red", 
             alpha = 0.7) +
  scale_color_viridis_d(name = "Bound Expansion Factor") +
  scale_shape_manual(values = c(16, 17), 
                     labels = c("Alternative", "Reference"),
                     name = "") +
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0.85, 0.95)) +
  labs(
    x = "Number of Grid Points (L)",
    y = "90% Coverage Probability",
    title = "B. Uncertainty Calibration"
  ) +
  theme_custom

# Panel C: Computational Cost
p3_time <- performance_metrics %>%
  ggplot(aes(x = factor(L), y = computation_time, 
             fill = factor(bound_expansion))) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  scale_fill_viridis_d(name = "Bound Expansion Factor") +
  labs(
    x = "Number of Grid Points (L)",
    y = "Computation Time (seconds)",
    title = "C. Computational Cost"
  ) +
  theme_custom

# Panel D: Accuracy-Coverage Trade-off
p4_tradeoff <- performance_metrics %>%
  ggplot(aes(x = rmse_mean, y = coverage_90)) +
  geom_point(aes(size = factor(L), color = factor(bound_expansion),
                 shape = is_reference), alpha = 0.8) +
  geom_hline(yintercept = 0.90, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  scale_size_manual(values = c(3, 5, 7), name = "Grid Points") +
  scale_color_viridis_d(name = "Bound Expansion") +
  scale_shape_manual(values = c(16, 17), 
                     labels = c("Alternative", "Reference"),
                     name = "") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "RMSE (Posterior Mean)",
    y = "90% Coverage Probability",
    title = "D. Accuracy-Coverage Trade-off"
  ) +
  theme_custom

# Combine all panels
combined_plot <- (p1_rmse | p2_coverage) / (p3_time | p4_tradeoff) +
  plot_annotation(
    title = "Grid Sensitivity Analysis: Impact of Resolution and Bounds on Inference Quality",
    subtitle = "Twin towers scenario (K=1,500, I=0.7, R=9) using Variational Bayes",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save the main figure
ggsave(
  filename = file.path(output_dir, "Figure_C1_Grid_Sensitivity.pdf"),
  plot = combined_plot,
  width = 14,
  height = 10,
  dpi = 300
)

cat("  Figure C1 saved.\n")

# ==============================================================================
# 8. CREATE SUMMARY TABLE
# ==============================================================================

cat("\nCreating summary comparison table...\n")

# Create detailed comparison table
comparison_table <- performance_metrics %>%
  # Get reference values
  mutate(
    ref_rmse = rmse_mean[is_reference],
    ref_coverage = coverage_90[is_reference],
    ref_time = computation_time[is_reference]
  ) %>%
  # Calculate relative changes
  mutate(
    rmse_change_pct = (rmse_mean - ref_rmse) / ref_rmse * 100,
    coverage_diff_pp = (coverage_90 - ref_coverage) * 100,  # percentage points
    time_ratio = computation_time / ref_time
  ) %>%
  select(
    L, bound_expansion, is_reference,
    rmse_mean, rmse_change_pct,
    coverage_90, coverage_diff_pp,
    computation_time, time_ratio
  ) %>%
  arrange(L, bound_expansion)

# Print formatted table
cat("\n========== GRID SENSITIVITY ANALYSIS RESULTS ==========\n")
cat("\nReference Configuration: L=101, Bound Expansion=0.5\n")
cat("Test Scenario: Twin Towers, I=0.7, R=9\n\n")

print(comparison_table %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>%
        knitr::kable(format = "simple"))

# Save results as CSV
write.csv(
  performance_metrics,
  file = file.path(output_dir, "Table_C1_Grid_Sensitivity_Metrics.csv"),
  row.names = FALSE
)

write.csv(
  comparison_table,
  file = file.path(output_dir, "Table_C2_Grid_Sensitivity_Comparison.csv"),
  row.names = FALSE
)

# ==============================================================================
# 9. ADDITIONAL DIAGNOSTIC PLOTS
# ==============================================================================

cat("\nCreating diagnostic plots...\n")

# Select three representative configurations for detailed comparison
selected_configs <- c(1, 5, 9)  # L=51/0.25, L=101/0.5 (reference), L=201/1.0

# Site-level estimation errors
p_site_errors <- results_all %>%
  filter(config_id %in% selected_configs) %>%
  mutate(
    error = theta_mean_est - theta_true,
    in_ci = theta_true >= ci_lower & theta_true <= ci_upper
  ) %>%
  ggplot(aes(x = theta_true, y = error, color = in_ci)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ config_label, scales = "free") +
  scale_color_manual(values = c("red", "darkgreen"),
                     labels = c("Outside CI", "Inside CI"),
                     name = "Coverage") +
  labs(
    x = "True Effect Size",
    y = "Estimation Error",
    title = "Figure C2. Site-level Estimation Errors Across Selected Grid Configurations",
    subtitle = "Comparison of coarse (L=51), reference (L=101), and fine (L=201) grids"
  ) +
  theme_custom

# Save diagnostic plot
ggsave(
  filename = file.path(output_dir, "Figure_C2_Site_Errors.pdf"),
  plot = p_site_errors,
  width = 12,
  height = 6,
  dpi = 300
)

cat("  Figure C2 saved.\n")

# ==============================================================================
# 10. SESSION INFORMATION
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Grid sensitivity analysis completed successfully!\n")
cat("\nKey findings:\n")
cat(sprintf("  - RMSE variation across configurations: %.1f%%\n",
            100 * (max(performance_metrics$rmse_mean) - 
                     min(performance_metrics$rmse_mean)) / 
              min(performance_metrics$rmse_mean)))
cat(sprintf("  - Coverage range: %.1f%% to %.1f%%\n",
            100 * min(performance_metrics$coverage_90),
            100 * max(performance_metrics$coverage_90)))
cat(sprintf("  - Computation time range: %.1f to %.1f seconds\n",
            min(performance_metrics$computation_time),
            max(performance_metrics$computation_time)))

cat("\nOutput files saved to:", output_dir, "\n")
cat("  - Figure_C1_Grid_Sensitivity.pdf\n")
cat("  - Figure_C2_Site_Errors.pdf\n")
cat("  - Table_C1_Grid_Sensitivity_Metrics.csv\n")
cat("  - Table_C2_Grid_Sensitivity_Comparison.csv\n")

cat("\nSession Information:\n")
print(sessionInfo())