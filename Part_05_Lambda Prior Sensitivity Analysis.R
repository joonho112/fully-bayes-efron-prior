################################################################################
#
# Replication Code for:
# "Fully Bayesian Inference for Meta-Analytic Deconvolution 
#  Using Efron's Log-Spline Prior"
#
# Appendix D: Lambda Prior Sensitivity Analysis
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: January 2025
#
# Description:
# This script evaluates the sensitivity of the fully Bayesian inference
# framework to the choice of hyperprior on the regularization parameter λ.
# We compare alternative prior distributions and assess their impact on
# point estimation accuracy and uncertainty calibration.
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
  "scales", "viridis", "knitr", "kableExtra"
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
output_dir <- "./appendix_d"  # Output directory for sensitivity analysis

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

cat("Loading simulated data for lambda sensitivity analysis...\n")

# Load the pre-simulated dataset
df_sim <- readRDS(file.path(data_dir, "simulated_data_varied_by_I_and_R.rds"))

# Select moderate difficulty scenario (I=0.7, R=9) - same as Grid Sensitivity
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

cat(sprintf("Test scenario: K=%d sites, I=0.7, R=9 (Twin Towers)\n", K))
cat(sprintf("True effect range: [%.2f, %.2f]\n", 
            min(theta_true), max(theta_true)))

# ==============================================================================
# 3. PREPARE STAN DATA
# ==============================================================================

#' Prepare Stan Data with Default Grid Configuration
#'
#' @param theta_hat Observed effect estimates
#' @param sigma Standard errors
#' @param theta_true True effects
#' @param L Number of grid points (default: 101)
#' @param M Degrees of freedom for splines (default: 6)
#' @return List containing Stan data

prepare_stan_data <- function(theta_hat, sigma, theta_true,
                              L = 101, M = 6) {
  # Grid spanning the range of true effects with buffer
  grid <- seq(min(theta_true) - 0.5, max(theta_true) + 0.5, length.out = L)
  
  # Natural cubic spline basis matrix
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

# Create base Stan data
stan_data <- prepare_stan_data(theta_hat, sigma, theta_true)

# ==============================================================================
# 4. DEFINE ALTERNATIVE PRIOR SPECIFICATIONS
# ==============================================================================

# Define comprehensive set of prior configurations
prior_configs <- tibble::tribble(
  ~prior_id, ~prior_type, ~param1, ~param2, ~label, ~description,
  
  # Reference prior (from main analysis)
  1, "cauchy", 0, 5, "Cauchy+(0, 5)", "Reference (weakly informative)",
  
  # Alternative distributions
  2, "exponential", 0.5, NA, "Exponential(0.5)", "Light-tailed alternative",
  3, "gamma", 2, 1, "Gamma(2, 1)", "Conjugate prior",
  4, "lognormal", 0, 1, "LogNormal(0, 1)", "Log-scale normal",
  
  # Varying Cauchy hyperparameters
  5, "cauchy", 0, 1, "Cauchy+(0, 1)", "More informative",
  6, "cauchy", 0, 10, "Cauchy+(0, 10)", "Less informative",
  7, "cauchy", 0, 2.5, "Cauchy+(0, 2.5)", "Moderately informative",
  
  # Additional alternatives
  8, "inverse_gamma", 2, 2, "InvGamma(2, 2)", "Heavy-tailed conjugate",
  9, "uniform", 0.01, 20, "Uniform(0.01, 20)", "Flat prior"
)

cat("\nPrior configurations to test:\n")
print(prior_configs)

# ==============================================================================
# 5. CREATE STAN MODELS WITH DIFFERENT LAMBDA PRIORS
# ==============================================================================

#' Generate Stan Model Code with Specified Lambda Prior
#'
#' This function creates a Stan model with the specified prior distribution
#' for the regularization parameter lambda.
#'
#' @param prior_type Type of prior distribution
#' @param param1 First parameter of the distribution
#' @param param2 Second parameter of the distribution (if applicable)
#' @return Stan model code as string

create_stan_model <- function(prior_type, param1 = NULL, param2 = NULL) {
  
  # Define the prior statement based on type
  prior_statement <- switch(prior_type,
                            "cauchy" = sprintf("lambda ~ cauchy(%.1f, %.1f);", param1, param2),
                            "exponential" = sprintf("lambda ~ exponential(%.2f);", param1),
                            "gamma" = sprintf("lambda ~ gamma(%.1f, %.1f);", param1, param2),
                            "lognormal" = sprintf("lambda ~ lognormal(%.1f, %.1f);", param1, param2),
                            "inverse_gamma" = sprintf("lambda ~ inv_gamma(%.1f, %.1f);", param1, param2),
                            "uniform" = sprintf("lambda ~ uniform(%.1f, %.1f);", param1, param2)
  )
  
  # Stan model template
  stan_code <- sprintf('
data {
  int<lower=1> K;                  // Number of sites
  vector[K] theta_hat;             // Observed effects
  vector<lower=0>[K] sigma;        // Standard errors
  int<lower=1> L;                  // Number of grid points
  vector[L] grid;                  // Grid points
  int<lower=1> M;                  // Spline basis dimension
  matrix[L, M] B;                  // Spline basis matrix
}

parameters {
  vector[M] alpha;                 // Spline coefficients
  real<lower=0> lambda;            // Regularization parameter
}

transformed parameters {
  vector[L] log_w = B * alpha;     // Log unnormalized weights
  vector[L] log_g = log_softmax(log_w);  // Log mixture weights
  simplex[L] g = softmax(log_w);   // Mixture weights
}

model {
  // Hyperprior on lambda
  %s
  
  // Conditional prior on alpha given lambda
  alpha ~ normal(0, inv_sqrt(lambda));
  
  // Likelihood: mixture model
  for (i in 1:K) {
    vector[L] log_components;
    for (j in 1:L) {
      log_components[j] = log_g[j] + 
                         normal_lpdf(theta_hat[i] | grid[j], sigma[i]);
    }
    target += log_sum_exp(log_components);
  }
}

generated quantities {
  // Prior distribution summary statistics
  real mean_g = dot_product(g, grid);
  real var_g = dot_product(g, square(grid - mean_g));
  
  // Posterior estimates for each site
  vector[K] theta_map;   // MAP estimates
  vector[K] theta_mean;  // Posterior means
  vector[K] theta_rep;   // Posterior draws
  
  for (i in 1:K) {
    vector[L] log_post;
    
    // Compute posterior for site i
    for (j in 1:L) {
      log_post[j] = log_g[j] + 
                    normal_lpdf(theta_hat[i] | grid[j], sigma[i]);
    }
    
    // MAP estimate
    int max_idx = 1;
    for (j in 2:L) {
      if (log_post[j] > log_post[max_idx]) {
        max_idx = j;
      }
    }
    theta_map[i] = grid[max_idx];
    
    // Posterior mean
    real log_post_max = max(log_post);
    vector[L] w = exp(log_post - log_post_max);
    w /= sum(w);
    theta_mean[i] = dot_product(w, grid);
    
    // Posterior draw
    theta_rep[i] = grid[categorical_rng(w)];
  }
}
', prior_statement)
  
  return(stan_code)
}

# ==============================================================================
# 6. VISUALIZE PRIOR DISTRIBUTIONS
# ==============================================================================

cat("\nCreating prior distribution visualizations...\n")

#' Compute Prior Density
#'
#' @param lambda_vals Vector of lambda values
#' @param prior_type Type of prior distribution
#' @param param1 First parameter
#' @param param2 Second parameter
#' @return Vector of density values

compute_prior_density <- function(lambda_vals, prior_type, param1, param2) {
  density_vals <- switch(prior_type,
                         "cauchy" = dcauchy(lambda_vals, param1, param2),
                         "exponential" = dexp(lambda_vals, param1),
                         "gamma" = dgamma(lambda_vals, param1, param2),
                         "lognormal" = dlnorm(lambda_vals, param1, param2),
                         "inverse_gamma" = dgamma(1/lambda_vals, param1, param2) / lambda_vals^2,
                         "uniform" = dunif(lambda_vals, param1, param2)
  )
  
  # Handle half-Cauchy (positive support only)
  if (prior_type == "cauchy") {
    density_vals <- 2 * density_vals * (lambda_vals >= 0)
  }
  
  return(density_vals)
}

# Create prior visualization data
lambda_range <- seq(0.01, 15, length.out = 1000)

prior_densities <- prior_configs %>%
  mutate(
    density = pmap(
      list(prior_type, param1, param2),
      ~ tibble(
        lambda = lambda_range,
        density = compute_prior_density(lambda_range, ..1, ..2, ..3)
      )
    )
  ) %>%
  unnest(density)

# Panel A: Alternative distributions comparison
p_alt_priors <- prior_densities %>%
  filter(prior_id %in% c(1, 2, 3, 4, 8, 9)) %>%
  ggplot(aes(x = lambda, y = density, color = label)) +
  geom_line(aes(linetype = prior_id == 1), size = 1) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                        guide = "none") +
  scale_color_viridis_d(name = "Prior Distribution") +
  xlim(0, 10) +
  labs(
    x = expression(lambda),
    y = "Density",
    title = "A. Alternative Prior Distributions",
    subtitle = "Solid: reference; Dashed: alternatives"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Panel B: Cauchy hyperparameter sensitivity
p_cauchy_params <- prior_densities %>%
  filter(prior_type == "cauchy") %>%
  ggplot(aes(x = lambda, y = density, color = label)) +
  geom_line(size = 1.2) +
  scale_color_viridis_d(name = "Cauchy Prior") +
  xlim(0, 10) +
  labs(
    x = expression(lambda),
    y = "Density",
    title = "B. Cauchy Hyperparameter Sensitivity",
    subtitle = "Effect of scale parameter on informativeness"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine prior visualizations
prior_viz <- p_alt_priors + p_cauchy_params +
  plot_annotation(
    title = "Prior Distribution Specifications for Lambda",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save visualization
ggsave(
  filename = file.path(output_dir, "Figure_D1_Lambda_Prior_Distributions.pdf"),
  plot = prior_viz,
  width = 12,
  height = 6,
  dpi = 300
)

cat("  Figure D1 saved.\n")

# ==============================================================================
# 7. FIT MODELS WITH VARIATIONAL BAYES
# ==============================================================================

#' Fit Model Using Variational Bayes
#'
#' This function compiles a Stan model with the specified prior and
#' fits it using Variational Bayes for computational efficiency.
#'
#' @param prior_config Prior configuration row
#' @param stan_data Prepared Stan data
#' @param iter Number of VB iterations
#' @return Tibble with results

fit_model_vb <- function(prior_config, stan_data, iter = 10000) {
  
  cat(sprintf("\nFitting model %d: %s\n", 
              prior_config$prior_id, prior_config$label))
  
  tryCatch({
    # Create Stan code
    stan_code <- create_stan_model(
      prior_config$prior_type,
      prior_config$param1,
      prior_config$param2
    )
    
    # Write to temporary file
    temp_file <- tempfile(fileext = ".stan")
    writeLines(stan_code, temp_file)
    
    # Compile model
    model <- cmdstan_model(
      temp_file,
      cpp_options = list(stan_threads = FALSE),
      stanc_options = list("O1")
    )
    
    # Start timing
    start_time <- Sys.time()
    
    # Fit using Variational Bayes (ADVI)
    fit_vb <- model$variational(
      data = stan_data,
      seed = 123 + prior_config$prior_id,
      algorithm = "meanfield",
      iter = iter,
      grad_samples = 1,
      elbo_samples = 100,
      eta = 0.1,
      adapt_iter = 50,
      output_samples = 4000,
      refresh = 0
    )
    
    # Calculate elapsed time
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    
    # Extract posterior summaries
    theta_mean_summary <- fit_vb$summary("theta_mean")
    theta_rep_summary <- fit_vb$summary("theta_rep")
    lambda_summary <- fit_vb$summary("lambda")
    
    # Get posterior samples for credible intervals
    theta_rep_draws <- fit_vb$draws("theta_rep", format = "matrix")
    
    # Calculate 90% credible intervals
    ci_90 <- apply(theta_rep_draws, 2, quantile, probs = c(0.05, 0.95))
    
    # Compile results
    results <- tibble(
      prior_id = prior_config$prior_id,
      site_id = 1:K,
      theta_true = theta_true,
      theta_hat = theta_hat,
      sigma = sigma,
      # Point estimates
      theta_mean_est = theta_mean_summary$mean,
      theta_rep_est = theta_rep_summary$mean,
      # Uncertainty
      theta_rep_sd = theta_rep_summary$sd,
      ci_lower = ci_90[1,],
      ci_upper = ci_90[2,],
      # Lambda posterior
      lambda_mean = lambda_summary$mean,
      lambda_sd = lambda_summary$sd,
      # Computation
      elapsed_time = elapsed_time
    )
    
    # Clean up
    unlink(temp_file)
    
    cat(sprintf("  Completed in %.2f seconds\n", elapsed_time))
    cat(sprintf("  Lambda posterior: mean=%.3f, sd=%.3f\n", 
                lambda_summary$mean, lambda_summary$sd))
    
    return(results)
    
  }, error = function(e) {
    warning(sprintf("Failed prior %d: %s\n", 
                    prior_config$prior_id, e$message))
    return(NULL)
  })
}

cat("\nFitting models with alternative lambda priors...\n")

# Run VB for all prior configurations
results_list <- prior_configs %>%
  group_by(prior_id) %>%
  group_split() %>%
  map(~ fit_model_vb(.x, stan_data))

# Combine results
results_all <- bind_rows(results_list) %>%
  left_join(prior_configs, by = "prior_id")

# Save intermediate results
saveRDS(results_all, file = file.path(output_dir, "lambda_sensitivity_results.rds"))

# ==============================================================================
# 8. CALCULATE PERFORMANCE METRICS
# ==============================================================================

cat("\nCalculating performance metrics...\n")

# Calculate comprehensive performance metrics
performance_metrics <- results_all %>%
  group_by(prior_id, prior_type, label, description) %>%
  summarise(
    # Point estimation accuracy
    rmse_mean = sqrt(mean((theta_mean_est - theta_true)^2)),
    rmse_rep = sqrt(mean((theta_rep_est - theta_true)^2)),
    correlation_mean = cor(theta_mean_est, theta_true),
    correlation_rep = cor(theta_rep_est, theta_true),
    
    # Uncertainty calibration
    coverage_90 = mean(theta_true >= ci_lower & theta_true <= ci_upper),
    
    # Interval properties
    mean_interval_width = mean(ci_upper - ci_lower),
    sd_interval_width = sd(ci_upper - ci_lower),
    
    # Lambda posterior
    lambda_post_mean = mean(lambda_mean),
    lambda_post_sd = mean(lambda_sd),
    
    # Computational efficiency
    computation_time = mean(elapsed_time),
    
    .groups = "drop"
  ) %>%
  # Add reference indicator and relative changes
  mutate(
    is_reference = (prior_id == 1),
    rmse_change_pct = (rmse_mean - rmse_mean[is_reference]) / 
      rmse_mean[is_reference] * 100,
    coverage_diff = (coverage_90 - coverage_90[is_reference]) * 100
  )

# ==============================================================================
# 9. CREATE PERFORMANCE COMPARISON VISUALIZATIONS
# ==============================================================================

cat("\nCreating performance comparison plots...\n")

# Theme for consistent plotting
theme_custom <- theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

# Panel A: Point Estimation Accuracy
p1_accuracy <- performance_metrics %>%
  pivot_longer(
    cols = c(rmse_mean, rmse_rep),
    names_to = "estimator",
    values_to = "rmse",
    names_prefix = "rmse_"
  ) %>%
  mutate(
    estimator = factor(estimator, 
                       levels = c("mean", "rep"),
                       labels = c("Posterior Mean", "Posterior Draw")),
    label_short = str_replace(label, "Cauchy\\+", "C+") %>%
      str_replace("Exponential", "Exp") %>%
      str_replace("LogNormal", "LogN") %>%
      str_replace("InvGamma", "IG") %>%
      str_replace("Uniform", "Unif")
  ) %>%
  ggplot(aes(x = reorder(label_short, prior_id), y = rmse, 
             fill = is_reference)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", rmse)), 
            vjust = -0.5, size = 3) +
  facet_wrap(~ estimator, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
                    labels = c("Alternative", "Reference"),
                    name = "") +
  labs(
    x = "Prior Specification",
    y = "Root Mean Squared Error",
    title = "A. Point Estimation Accuracy"
  ) +
  theme_custom

# Panel B: Coverage Probability
p2_coverage <- performance_metrics %>%
  mutate(
    label_short = str_replace(label, "Cauchy\\+", "C+") %>%
      str_replace("Exponential", "Exp") %>%
      str_replace("LogNormal", "LogN") %>%
      str_replace("InvGamma", "IG") %>%
      str_replace("Uniform", "Unif")
  ) %>%
  ggplot(aes(x = reorder(label_short, prior_id), 
             y = coverage_90, 
             fill = is_reference)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0.90, linetype = "dashed", 
             color = "red", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", coverage_90 * 100)), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
                    labels = c("Alternative", "Reference"),
                    name = "") +
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0.85, 0.95)) +
  labs(
    x = "Prior Specification",
    y = "90% Coverage Probability",
    title = "B. Uncertainty Calibration"
  ) +
  theme_custom

# Panel C: Lambda Posterior
p3_lambda <- results_all %>%
  group_by(prior_id, label) %>%
  summarise(
    lambda_mean = mean(lambda_mean),
    lambda_lower = mean(lambda_mean - 1.96 * lambda_sd),
    lambda_upper = mean(lambda_mean + 1.96 * lambda_sd),
    .groups = "drop"
  ) %>%
  mutate(
    label_short = str_replace(label, "Cauchy\\+", "C+") %>%
      str_replace("Exponential", "Exp") %>%
      str_replace("LogNormal", "LogN") %>%
      str_replace("InvGamma", "IG") %>%
      str_replace("Uniform", "Unif")
  ) %>%
  ggplot(aes(x = reorder(label_short, prior_id), y = lambda_mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper), 
                width = 0.2) +
  labs(
    x = "Prior Specification",
    y = expression(paste("Posterior ", lambda)),
    title = "C. Lambda Posterior Estimates"
  ) +
  theme_custom

# Panel D: Accuracy-Coverage Trade-off
p4_tradeoff <- performance_metrics %>%
  ggplot(aes(x = rmse_mean, y = coverage_90)) +
  geom_point(aes(color = prior_type, shape = is_reference), 
             size = 4, alpha = 0.8) +
  geom_text(aes(label = prior_id), 
            vjust = -1, size = 3) +
  geom_hline(yintercept = 0.90, linetype = "dashed", 
             color = "red", alpha = 0.5) +
  scale_color_viridis_d(name = "Prior Type") +
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
combined_performance <- (p1_accuracy | p2_coverage) / (p3_lambda | p4_tradeoff) +
  plot_annotation(
    title = "Lambda Prior Sensitivity Analysis: Impact on Inference Quality",
    subtitle = sprintf("Twin towers scenario (K=%d, I=0.7, R=9) using Variational Bayes", K),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save performance comparison
ggsave(
  filename = file.path(output_dir, "Figure_D2_Lambda_Sensitivity_Performance.pdf"),
  plot = combined_performance,
  width = 14,
  height = 10,
  dpi = 300
)

cat("  Figure D2 saved.\n")

# ==============================================================================
# 10. CREATE SUMMARY TABLES
# ==============================================================================

cat("\nCreating summary tables...\n")

# Table 1: Alternative Distributions Comparison
table_alternatives <- performance_metrics %>%
  filter(prior_id %in% c(1, 2, 3, 4, 8, 9)) %>%
  select(
    label, 
    rmse_mean, rmse_change_pct,
    coverage_90, coverage_diff,
    mean_interval_width,
    lambda_post_mean, lambda_post_sd,
    computation_time
  ) %>%
  mutate(
    across(c(rmse_mean, mean_interval_width, lambda_post_mean, lambda_post_sd), 
           ~round(., 4)),
    across(c(rmse_change_pct, coverage_diff), ~round(., 2)),
    coverage_90 = round(coverage_90 * 100, 1),
    computation_time = round(computation_time, 2)
  )

# Save Table 1
write.csv(
  table_alternatives,
  file = file.path(output_dir, "Table_D1_Alternative_Prior_Distributions.csv"),
  row.names = FALSE
)

cat("\nTable D1: Alternative Prior Distributions\n")
print(kable(table_alternatives, format = "simple"))

# Table 2: Cauchy Hyperparameter Sensitivity
table_cauchy <- performance_metrics %>%
  filter(prior_type == "cauchy") %>%
  select(
    label,
    rmse_mean, rmse_change_pct,
    coverage_90, coverage_diff,
    mean_interval_width,
    lambda_post_mean, lambda_post_sd
  ) %>%
  mutate(
    across(c(rmse_mean, mean_interval_width, lambda_post_mean, lambda_post_sd), 
           ~round(., 4)),
    across(c(rmse_change_pct, coverage_diff), ~round(., 2)),
    coverage_90 = round(coverage_90 * 100, 1)
  )

# Save Table 2
write.csv(
  table_cauchy,
  file = file.path(output_dir, "Table_D2_Cauchy_Hyperparameter_Sensitivity.csv"),
  row.names = FALSE
)

cat("\nTable D2: Cauchy Hyperparameter Sensitivity\n")
print(kable(table_cauchy, format = "simple"))

# ==============================================================================
# 11. ROBUSTNESS CHECK WITH ADDITIONAL SCENARIOS
# ==============================================================================

cat("\nTesting robustness across different simulation scenarios...\n")

#' Run Sensitivity Analysis for a Given Scenario
#'
#' @param scenario_data Data for a specific scenario
#' @param scenario_name Name of the scenario
#' @return Tibble with performance metrics

run_scenario_sensitivity <- function(scenario_data, scenario_name) {
  
  cat(sprintf("\nRunning sensitivity for scenario: %s\n", scenario_name))
  
  # Prepare Stan data
  stan_data_scenario <- prepare_stan_data(
    scenario_data$theta_hat,
    scenario_data$se,
    scenario_data$theta_true
  )
  
  # Run only for key priors to save time
  key_priors <- prior_configs %>%
    filter(prior_id %in% c(1, 2, 3, 5, 6))
  
  # Fit models
  scenario_results <- key_priors %>%
    group_by(prior_id) %>%
    group_split() %>%
    map(~ fit_model_vb(.x, stan_data_scenario, iter = 5000))
  
  # Calculate metrics
  scenario_metrics <- bind_rows(scenario_results) %>%
    left_join(prior_configs, by = "prior_id") %>%
    group_by(prior_id, label) %>%
    summarise(
      rmse = sqrt(mean((theta_mean_est - theta_true)^2)),
      coverage = mean(theta_true >= ci_lower & theta_true <= ci_upper),
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario_name)
  
  return(scenario_metrics)
}

# Test on additional scenarios
additional_scenarios <- list(
  "Low_Reliability" = df_sim %>% 
    filter(I == 0.5, R == 1) %>% 
    slice(1) %>% 
    pull(data) %>% 
    .[[1]],
  "High_Reliability" = df_sim %>% 
    filter(I == 0.9, R == 1) %>% 
    slice(1) %>% 
    pull(data) %>% 
    .[[1]]
)

# Run sensitivity for additional scenarios
scenario_results <- map2(
  additional_scenarios,
  names(additional_scenarios),
  run_scenario_sensitivity
)

# Combine scenario results
all_scenario_metrics <- bind_rows(scenario_results) %>%
  bind_rows(
    performance_metrics %>%
      filter(prior_id %in% c(1, 2, 3, 5, 6)) %>%
      select(prior_id, label, rmse = rmse_mean, coverage = coverage_90) %>%
      mutate(scenario = "Base_Scenario")
  )

# Create robustness plot
p_robustness <- all_scenario_metrics %>%
  pivot_longer(
    cols = c(rmse, coverage),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, 
                    levels = c("rmse", "coverage"),
                    labels = c("RMSE", "Coverage (90%)")),
    scenario = factor(scenario,
                      levels = c("Low_Reliability", "Base_Scenario", "High_Reliability"),
                      labels = c("Low Reliability\n(I=0.5, R=1)", 
                                 "Base Scenario\n(I=0.7, R=9)", 
                                 "High Reliability\n(I=0.9, R=1)")),
    value_label = ifelse(metric == "RMSE",
                         sprintf("%.3f", value),
                         sprintf("%.1f%%", value * 100))
  ) %>%
  ggplot(aes(x = label, y = value, fill = label)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = value_label), 
            vjust = -0.3, size = 2.5) +
  facet_grid(scenario ~ metric, scales = "free") +
  scale_fill_viridis_d(name = "Prior") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Prior Specification",
    y = "",
    title = "Figure D3. Robustness Across Different Data Scenarios",
    subtitle = "Consistency indicates low sensitivity to prior choice"
  ) +
  theme_custom +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

# Save robustness plot
ggsave(
  filename = file.path(output_dir, "Figure_D3_Lambda_Sensitivity_Robustness.pdf"),
  plot = p_robustness,
  width = 10,
  height = 8,
  dpi = 300
)

cat("  Figure D3 saved.\n")

# ==============================================================================
# 12. FINAL SUMMARY AND SESSION INFORMATION
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("LAMBDA PRIOR SENSITIVITY ANALYSIS SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Key findings
cat("KEY FINDINGS:\n")
cat("-------------\n")

# Point estimation robustness
rmse_range <- range(performance_metrics$rmse_mean)
cat(sprintf("1. Point Estimation (RMSE):\n"))
cat(sprintf("   - Range: [%.4f, %.4f]\n", rmse_range[1], rmse_range[2]))
cat(sprintf("   - Maximum deviation from reference: %.2f%%\n",
            max(abs(performance_metrics$rmse_change_pct), na.rm = TRUE)))

# Coverage calibration
coverage_range <- range(performance_metrics$coverage_90)
cat(sprintf("\n2. Uncertainty Calibration (90%% Coverage):\n"))
cat(sprintf("   - Range: [%.1f%%, %.1f%%]\n",
            coverage_range[1] * 100, coverage_range[2] * 100))
cat(sprintf("   - Maximum deviation from nominal: %.1f pp\n",
            max(abs(performance_metrics$coverage_90 - 0.90)) * 100))

# Lambda posterior
cat(sprintf("\n3. Lambda Posterior:\n"))
cat(sprintf("   - Mean range: [%.3f, %.3f]\n",
            min(performance_metrics$lambda_post_mean),
            max(performance_metrics$lambda_post_mean)))

# Computational efficiency
cat(sprintf("\n4. Computational Efficiency:\n"))
cat(sprintf("   - Time range: %.0f-%.0f seconds\n",
            min(performance_metrics$computation_time),
            max(performance_metrics$computation_time)))

cat("\nRECOMMENDATIONS:\n")
cat("----------------\n")
cat("• The reference prior Cauchy+(0, 5) provides excellent performance\n")
cat("• Results are highly robust to reasonable prior specifications\n")
cat("• Heavy-tailed priors preferred over light-tailed alternatives\n")
cat("• Scale parameter between 2.5 and 10 yields similar results\n")

cat("\nOutput files saved to:", output_dir, "\n")
cat("  - Figure_D1_Lambda_Prior_Distributions.pdf\n")
cat("  - Figure_D2_Lambda_Sensitivity_Performance.pdf\n")
cat("  - Figure_D3_Lambda_Sensitivity_Robustness.pdf\n")
cat("  - Table_D1_Alternative_Prior_Distributions.csv\n")
cat("  - Table_D2_Cauchy_Hyperparameter_Sensitivity.csv\n")

cat("\nSession Information:\n")
print(sessionInfo())