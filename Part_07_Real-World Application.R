################################################################################
#
# Replication Code for:
# "Fully Bayesian Inference for Meta-Analytic Deconvolution 
#  Using Efron's Log-Spline Prior"
#
# Section 5: Real-World Application - Firm-Level Labor Market Discrimination
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: January 2025
#
# Description:
# This script applies the fully Bayesian framework to analyze firm-level
# discrimination in hiring using audit study data from Kline, Rose, and 
# Walters (2024). We estimate both racial (White vs Black) and gender 
# (Male vs Female) discrimination parameters with proper uncertainty quantification.
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
  "tidyverse", "cmdstanr", "splines", "patchwork", "scales", 
  "viridis", "kableExtra", "ggrepel", "posterior", "bayesplot",
  "fixest", "haven"
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
output_dir <- "./results"  # Output directory for results

# Create directories if they don't exist
for (dir in c(data_dir, stan_dir, output_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

setwd(work_dir)

# ==============================================================================
# 2. DATA PREPARATION
# ==============================================================================

cat("Loading and preparing discrimination audit study data...\n")

#' Load and Prepare Discrimination Data
#'
#' This function loads the Kline, Rose, and Walters (2024) data and
#' calculates firm-level discrimination estimates using OLS.
#'
#' @param data_path Path to the data file
#' @return Data frame with firm-level discrimination estimates

prepare_discrimination_data <- function(data_path) {
  
  # Load data (assuming Stata format - adjust as needed)
  df <- haven::read_dta(data_path)
  
  # Step 1: Calculate firm-level callback rates and filter firms
  # Following KRW criteria: callback rate >= 3% and at least 40 jobs
  firm_stats <- df %>%
    group_by(firm_id) %>%
    summarise(
      p = mean(callback),
      n_jobs = n_distinct(job_id),
      firm_name = first(firm_name),
      .groups = "drop"
    ) %>%
    filter(p >= 0.03 & n_jobs >= 40)
  
  # Filter main dataset to qualifying firms
  df_filtered <- df %>%
    filter(firm_id %in% firm_stats$firm_id)
  
  # Create sequential firm identifier
  firm_stats <- firm_stats %>%
    arrange(firm_id) %>%
    mutate(j = row_number())
  
  # Add sequential identifier to main data
  df_filtered <- df_filtered %>%
    left_join(firm_stats %>% select(firm_id, j), by = "firm_id")
  
  cat(sprintf("Number of firms after filtering: %d\n", max(df_filtered$j)))
  cat(sprintf("Number of observations: %d\n", nrow(df_filtered)))
  
  return(df_filtered)
}

#' Calculate Firm-Level Discrimination Gaps
#'
#' Estimates firm-specific discrimination parameters using OLS with
#' job-level clustering for standard errors.
#'
#' @param data Filtered dataset
#' @return Data frame with firm-level gap estimates

calculate_discrimination_gaps <- function(data) {
  
  results_list <- list()
  
  for (firm in unique(data$j)) {
    firm_data <- data %>% filter(j == firm)
    
    # White vs Black gap regression
    white_reg <- fixest::feols(
      callback ~ white, 
      data = firm_data, 
      cluster = ~job_id
    )
    
    # Male vs Female gap regression  
    male_reg <- fixest::feols(
      callback ~ male,
      data = firm_data,
      cluster = ~job_id
    )
    
    results_list[[firm]] <- data.frame(
      j = firm,
      firm_id = unique(firm_data$firm_id),
      firm_name = unique(firm_data$firm_name),
      theta_white = coef(white_reg)["white"],
      s_white = se(white_reg)["white"],
      theta_male = coef(male_reg)["male"],
      s_male = se(male_reg)["male"],
      n_jobs = n_distinct(firm_data$job_id),
      n_obs = nrow(firm_data)
    )
  }
  
  results_df <- bind_rows(results_list)
  return(results_df)
}

# Load and prepare data
# Note: Users need to provide the path to KRW data
data_path <- file.path(data_dir, "krw_data.dta")

if (file.exists(data_path)) {
  df_filtered <- prepare_discrimination_data(data_path)
  firm_gaps <- calculate_discrimination_gaps(df_filtered)
} else {
  # Create simulated example data for demonstration if real data not available
  cat("Real data not found. Creating simulated example data...\n")
  
  set.seed(2025)
  n_firms <- 100
  
  # Simulate firm-level discrimination estimates
  firm_gaps <- tibble(
    j = 1:n_firms,
    firm_id = 1000 + 1:n_firms,
    firm_name = paste("Firm", 1:n_firms),
    # Racial discrimination: bimodal distribution
    theta_white = c(rnorm(n_firms/2, -0.02, 0.03), 
                    rnorm(n_firms/2, 0.05, 0.04)),
    s_white = runif(n_firms, 0.01, 0.05),
    # Gender discrimination: unimodal distribution
    theta_male = rnorm(n_firms, 0.02, 0.03),
    s_male = runif(n_firms, 0.01, 0.04),
    n_jobs = sample(40:200, n_firms, replace = TRUE),
    n_obs = sample(400:2000, n_firms, replace = TRUE)
  )
}

# Save prepared data
saveRDS(firm_gaps, file = file.path(data_dir, "firm_discrimination_estimates.rds"))

# Display summary statistics
cat("\nSummary of discrimination estimates:\n")
summary_stats <- firm_gaps %>%
  summarise(
    across(starts_with("theta"), 
           list(mean = ~mean(., na.rm = TRUE),
                sd = ~sd(., na.rm = TRUE),
                min = ~min(., na.rm = TRUE),
                max = ~max(., na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  )
print(summary_stats)

# ==============================================================================
# 3. PREPARE STAN DATA
# ==============================================================================

#' Prepare Stan Data for Discrimination Analysis
#'
#' @param theta_hat Vector of firm-level discrimination estimates
#' @param sigma Vector of standard errors
#' @param L Number of grid points (default: 101)
#' @param M Degrees of freedom for splines (default: 6)
#' @return List containing Stan data

prepare_discrimination_stan_data <- function(theta_hat, sigma, L = 101, M = 6) {
  
  # Remove any NA values
  valid_idx <- !is.na(theta_hat) & !is.na(sigma) & sigma > 0
  theta_hat <- theta_hat[valid_idx]
  sigma <- sigma[valid_idx]
  
  K <- length(theta_hat)
  
  # Create grid spanning the range with expansion
  range_theta <- range(theta_hat, na.rm = TRUE)
  expansion <- 0.5 * diff(range_theta)
  
  grid <- seq(
    from = range_theta[1] - expansion,
    to = range_theta[2] + expansion,
    length.out = L
  )
  
  # Create spline basis matrix
  B <- splines::ns(grid, df = M, intercept = FALSE)
  
  # Return Stan data list
  list(
    K = K,
    theta_hat = theta_hat,
    sigma = sigma,
    L = L,
    grid = grid,
    M = M,
    B = B
  )
}

# Prepare data for racial discrimination
stan_data_race <- prepare_discrimination_stan_data(
  theta_hat = firm_gaps$theta_white,
  sigma = firm_gaps$s_white
)

# Prepare data for gender discrimination  
stan_data_gender <- prepare_discrimination_stan_data(
  theta_hat = firm_gaps$theta_male,
  sigma = firm_gaps$s_male
)

cat(sprintf("\nRacial discrimination analysis: K = %d firms\n", stan_data_race$K))
cat(sprintf("Gender discrimination analysis: K = %d firms\n", stan_data_gender$K))

# ==============================================================================
# 4. COMPILE AND FIT STAN MODELS
# ==============================================================================

cat("\nCompiling Stan model...\n")

# Check if Stan model exists
stan_file <- file.path(stan_dir, "efron_log_spline_model.stan")

if (!file.exists(stan_file)) {
  stop("Stan model file not found. Please ensure 'efron_log_spline_model.stan' ",
       "is in the stan directory.")
}

# Compile with threading for MCMC
compiled_model_mcmc <- cmdstan_model(
  stan_file,
  cpp_options = list(stan_threads = TRUE),
  stanc_options = list("O1")
)

# Compile without threading for optimization methods
compiled_model_opt <- cmdstan_model(
  stan_file,
  cpp_options = list(stan_threads = FALSE),
  stanc_options = list("O1")
)

cat("Models compiled successfully.\n")

#' Fit All Estimation Methods
#'
#' Fits the model using multiple methods: MLE, MAP, VB, Laplace, 
#' Pathfinder, and full MCMC.
#'
#' @param stan_data Prepared Stan data
#' @param compiled_model_opt Model for optimization methods
#' @param compiled_model_mcmc Model for MCMC
#' @param seed Random seed
#' @param label Analysis label
#' @return List of fitted models

fit_all_methods <- function(stan_data, compiled_model_opt, compiled_model_mcmc, 
                            seed = 123, label = "analysis") {
  
  cat(sprintf("\n===== Fitting %s =====\n", label))
  results <- list()
  
  # 1. Empirical Bayes (MLE)
  cat("  1. Fitting Empirical Bayes (MLE)...")
  start_time <- Sys.time()
  tryCatch({
    results$mle <- compiled_model_opt$optimize(
      data = stan_data,
      seed = seed,
      algorithm = "lbfgs",
      init = 0,
      jacobian = FALSE
    )
    cat(sprintf(" Done (%.1f sec)\n", 
                difftime(Sys.time(), start_time, units = "secs")))
  }, error = function(e) {
    cat(" Failed\n")
    results$mle <- NULL
  })
  
  # 2. MAP
  cat("  2. Fitting MAP...")
  start_time <- Sys.time()
  tryCatch({
    results$map <- compiled_model_opt$optimize(
      data = stan_data,
      seed = seed,
      algorithm = "lbfgs",
      init = 0,
      jacobian = TRUE
    )
    cat(sprintf(" Done (%.1f sec)\n", 
                difftime(Sys.time(), start_time, units = "secs")))
  }, error = function(e) {
    cat(" Failed\n")
    results$map <- NULL
  })
  
  # 3. Variational Bayes
  cat("  3. Fitting Variational Bayes...")
  start_time <- Sys.time()
  tryCatch({
    results$vb <- compiled_model_opt$variational(
      data = stan_data,
      seed = seed,
      algorithm = "meanfield",
      iter = 10000,
      grad_samples = 1,
      elbo_samples = 100,
      eta = 0.1,
      adapt_iter = 50,
      output_samples = 4000,
      refresh = 0
    )
    cat(sprintf(" Done (%.1f sec)\n", 
                difftime(Sys.time(), start_time, units = "secs")))
  }, error = function(e) {
    cat(" Failed\n")
    results$vb <- NULL
  })
  
  # 4. Full MCMC (Fully Bayesian)
  cat("  4. Fitting MCMC (this may take a few minutes)...")
  start_time <- Sys.time()
  tryCatch({
    results$mcmc <- compiled_model_mcmc$sample(
      data = stan_data,
      seed = seed,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 3000,
      refresh = 0,
      adapt_delta = 0.9,
      threads_per_chain = 2
    )
    cat(sprintf(" Done (%.1f sec)\n", 
                difftime(Sys.time(), start_time, units = "secs")))
  }, error = function(e) {
    cat(" Failed\n")
    results$mcmc <- NULL
  })
  
  return(results)
}

# Fit models for racial discrimination
fits_race <- fit_all_methods(
  stan_data_race, 
  compiled_model_opt, 
  compiled_model_mcmc,
  label = "Racial Discrimination"
)

# Fit models for gender discrimination
fits_gender <- fit_all_methods(
  stan_data_gender,
  compiled_model_opt,
  compiled_model_mcmc,
  label = "Gender Discrimination"
)

# Save fitted objects
saveRDS(fits_race, file = file.path(output_dir, "fits_race_discrimination.rds"))
saveRDS(fits_gender, file = file.path(output_dir, "fits_gender_discrimination.rds"))

# ==============================================================================
# 5. EXTRACT POSTERIOR DISTRIBUTIONS
# ==============================================================================

cat("\nExtracting posterior distributions...\n")

#' Extract Posterior Distributions
#'
#' Extracts posterior distributions for both theta_rep (FB) and 
#' theta_mean (EB) parameters.
#'
#' @param mcmc_fit MCMC fit object
#' @param firm_data Original firm data
#' @return List with summary statistics and posterior draws

extract_posterior_distributions <- function(mcmc_fit, firm_data) {
  
  # Extract theta_rep (FB approach with parameter uncertainty)
  theta_rep_draws <- mcmc_fit$draws("theta_rep", format = "matrix")
  
  # Extract theta_mean (EB approach with estimator uncertainty)  
  theta_mean_draws <- mcmc_fit$draws("theta_mean", format = "matrix")
  
  # Calculate summaries for theta_rep
  theta_rep_summary <- tibble(
    firm_id = 1:ncol(theta_rep_draws),
    theta_rep_mean = colMeans(theta_rep_draws),
    theta_rep_q5 = apply(theta_rep_draws, 2, quantile, 0.05),
    theta_rep_q95 = apply(theta_rep_draws, 2, quantile, 0.95),
    theta_rep_sd = apply(theta_rep_draws, 2, sd),
    prob_positive_rep = colMeans(theta_rep_draws > 0)
  )
  
  # Calculate summaries for theta_mean
  theta_mean_summary <- tibble(
    firm_id = 1:ncol(theta_mean_draws),
    theta_mean_mean = colMeans(theta_mean_draws),
    theta_mean_q5 = apply(theta_mean_draws, 2, quantile, 0.05),
    theta_mean_q95 = apply(theta_mean_draws, 2, quantile, 0.95),
    theta_mean_sd = apply(theta_mean_draws, 2, sd),
    prob_positive_mean = colMeans(theta_mean_draws > 0)
  )
  
  # Combine with firm data
  combined_summary <- theta_rep_summary %>%
    left_join(theta_mean_summary, by = "firm_id") %>%
    bind_cols(firm_data) %>%
    mutate(
      significant_rep = (theta_rep_q5 > 0) | (theta_rep_q95 < 0),
      significant_mean = (theta_mean_q5 > 0) | (theta_mean_q95 < 0)
    )
  
  return(list(
    summary = combined_summary,
    theta_rep_draws = theta_rep_draws,
    theta_mean_draws = theta_mean_draws
  ))
}

# Extract distributions
race_distributions <- extract_posterior_distributions(fits_race$mcmc, firm_gaps)
gender_distributions <- extract_posterior_distributions(fits_gender$mcmc, firm_gaps)

# Save extracted distributions
saveRDS(race_distributions, 
        file = file.path(output_dir, "race_posterior_distributions.rds"))
saveRDS(gender_distributions, 
        file = file.path(output_dir, "gender_posterior_distributions.rds"))

# ==============================================================================
# 6. CREATE VISUALIZATIONS
# ==============================================================================

cat("\nCreating visualizations...\n")

# Theme for consistent plotting
theme_custom <- theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

#' Create Caterpillar Plot Comparing FB vs EB
#'
#' @param distributions Posterior distributions
#' @param discrimination_type Type of discrimination
#' @return ggplot object

create_caterpillar_comparison <- function(distributions, 
                                          discrimination_type = "Racial") {
  
  # Prepare data
  plot_data <- distributions$summary %>%
    arrange(theta_rep_mean) %>%
    mutate(rank = row_number())
  
  # Panel A: Fully Bayesian (theta_rep)
  panel_a <- plot_data %>%
    ggplot(aes(x = rank, y = theta_rep_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = theta_rep_q5, ymax = theta_rep_q95,
                      color = significant_rep), 
                  width = 0, alpha = 0.7) +
    geom_point(aes(color = significant_rep), size = 0.8) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray60"),
                       guide = "none") +
    labs(x = NULL, 
         y = "Discrimination Parameter",
         title = "Fully Bayesian Approach",
         subtitle = "Parameter uncertainty (theta_rep)") +
    theme_custom
  
  # Panel B: Empirical Bayes (theta_mean)
  panel_b <- plot_data %>%
    ggplot(aes(x = rank, y = theta_mean_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = theta_mean_q5, ymax = theta_mean_q95,
                      color = significant_mean), 
                  width = 0, alpha = 0.7) +
    geom_point(aes(color = significant_mean), size = 0.8) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray60"),
                       guide = "none") +
    labs(x = "Firms (ranked by posterior mean)",
         y = "Discrimination Parameter",
         title = "Empirical Bayes Approach",
         subtitle = "Estimator uncertainty (theta_mean)") +
    theme_custom
  
  # Combine panels
  combined_plot <- panel_a / panel_b +
    plot_annotation(
      title = paste(discrimination_type, "Discrimination: FB vs EB Comparison"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  return(combined_plot)
}

# Create caterpillar plots
caterpillar_race <- create_caterpillar_comparison(
  race_distributions, "Racial"
)
caterpillar_gender <- create_caterpillar_comparison(
  gender_distributions, "Gender"
)

# Save plots
ggsave(
  filename = file.path(output_dir, "Figure_10_Caterpillar_Race.pdf"),
  plot = caterpillar_race,
  width = 8,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "Figure_11_Caterpillar_Gender.pdf"),
  plot = caterpillar_gender,
  width = 8,
  height = 10,
  dpi = 300
)

#' Create Density Comparison Plot
#'
#' @param dist_data Posterior distributions
#' @param discrimination_type Type of discrimination
#' @param theta_col Column name for raw estimates
#' @return ggplot object

create_density_comparison <- function(dist_data, discrimination_type = "Racial", 
                                      theta_col = "theta_white") {
  
  # Prepare data
  plot_data <- tibble(
    value = c(dist_data$summary[[theta_col]],
              dist_data$summary$theta_rep_mean),
    type = rep(c("Raw OLS Estimates", "FB Posterior Means"), 
               each = nrow(dist_data$summary))
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = value, fill = type)) +
    geom_vline(xintercept = 0, color = "red", linetype = "solid", 
               size = 1, alpha = 0.7) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Raw OLS Estimates" = "#E69F00",
                                 "FB Posterior Means" = "#56B4E9"),
                      name = "Estimate Type") +
    labs(x = "Discrimination Parameter",
         y = "Density",
         title = paste(discrimination_type, "Discrimination"),
         subtitle = "Comparison of raw estimates vs posterior means") +
    theme_custom
  
  return(p)
}

# Create density plots
density_race <- create_density_comparison(
  race_distributions, "Racial", "theta_white"
)
density_gender <- create_density_comparison(
  gender_distributions, "Gender", "theta_male"
)

# Combine density plots
density_combined <- density_race / density_gender

# Save density plot
ggsave(
  filename = file.path(output_dir, "Figure_12_Density_Comparison.pdf"),
  plot = density_combined,
  width = 8,
  height = 10,
  dpi = 300
)

cat("  Figures 10-12 saved.\n")

# ==============================================================================
# 7. SUMMARY STATISTICS AND TABLES
# ==============================================================================

cat("\nCalculating summary statistics...\n")

# Create summary table
summary_table <- tibble(
  Type = c("Racial", "Gender"),
  `Mean (Raw)` = c(mean(firm_gaps$theta_white), mean(firm_gaps$theta_male)),
  `SD (Raw)` = c(sd(firm_gaps$theta_white), sd(firm_gaps$theta_male)),
  `Mean (FB)` = c(mean(race_distributions$summary$theta_rep_mean),
                  mean(gender_distributions$summary$theta_rep_mean)),
  `SD (FB)` = c(sd(race_distributions$summary$theta_rep_mean),
                sd(gender_distributions$summary$theta_rep_mean)),
  `% Significant` = c(
    mean(race_distributions$summary$significant_rep) * 100,
    mean(gender_distributions$summary$significant_rep) * 100
  )
) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

# Print and save table
cat("\nTable 4: Summary of Discrimination Estimates\n")
print(summary_table)

write.csv(
  summary_table,
  file = file.path(output_dir, "Table_4_Summary_Statistics.csv"),
  row.names = FALSE
)

# ==============================================================================
# 8. SESSION INFORMATION
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Real-world data analysis completed successfully!\n")
cat("\nKey findings:\n")
cat(sprintf("  - Racial discrimination: %.1f%% of firms show significant effects\n",
            mean(race_distributions$summary$significant_rep) * 100))
cat(sprintf("  - Gender discrimination: %.1f%% of firms show significant effects\n",
            mean(gender_distributions$summary$significant_rep) * 100))
cat(sprintf("  - Shrinkage from raw to FB: %.1f%% (racial), %.1f%% (gender)\n",
            (1 - sd(race_distributions$summary$theta_rep_mean) / 
               sd(firm_gaps$theta_white)) * 100,
            (1 - sd(gender_distributions$summary$theta_rep_mean) / 
               sd(firm_gaps$theta_male)) * 100))

cat("\nOutput files saved to:", output_dir, "\n")
cat("  - Figures 10-12: Caterpillar and density plots\n")
cat("  - Table 4: Summary statistics\n")
cat("  - Posterior distributions and fitted models\n")

cat("\nSession Information:\n")
print(sessionInfo())