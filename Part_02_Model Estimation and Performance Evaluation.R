################################################################################
#
# Replication Code for:
# "Fully Bayesian Inference for Meta-Analytic Deconvolution 
#  Using Efron's Log-Spline Prior"
#
# Part 2: Model Estimation and Performance Evaluation
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: January 2025
#
# Description:
# This script fits the Bayesian log-spline model using various estimation
# methods (MCMC, Empirical Bayes, Variational Bayes, etc.) and evaluates
# their performance across different simulation scenarios.
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
  "tidyverse", "splines", "cmdstanr", "posterior", 
  "bayesplot", "tidybayes", "furrr", "future", 
  "rstan", "tictoc", "progressr"
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

# Install CmdStan if not already installed
if (!cmdstanr::cmdstan_path()) {
  cmdstanr::install_cmdstan()
}

# Set working directory
# Users should modify these paths according to their system
work_dir <- "."  # Set your working directory here
data_dir <- "./data"  # Set your data directory here
stan_dir <- "./stan"  # Set your Stan files directory here

# Create directories if they don't exist
for (dir in c(data_dir, stan_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

setwd(work_dir)

# ==============================================================================
# 2. LOAD SIMULATED DATA
# ==============================================================================

cat("Loading simulated data...\n")

# Load the pre-simulated dataset from Part 1
df <- readRDS(file.path(data_dir, "simulated_data_varied_by_I_and_R.rds"))

cat("Data loaded successfully:\n")
print(df)

# ==============================================================================
# 3. PREPARE STAN DATA
# ==============================================================================

#' Prepare Stan Data from Simulated Dataset
#'
#' This function converts the simulated data into the format required
#' by the Stan model, including the spline basis matrix.
#'
#' @param df_sub Tibble containing simulated data for one scenario
#' @return List containing all required Stan data elements

prepare_stan_data <- function(df_sub) {
  # Extract data vectors
  K <- nrow(df_sub)
  theta_true <- df_sub %>% pull(theta_true)
  se <- df_sub %>% pull(se)
  theta_hat <- df_sub %>% pull(theta_hat)
  
  # Create grid spanning the support of true effects
  # Grid extends slightly beyond observed range for stability
  grid <- seq(min(theta_true) - 0.5, max(theta_true) + 0.5, length.out = 101)
  L <- length(grid)
  
  # Create natural cubic spline basis matrix
  # M = 6 provides sufficient flexibility for bimodal distributions
  M <- 6  
  B <- splines::ns(grid, df = M, intercept = FALSE)
  
  # Return Stan data list
  list(
    K = K,                # Number of sites
    theta_hat = theta_hat,  # Observed site effects
    sigma = se,            # Site-specific standard errors
    L = L,                 # Number of grid points
    grid = grid,           # Grid points for discretization
    M = M,                 # Spline basis dimension
    B = B                  # Spline basis matrix
  )
}

# Apply function to create stan_data column
cat("\nPreparing Stan data for each scenario...\n")

df2 <- df %>%
  mutate(stan_data = map(data, prepare_stan_data))

cat("Stan data prepared for", nrow(df2), "scenarios.\n")

# ==============================================================================
# 4. COMPILE STAN MODEL
# ==============================================================================

cat("\nCompiling Stan model...\n")

# Note: The Stan model file should be placed in the stan directory
# Users need to ensure the Stan file is available
stan_file <- file.path(stan_dir, "efron_log_spline_model.stan")

if (!file.exists(stan_file)) {
  stop("Stan model file not found. Please ensure 'efron_log_spline_model.stan' ",
       "is in the stan directory.")
}

# Compile model with threading support for parallel chains
compiled_model <- cmdstan_model(
  stan_file,
  cpp_options = list(stan_threads = TRUE),
  stanc_options = list("O1")
)

cat("Model compiled successfully.\n")

# Also compile without threading for optimization methods
compiled_model_no_threads <- cmdstan_model(
  stan_file,
  cpp_options = list(stan_threads = FALSE),
  stanc_options = list("O1")
)

# ==============================================================================
# 5. FIT MODELS FOR SELECTED SCENARIOS
# ==============================================================================

# Select key scenarios for analysis (6 conditions: 3 reliability × 2 heterogeneity)
df2_sub <- df2 %>%
  filter(I %in% c(0.5, 0.7, 0.9)) %>%
  filter(R %in% c(1, 9)) %>%
  mutate(row_id = row_number())

cat("\nSelected scenarios for analysis:\n")
print(df2_sub %>% select(row_id, I, R))

# ==============================================================================
# 6. RUN MCMC SAMPLING
# ==============================================================================

#' Fit Stan Model and Save Results
#'
#' This function runs MCMC sampling for a single scenario and saves
#' the fitted model object.
#'
#' @param stan_data Stan data list
#' @param I_val Reliability value
#' @param R_val Heterogeneity ratio
#' @param row_id Row identifier
#' @param model Compiled Stan model
#' @return Tibble with fitting results and diagnostics

fit_and_save_stan <- function(stan_data, I_val, R_val, row_id, model) {
  # Generate filename for saving
  file_name <- sprintf("%02d_fit_StanFit_I_%02d_R_%d.rds",
                       row_id, as.integer(I_val * 10), R_val)
  
  cat(sprintf("\nFitting scenario %d (I=%.1f, R=%d)...\n", 
              row_id, I_val, R_val))
  
  tryCatch({
    # Start timing
    start_time <- Sys.time()
    
    # Run MCMC sampling
    fit <- model$sample(
      data = stan_data,
      seed = 123 + row_id,      # Unique seed for reproducibility
      chains = 4,                # Number of Markov chains
      parallel_chains = 4,       # Run chains in parallel
      iter_warmup = 1000,        # Warmup iterations
      iter_sampling = 3000,      # Sampling iterations
      refresh = 1000,            # Print progress every 1000 iterations
      adapt_delta = 0.9,         # Target acceptance rate
      threads_per_chain = 1      # Threads per chain
    )
    
    # Save fitted object
    fit$save_object(file = file.path(data_dir, file_name))
    
    # Calculate diagnostics
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    divergences <- sum(fit$sampler_diagnostics()[, , "divergent__"])
    max_rhat <- max(fit$summary()$rhat, na.rm = TRUE)
    
    cat(sprintf("  Completed in %.1f seconds (divergences: %d, max Rhat: %.3f)\n",
                elapsed, divergences, max_rhat))
    
    # Return results
    tibble(
      row_id = row_id,
      I = I_val,
      R = R_val,
      file_name = file_name,
      success = TRUE,
      elapsed_time = elapsed,
      divergences = divergences,
      max_rhat = max_rhat,
      stan_fit = list(fit)
    )
    
  }, error = function(e) {
    cat(sprintf("  Failed: %s\n", e$message))
    tibble(
      row_id = row_id,
      I = I_val,
      R = R_val,
      file_name = file_name,
      success = FALSE,
      elapsed_time = NA_real_,
      divergences = NA_integer_,
      max_rhat = NA_real_,
      stan_fit = list(NULL)
    )
  })
}

# Set up parallel processing
n_cores <- parallel::detectCores()
n_workers <- min(6, floor(n_cores / 4))  # Account for 4 chains per model
cat(sprintf("\nUsing %d parallel workers on %d available cores\n", 
            n_workers, n_cores))

plan(multisession, workers = n_workers)

# Run MCMC fitting for all scenarios
cat("\nRunning MCMC sampling for all scenarios...\n")
tic("Total MCMC fitting time")

results_mcmc <- df2_sub %>%
  mutate(
    fit_result = future_pmap(
      .l = list(stan_data, I, R, row_id),
      .f = ~ fit_and_save_stan(..1, ..2, ..3, ..4, compiled_model),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
  ) %>%
  unnest(fit_result)

toc()

# Reset to sequential processing
plan(sequential)

# ==============================================================================
# 7. EMPIRICAL BAYES WITH DELTA METHOD
# ==============================================================================

#' Empirical Bayes Estimation with Delta Method
#'
#' This function implements the empirical Bayes approach with delta method
#' for uncertainty quantification as described in the paper.
#'
#' @param stan_data Stan data list
#' @param theta_type Type of estimator ("theta_rep", "theta_mean", "theta_map")
#' @param conf_levels Confidence levels for intervals
#' @param stan_mod Compiled RStan model (not CmdStan)
#' @return Tibble with point estimates and confidence intervals

eb_delta_method <- function(stan_data, 
                            theta_type = c("theta_rep", "theta_mean", "theta_map"),
                            conf_levels = c(0.90, 0.95), 
                            stan_mod = NULL) {
  
  # Match theta type argument
  theta_type <- match.arg(theta_type)
  
  # Compile RStan model if not provided
  if (is.null(stan_mod)) {
    stan_mod <- rstan::stan_model(file = stan_file)
  }
  
  tryCatch({
    # 1. Run optimization to get MAP estimates with Hessian
    fit_opt <- rstan::optimizing(
      object = stan_mod, 
      data = stan_data, 
      hessian = TRUE
    )
    
    # 2. Extract data dimensions
    K <- stan_data$K
    theta_hat <- stan_data$theta_hat
    se <- stan_data$sigma
    grid <- stan_data$grid
    Q <- stan_data$B
    
    # 3. Extract estimated parameters
    alpha <- fit_opt$par[grep("^alpha\\[", names(fit_opt$par))]
    
    # Extract appropriate theta values
    theta_pattern <- paste0("^", theta_type, "\\[")
    theta_values <- fit_opt$par[grep(theta_pattern, names(fit_opt$par))]
    
    # 4. Compute covariance matrix of alpha (inverse negative Hessian)
    N_alpha <- length(alpha)
    H_alpha <- fit_opt$hessian[1:N_alpha, 1:N_alpha]
    Sigma_alpha <- -solve(H_alpha)
    
    # 5. Compute mixture weights via softmax
    eta <- as.vector(Q %*% alpha)
    w <- exp(eta - max(eta))  # Numerical stability
    w <- w / sum(w)
    
    # 6. Compute gradient of weights with respect to alpha
    dw_deta <- diag(w) - w %*% t(w)  # Softmax Jacobian
    dw_dalpha <- dw_deta %*% Q
    
    # 7. Initialize results
    n_levels <- length(conf_levels)
    ci_matrix <- matrix(NA, nrow = K, ncol = 2 * n_levels)
    se_vector <- numeric(K)
    
    # 8. Compute delta method variance for each site
    for (k in 1:K) {
      # Posterior weights for site k
      log_post_k <- log(w) + dnorm(theta_hat[k], mean = grid, sd = se[k], log = TRUE)
      w_k <- exp(log_post_k - max(log_post_k))
      w_k <- w_k / sum(w_k)
      
      # Gradient computation
      dw_k_dw <- (diag(as.vector(w_k)) * (1 / w)) - (w_k %*% t(w_k / w))
      
      # Gradient for theta estimator
      if (theta_type %in% c("theta_mean", "theta_rep")) {
        dtheta_k_dw <- t(grid) %*% dw_k_dw
        grad_theta_k <- dtheta_k_dw %*% dw_dalpha
      } else if (theta_type == "theta_map") {
        # Smooth approximation for MAP
        beta <- 100
        w_k_concentrated <- exp(beta * log(w_k))
        w_k_concentrated <- w_k_concentrated / sum(w_k_concentrated)
        dw_k_conc_dw_k <- beta * diag(as.vector(w_k_concentrated))
        dtheta_k_dw <- t(grid) %*% dw_k_conc_dw_k %*% dw_k_dw
        grad_theta_k <- dtheta_k_dw %*% dw_dalpha
      }
      
      # Compute variance and standard error
      var_theta_k <- as.numeric(grad_theta_k %*% Sigma_alpha %*% t(grad_theta_k))
      se_theta_k <- sqrt(var_theta_k)
      se_vector[k] <- se_theta_k
      
      # Compute confidence intervals
      for (p in seq_len(n_levels)) {
        level <- conf_levels[p]
        z_alpha <- qnorm((1 - level) / 2, lower.tail = FALSE)
        ci_matrix[k, (2*p - 1):(2*p)] <- c(
          theta_values[k] - z_alpha * se_theta_k,
          theta_values[k] + z_alpha * se_theta_k
        )
      }
    }
    
    # 9. Create result tibble
    result_tib <- tibble(
      theta_eb = theta_values,
      se_eb = se_vector
    )
    
    # Add confidence interval columns
    for (p in seq_len(n_levels)) {
      level_pct <- conf_levels[p] * 100
      lower_q <- (100 - level_pct) / 2
      upper_q <- 100 - lower_q
      result_tib[[paste0("q", lower_q)]] <- ci_matrix[, 2*p - 1]
      result_tib[[paste0("q", upper_q)]] <- ci_matrix[, 2*p]
    }
    
    attr(result_tib, "theta_type") <- theta_type
    return(result_tib)
    
  }, error = function(e) {
    message("EB Delta Method error: ", e$message)
    return(NULL)
  })
}

# Compile RStan model once for EB
cat("\nCompiling RStan model for Empirical Bayes...\n")
stan_mod_rstan <- rstan::stan_model(file = stan_file)

# Apply EB delta method to all scenarios
cat("Computing Empirical Bayes estimates with delta method...\n")

df3 <- results_mcmc %>%
  mutate(
    eb_theta_ci = map(
      stan_data, 
      ~ eb_delta_method(.x, 
                        theta_type = "theta_rep",
                        conf_levels = c(0.90, 0.95),
                        stan_mod = stan_mod_rstan)
    )
  )

# ==============================================================================
# 8. ALTERNATIVE ESTIMATION METHODS
# ==============================================================================

cat("\nRunning alternative estimation methods...\n")

#' Apply All Estimation Methods
#'
#' This function applies various approximation methods (MLE, MAP, Laplace,
#' Variational Bayes, Pathfinder) to a single dataset.
#'
#' @param stan_data Stan data list
#' @param model Compiled CmdStan model
#' @return List of fitted models

apply_all_methods <- function(stan_data, model) {
  
  # Penalized MLE (without Jacobian adjustment)
  fit_mle <- tryCatch({
    model$optimize(
      data = stan_data,
      algorithm = "lbfgs",
      init = 0,
      jacobian = FALSE
    )
  }, error = function(e) NULL)
  
  # MAP estimate (with Jacobian adjustment)
  fit_map <- tryCatch({
    model$optimize(
      data = stan_data,
      algorithm = "lbfgs",
      init = 0,
      jacobian = TRUE
    )
  }, error = function(e) NULL)
  
  # Laplace approximation
  fit_laplace <- if (!is.null(fit_map)) {
    tryCatch({
      model$laplace(
        data = stan_data,
        mode = fit_map
      )
    }, error = function(e) NULL)
  } else NULL
  
  # Variational Bayes (ADVI)
  fit_vb <- tryCatch({
    model$variational(
      data = stan_data,
      algorithm = "meanfield",
      iter = 10000
    )
  }, error = function(e) NULL)
  
  # Pathfinder
  fit_pf <- tryCatch({
    model$pathfinder(
      data = stan_data
    )
  }, error = function(e) NULL)
  
  # Return all fits as list
  list(
    fit_mle = fit_mle,
    fit_map = fit_map,
    fit_laplace = fit_laplace,
    fit_vb = fit_vb,
    fit_pf = fit_pf
  )
}

# Apply alternative methods to all scenarios
df4 <- df3 %>%
  mutate(
    fit_list = map(
      stan_data,
      ~ apply_all_methods(.x, compiled_model_no_threads)
    )
  )

cat("Alternative estimation methods completed.\n")

# ==============================================================================
# 9. EXTRACT AND ORGANIZE RESULTS
# ==============================================================================

cat("\nExtracting posterior summaries...\n")

#' Extract Theta Summary Statistics
#'
#' This function extracts summary statistics for theta parameters
#' from a fitted Stan model.
#'
#' @param fit CmdStanFit object
#' @return Tibble with summary statistics

extract_theta_summary <- function(fit) {
  tryCatch({
    # Extract theta_rep summary
    df_theta_rep <- fit$summary(
      "theta_rep",
      posterior::default_summary_measures(),
      quantiles = ~ quantile2(., probs = c(0.025, 0.975)),
      posterior::default_convergence_measures()
    )
    
    # Extract theta_map (mode)
    df_theta_map <- fit$summary(
      "theta_map",
      posterior::default_summary_measures()
    ) %>%
      select(mean) %>%
      rename(mode = mean)
    
    # Combine results
    bind_cols(df_theta_rep, df_theta_map) %>%
      relocate(mode, .after = "median")
    
  }, error = function(e) {
    warning(paste("Error extracting theta summary:", e$message))
    return(NULL)
  })
}

# Extract summaries for all MCMC fits
df5 <- df4 %>%
  mutate(
    sum_theta = map(stan_fit, extract_theta_summary)
  )

# ==============================================================================
# 10. SAVE RESULTS
# ==============================================================================

cat("\nSaving results...\n")

# Save complete results
saveRDS(df5, file = file.path(data_dir, "simulation_results_complete.rds"))

# Create summary table of key metrics
summary_table <- df5 %>%
  select(I, R, divergences, max_rhat, elapsed_time) %>%
  mutate(
    convergence = ifelse(max_rhat < 1.01 & divergences == 0, "Good", "Check")
  )

cat("\nSummary of MCMC runs:\n")
print(summary_table)

# ==============================================================================
# 11. SESSION INFORMATION
# ==============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Analysis completed successfully!\n")
cat("\nKey outputs saved:\n")
cat("  - simulation_results_complete.rds\n")
cat("  - Individual Stan fit files (per scenario)\n")
cat("\nSession Information:\n")
print(sessionInfo())