################################################################################
#
# Replication Code for:
# "Fully Bayesian Inference for Meta-Analytic Deconvolution 
#  Using Efron's Log-Spline Prior"
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: January 2025
#
# Description:
# This script generates simulated data and produces Figures 1-2 from the paper.
# The simulation creates heteroscedastic and homoscedastic scenarios with 
# varying reliability levels to test the performance of Bayesian deconvolution.
#
################################################################################

# ==============================================================================
# 1. SETUP AND PACKAGE LOADING
# ==============================================================================

# Clear workspace and memory
rm(list = ls(all.names = TRUE))
gc()

# Load required packages
# Note: Install packages if not already available
required_packages <- c("tidyverse", "deconvolveR", "ggplot2", "dplyr", 
                       "tidyr", "purrr", "tibble", "viridis")

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

# Create data directory if it doesn't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

setwd(work_dir)

# ==============================================================================
# 2. DATA SIMULATION FUNCTIONS
# ==============================================================================

#' Simulate Site-Specific Observations Under Heteroscedastic or Homoscedastic SEs
#'
#' This function generates observed site-specific effects from true effects
#' under specified reliability (I) and heterogeneity (R) conditions.
#'
#' @param theta_true Numeric vector of true site effects (length K)
#' @param I Target reliability in (0,1). Higher values indicate more informative data
#' @param R Heterogeneity ratio. R = 1 for homoscedastic; R > 1 for heteroscedastic
#' @param shuffle Logical. If TRUE (default), randomly permute the variances across sites
#' @param seed Optional RNG seed for reproducibility
#' @return Data frame with columns: theta_true, se2, se, theta_hat, I_level, R_level

simulate_theta_hat <- function(theta_true, I, R = 1, shuffle = TRUE, seed = NULL) {
  
  # Input validation
  if (!is.null(seed)) set.seed(seed)
  if (I <= 0 || I >= 1) stop("I must be between 0 and 1 (exclusive).")
  if (R <= 0) stop("R must be positive.")
  
  K <- length(theta_true)
  
  # Calculate geometric mean of SE^2 implied by reliability I
  # Based on: I = 1 / (1 + GM(se^2))
  gm_se2 <- (1 - I) / I
  
  if (abs(R - 1) < 1e-8) {
    # Homoscedastic case: all sites have same variance
    se2 <- rep(gm_se2, K)
  } else {
    # Heteroscedastic case: variances vary across sites
    se2_max <- gm_se2 * R  # Maximum variance
    se2_min <- gm_se2 / R  # Minimum variance
    
    # Generate K variances equally spaced on log-scale
    se2 <- exp(seq(log(se2_min), log(se2_max), length.out = K))
    
    # Randomly assign variances to sites if shuffle = TRUE
    if (shuffle) se2 <- sample(se2)
  }
  
  # Calculate standard errors and generate observed effects
  se <- sqrt(se2)
  theta_hat <- rnorm(K, mean = theta_true, sd = se)
  
  # Return results as tibble
  tibble(
    theta_true = theta_true,
    se2 = se2,
    se = se,
    theta_hat = theta_hat,
    I_level = I,
    R_level = R
  )
}

# ==============================================================================
# 3. LOAD TRUE EFFECTS DATA
# ==============================================================================

# Load the "Twin Towers" bimodal distribution from deconvolveR package
# This represents a challenging test case with two distinct modes
data(disjointTheta)
theta_true <- disjointTheta
K <- length(theta_true)

cat("Dataset loaded successfully:\n")
cat("  Number of sites: ", K, "\n")
cat("  Range of true effects: [", round(min(theta_true), 2), 
    ", ", round(max(theta_true), 2), "]\n")
cat("  Variance of true effects: ", round(var(theta_true), 3), "\n\n")

# ==============================================================================
# 4. GENERATE SIMULATION SCENARIOS
# ==============================================================================

# Define simulation parameters
# I: Average reliability (signal-to-noise ratio)
# R: Heterogeneity ratio (1 = homoscedastic, >1 = heteroscedastic)
I_levels <- c(0.5, 0.6, 0.7, 0.8, 0.9)
R_levels <- c(1, 5, 9)

# Set global seed for reproducibility
set.seed(12345)

# Create all combinations of I and R
sim_conditions <- expand_grid(I = I_levels, R = R_levels)

cat("Simulation design:\n")
cat("  I levels (reliability): ", paste(I_levels, collapse = ", "), "\n")
cat("  R levels (heterogeneity): ", paste(R_levels, collapse = ", "), "\n")
cat("  Total scenarios: ", nrow(sim_conditions), "\n\n")

# Generate simulated datasets for each combination
sim_data <- sim_conditions %>%
  mutate(
    data = pmap(
      list(I, R),
      ~ simulate_theta_hat(
        theta_true = theta_true,
        I = ..1,
        R = ..2,
        shuffle = TRUE
      )
    )
  )

# Display summary of generated data
cat("Simulation data generated successfully:\n")
print(sim_data %>% select(I, R))

# ==============================================================================
# 5. CREATE FIGURE 1: DENSITY COMPARISONS
# ==============================================================================

cat("\nGenerating Figure 1: True vs. Observed Site-Specific Effects...\n")

# Select conditions for visualization (6 scenarios for 3x2 grid)
sel_conditions <- tibble(
  I = c(0.5, 0.5, 0.7, 0.7, 0.9, 0.9),
  R = c(1, 9, 1, 9, 1, 9)
)

# Prepare data for plotting
sim_sel <- sim_data %>%
  right_join(sel_conditions, by = c("I", "R")) %>%
  unnest(data) %>%
  group_by(I, R) %>%
  mutate(
    mean_se = sqrt(exp(mean(log(se2)))),  # Geometric mean of SE
    panel_lab = paste0("I = ", I, ", R = ", R, "\nMean(SE) = ", round(mean_se, 2)),
    hetero_lab = if_else(R == 1, "Homoscedastic", "Heteroscedastic (R = 9)")
  ) %>%
  ungroup()

# Create Figure 1
bin_width <- 0.25

fig1 <- ggplot(sim_sel, aes(x = theta_true)) +
  # True distribution (red)
  geom_histogram(
    aes(y = after_stat(density)),
    fill = "firebrick", 
    colour = "firebrick4",
    alpha = 0.4, 
    binwidth = bin_width
  ) +
  # Observed distribution (gray)
  geom_histogram(
    aes(x = theta_hat, y = after_stat(density)),
    fill = "grey70", 
    colour = "grey50",
    alpha = 0.35, 
    binwidth = bin_width
  ) +
  # Density curve for observed effects
  geom_density(
    aes(x = theta_hat), 
    colour = "black", 
    size = 0.7
  ) +
  # Facet by reliability and heteroscedasticity
  facet_grid(
    rows = vars(I), 
    cols = vars(hetero_lab),
    labeller = labeller(I = as_labeller(function(x) paste("I =", x)))
  ) +
  labs(
    x = expression(theta),
    y = "Density",
    title = "Figure 1. True vs. Observed Site-Specific Effects\nacross Informativeness (I) and Heteroscedasticity (R) Conditions"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  # Add panel labels with Mean(SE)
  geom_text(
    data = sim_sel %>% distinct(I, hetero_lab, panel_lab),
    aes(label = panel_lab),
    x = -Inf, y = Inf, 
    hjust = -0.05, vjust = 1.15,
    size = 3.3, 
    colour = "black"
  )

print(fig1)

# Save Figure 1
ggsave(
  filename = file.path(data_dir, "Figure_01_Simulation_Scenario.pdf"),
  plot = fig1,
  width = 8,
  height = 7,
  dpi = 300
)

cat("Figure 1 saved successfully.\n")

# ==============================================================================
# 6. CREATE FIGURE 2: SCATTERPLOT WITH HETEROSCEDASTICITY
# ==============================================================================

cat("\nGenerating Figure 2: True vs. Observed Effects with Error Heterogeneity...\n")

# Create scatterplot with color-coded standard errors
fig2 <- ggplot(sim_sel, aes(x = theta_true, y = theta_hat)) +
  # Add diagonal reference line
  geom_abline(
    intercept = 0, 
    slope = 1,
    linetype = "dashed", 
    color = "black", 
    alpha = 0.7
  ) +
  # Points colored by standard error
  geom_point(
    aes(color = se), 
    alpha = 0.6, 
    size = 1
  ) +
  # Color scale for SE
  scale_color_viridis_c(
    option = "plasma",
    name = "Standard\nError",
    limits = c(0, 3)
  ) +
  # Facet by reliability and heteroscedasticity
  facet_grid(
    rows = vars(I), 
    cols = vars(hetero_lab),
    labeller = labeller(I = as_labeller(function(x) paste("I =", x)))
  ) +
  labs(
    x = expression("True"~theta),
    y = expression("Observed"~hat(theta)),
    title = "Figure 2. True vs Observed Effects: Color = Standard Error",
    subtitle = "Color intensity shows measurement error heterogeneity"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "right"
  ) +
  # Add panel labels
  geom_text(
    data = sim_sel %>% distinct(I, hetero_lab, panel_lab),
    aes(label = panel_lab),
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.15,
    size = 3.3,
    colour = "black"
  )

print(fig2)

# Save Figure 2
ggsave(
  filename = file.path(data_dir, "Figure_02_Simulation_Scenario_Scatterplot.pdf"),
  plot = fig2,
  width = 8,
  height = 7,
  dpi = 300
)

cat("Figure 2 saved successfully.\n")

# ==============================================================================
# 7. SAVE SIMULATION DATA
# ==============================================================================

cat("\nSaving simulation data...\n")

# Save the complete simulation data
saveRDS(
  sim_data,
  file = file.path(data_dir, "simulated_data_varied_by_I_and_R.rds")
)

cat("Data saved successfully to:", file.path(data_dir, "simulated_data_varied_by_I_and_R.rds"), "\n")

# ==============================================================================
# 8. SESSION INFORMATION
# ==============================================================================

cat("\n" , paste(rep("=", 80), collapse = ""), "\n")
cat("Simulation completed successfully!\n")
cat("Generated files:\n")
cat("  - Figure_01_Simulation_Scenario.pdf\n")
cat("  - Figure_02_Simulation_Scenario_Scatterplot.pdf\n")
cat("  - simulated_data_varied_by_I_and_R.rds\n")
cat("\nSession Information:\n")
print(sessionInfo())
