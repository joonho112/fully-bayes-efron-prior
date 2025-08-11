# Replication Package for "Fully Bayesian Inference for Meta-Analytic Deconvolution Using Efron's Log-Spline Prior"

## Overview

This repository contains the complete replication package for the paper:

**Title:** Fully Bayesian Inference for Meta-Analytic Deconvolution Using Efron's Log-Spline Prior  
**Authors:** JoonHo Lee and Daihe Sui  
**Journal:** Mathematics  
**Year:** 2025  

This package enables full reproduction of all figures, tables, and results presented in the main manuscript and appendices.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Installation Guide](#installation-guide)
3. [Repository Structure](#repository-structure)
4. [Data Description](#data-description)
5. [Reproduction Instructions](#reproduction-instructions)
6. [Expected Outputs](#expected-outputs)
7. [Troubleshooting](#troubleshooting)
8. [Contact Information](#contact-information)

## System Requirements

### Software Dependencies
- **R** (≥ 4.2.0)
- **Stan** (via CmdStan interface)
- **Operating System:** Windows 10+, macOS 10.15+, or Linux

### Hardware Requirements
- **RAM:** Minimum 8GB (16GB recommended for full MCMC analysis)
- **CPU:** Multi-core processor recommended (4+ cores for parallel processing)
- **Storage:** ~2GB free space for outputs

### Estimated Runtime
- Full replication: 4-6 hours
- Quick replication (reduced iterations): 1-2 hours

## Installation Guide

### Step 1: Clone or Download Repository
```bash
git clone https://github.com/[username]/efron-log-spline-replication.git
cd efron-log-spline-replication
```

### Step 2: Install Required R Packages
Open R and run:
```r
# Install CRAN packages
required_packages <- c(
  "tidyverse", "splines", "patchwork", "scales", "viridis",
  "knitr", "kableExtra", "ggrepel", "posterior", "bayesplot",
  "fixest", "haven", "deconvolveR", "furrr", "progressr"
)

install.packages(required_packages)

# Install CmdStan interface
install.packages("cmdstanr", 
                 repos = c("https://mc-stan.org/r-packages/", 
                          getOption("repos")))
```

### Step 3: Install CmdStan
```r
library(cmdstanr)
install_cmdstan(cores = 4)
```

This will download and install CmdStan automatically (~200MB download).

## Repository Structure

```
efron-log-spline-replication/
│
├── README.md                           # This file
├── stan/                              # Stan model files
│   └── efron_log_spline_model.stan   # Main Stan model
│
├── data/                              # Data files (provided)
│   ├── simulated_data_varied_by_I_and_R.rds
│   ├── krw_data.rds                  # Real-world discrimination data
│   └── firm_discrimination_estimates.rds
│
├── code/                              # R scripts for analysis
│   ├── Part_01_Data_Generation.R
│   ├── Part_02_Model_Estimation_and_Performance_Evaluation.R
│   ├── Part_03_Figures_and_Tables_for_Publication.R
│   ├── Part_04_Grid_Resolution_and_Bounds_Sensitivity_Analysis.R
│   ├── Part_05_Lambda_Prior_Sensitivity_Analysis.R
│   ├── Part_06_Small-K_Scenarios_Analysis.R
│   └── Part_07_Real-World_Application.R
│
├── figures/                           # Output directory for figures (created)
├── results/                           # Output directory for results (created)
├── appendix_c/                        # Appendix C outputs (created)
├── appendix_d/                        # Appendix D outputs (created)
└── appendix_e/                        # Appendix E outputs (created)
```

## Data Description

### 1. `simulated_data_varied_by_I_and_R.rds`
- **Content:** Simulated "twin towers" bimodal datasets
- **Structure:** Tibble with columns:
  - `I`: Reliability parameter (0.5, 0.6, 0.7, 0.8, 0.9)
  - `R`: Heterogeneity ratio (1, 5, 9)
  - `data`: Nested data frames containing theta_true, theta_hat, se
- **Size:** 1.1MB
- **Usage:** Main simulation studies (Sections 3-4)

### 2. `krw_data.rds`
- **Content:** Firm-level labor market discrimination audit study data
- **Source:** Kline, Rose, and Walters (2024)
- **Variables:** firm_id, job_id, callback, white, male, firm_name
- **Usage:** Real-world application (Section 5)
- **Note:** If unavailable, code generates simulated example data

### 3. `firm_discrimination_estimates.rds`
- **Content:** Pre-computed firm-level discrimination estimates
- **Structure:** Firm-level OLS estimates with standard errors
- **Usage:** Alternative starting point for Section 5 analysis

## Reproduction Instructions

### Quick Start (Essential Results Only)

To reproduce the main results quickly (~30 minutes):

```r
# Set working directory to repository root
setwd("path/to/efron-log-spline-replication")

# Run essential analyses
source("code/Part_01_Data_Generation.R")           # Figures 1-2
source("code/Part_03_Figures_and_Tables_for_Publication.R")  # Figures 3-9, Tables 1-2
source("code/Part_07_Real-World_Application.R")     # Figures 10-12, Table 4
```

### Complete Replication (All Results)

For full reproduction including all appendices:

#### Step 1: Generate Simulation Data
```r
source("code/Part_01_Data_Generation.R")
```
**Outputs:**
- `Figure_01_Simulation_Scenario.pdf`
- `Figure_02_Simulation_Scenario_Scatterplot.pdf`
- `simulated_data_varied_by_I_and_R.rds`

**Runtime:** ~5 minutes

#### Step 2: Fit Models and Evaluate Performance
```r
source("code/Part_02_Model_Estimation_and_Performance_Evaluation.R")
```
**Outputs:**
- Stan model fits for each scenario
- `simulation_results_complete.rds`

**Runtime:** ~2-3 hours (MCMC sampling)

**Note:** To reduce runtime, modify the MCMC settings:
```r
# In Part_02, line ~360, change:
iter_warmup = 500,      # Reduced from 1000
iter_sampling = 1000,   # Reduced from 3000
```

#### Step 3: Generate Main Figures and Tables
```r
source("code/Part_03_Figures_and_Tables_for_Publication.R")
```
**Outputs:**
- Figures 3-9 (PDF format)
- Tables 1-2 (CSV format)

**Runtime:** ~10 minutes

#### Step 4: Appendix C - Grid Sensitivity Analysis
```r
source("code/Part_04_Grid_Resolution_and_Bounds_Sensitivity_Analysis.R")
```
**Outputs in `appendix_c/`:**
- `Figure_C1_Grid_Sensitivity.pdf`
- `Figure_C2_Site_Errors.pdf`
- `Table_C1_Grid_Sensitivity_Metrics.csv`
- `Table_C2_Grid_Sensitivity_Comparison.csv`

**Runtime:** ~20 minutes

#### Step 5: Appendix D - Lambda Prior Sensitivity
```r
source("code/Part_05_Lambda_Prior_Sensitivity_Analysis.R")
```
**Outputs in `appendix_d/`:**
- `Figure_D1_Lambda_Prior_Distributions.pdf`
- `Figure_D2_Lambda_Sensitivity_Performance.pdf`
- `Figure_D3_Lambda_Sensitivity_Robustness.pdf`
- Tables D1-D2 (CSV format)

**Runtime:** ~30 minutes

#### Step 6: Appendix E - Small-K Scenarios
```r
source("code/Part_06_Small-K_Scenarios_Analysis.R")
```
**Outputs in `appendix_e/`:**
- `Figure_E1_Small_K_Performance.pdf`
- `Figure_E2_K50_Detail.pdf`
- Tables E1-E2 (CSV format)

**Runtime:** ~45 minutes

#### Step 7: Real-World Application
```r
source("code/Part_07_Real-World_Application.R")
```
**Outputs in `results/`:**
- `Figure_10_Caterpillar_Race.pdf`
- `Figure_11_Caterpillar_Gender.pdf`
- `Figure_12_Density_Comparison.pdf`
- `Table_4_Summary_Statistics.csv`

**Runtime:** ~15 minutes

## Expected Outputs

### Main Manuscript
| Figure/Table | Description | File |
|-------------|-------------|------|
| Figure 1 | Simulation scenarios | `Figure_01_Simulation_Scenario.pdf` |
| Figure 2 | True vs observed effects | `Figure_02_Simulation_Scenario_Scatterplot.pdf` |
| Figures 3-9 | Performance evaluation | `Figure_03-09_*.pdf` |
| Figures 10-12 | Real-world application | `Figure_10-12_*.pdf` |
| Table 1 | Performance metrics | `Table_1_Performance_Metrics.csv` |
| Table 2 | Coverage rates | `Table_2_Coverage_Rates.csv` |
| Table 4 | Discrimination summary | `Table_4_Summary_Statistics.csv` |

### Appendices
| Appendix | Content | Number of Outputs |
|----------|---------|-------------------|
| C | Grid sensitivity | 2 figures, 2 tables |
| D | Lambda prior sensitivity | 3 figures, 2 tables |
| E | Small-K scenarios | 2 figures, 2 tables |

## Troubleshooting

### Common Issues and Solutions

#### 1. Stan Compilation Error
```
Error: Stan model does not compile
```
**Solution:** Ensure CmdStan is properly installed:
```r
cmdstanr::cmdstan_path()  # Should return path
cmdstanr::rebuild_cmdstan()  # Rebuild if needed
```

#### 2. Memory Issues
```
Error: cannot allocate vector of size...
```
**Solution:** 
- Close other applications
- Reduce parallel workers: `plan(multisession, workers = 2)`
- Reduce MCMC iterations (see Quick Start options)

#### 3. Missing Data File
```
Error: krw_data.rds not found
```
**Solution:** The code will automatically generate simulated example data. To use real data, ensure `krw_data.rds` is in the `data/` folder.

#### 4. Package Installation Issues
**Solution:** Install packages individually:
```r
# If tidyverse fails, install components:
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble"))
```

#### 5. Long Runtime
To reduce computation time:
- Skip Part_02 and use pre-computed results (if provided)
- Reduce VB iterations: Change `iter = 10000` to `iter = 5000`
- Run only selected scenarios instead of all

### Platform-Specific Notes

**Windows Users:**
- Use forward slashes in paths: `C:/Users/...` not `C:\Users\...`
- May need to install Rtools for Stan compilation

**Mac Users:**
- May need to install Xcode Command Line Tools
- Run in Terminal: `xcode-select --install`

**Linux Users:**
- Ensure g++ compiler is installed: `sudo apt-get install g++`

## Reproducibility Guarantee

This package has been tested on:
- Windows 11 (R 4.3.1, CmdStan 2.33.1)
- macOS Ventura 13.5 (R 4.3.0, CmdStan 2.33.1)
- Ubuntu 22.04 LTS (R 4.2.3, CmdStan 2.33.1)

Results may vary slightly due to:
- Different random number generators across platforms
- Floating-point arithmetic differences
- MCMC sampling variability

Expected variation: <1% in point estimates, <2% in coverage rates

## Citation

If you use this code or methodology, please cite:

```bibtex
@article{lee2025fully,
  title={Fully Bayesian Inference for Meta-Analytic Deconvolution Using Efron's Log-Spline Prior},
  author={Lee, JoonHo and Sui, Daihe},
  journal={Mathematics},
  year={2025},
  volume={XX},
  pages={XXX--XXX}
}
```

## Contact Information

**Corresponding Author:**  
JoonHo Lee  
Email: jlee296@ua.edu  
The University of Alabama  

For questions or bug reports, please:
1. Check the Troubleshooting section
2. Open an issue on GitHub
3. Contact the corresponding author

## License

This replication package is released under the MIT License. See LICENSE file for details.

## Acknowledgments

The authors would like to thank the Institute of Education Sciences for their support of this research. We are also grateful to the developers of the `deconvolveR` R package for making their simulated data publicly available, which facilitated the comparative simulation analysis presented in this work. We thank Dr. Christopher Walters for generously providing public access to the replication data from Kline et al. (2022, 2024), which enabled the real-data application examining firm-level labor market discrimination.

---

Last Updated: 2025-08-10
Version: 1.0.0