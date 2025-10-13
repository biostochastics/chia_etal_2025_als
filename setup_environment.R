#!/usr/bin/env Rscript
# Setup script for ALS biomarker investigation project
# Initializes renv and installs required packages

cat("============================================\n")
cat("ALS Biomarker Investigation - Setup\n")
cat("============================================\n\n")

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  cat("Installing renv...\n")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialize renv (bare = TRUE for clean start)
cat("\nInitializing renv environment...\n")
renv::init(bare = TRUE, restart = FALSE)

# Define required packages
packages <- c(
  # Data manipulation
  "tidyverse",
  "data.table",
  "dtplyr",

  # Pipeline management
  "targets",
  "tarchetypes",

  # Statistical analysis
  "broom",
  "DescTools", # For Cramér's V
  "vcd", # For association measures

  # Machine learning
  "caret",
  "ranger", # Random forest
  "glmnet", # Regularized regression
  "DALEX", # Model explainability
  "pROC", # ROC curves

  # Visualization
  "ggplot2",
  "patchwork",
  "ggpubr",
  "corrplot",

  # Reporting
  "knitr",
  "rmarkdown",

  # Testing
  "testthat",
  "covr"
)

# Install packages
cat("\nInstalling required packages...\n")
cat("This may take several minutes...\n\n")

for (pkg in packages) {
  cat(sprintf("Installing %s...\n", pkg))
  renv::install(pkg)
}

# Take snapshot
cat("\nCreating renv.lock snapshot...\n")
renv::snapshot(prompt = FALSE)

cat("\n============================================\n")
cat("Setup complete!\n")
cat("============================================\n")
cat("\nNext steps:\n")
cat("1. Source this file: source('setup_environment.R')\n")
cat("2. Run targets pipeline: targets::tar_make()\n")
cat("3. View pipeline: targets::tar_visnetwork()\n\n")
