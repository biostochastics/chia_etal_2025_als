# Install limma if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}

cat("limma installation complete\n")
