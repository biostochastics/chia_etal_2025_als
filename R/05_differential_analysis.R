#' Stratified Differential Protein Expression Analysis
#'
#' Functions for comparing ALS vs controls within country strata
#' to identify tube-robust biomarkers.
#'
#' @description
#' Uses limma to perform differential expression analysis:
#' - Italy: ALS vs Healthy controls (HEPARIN tube)
#' - US: ALS vs Healthy + Neurological controls (EDTA tube)
#' - Meta-analysis: Concordance between countries
#'
#' Proteins significant in BOTH countries (concordant effects) are
#' candidates for tube-robust biomarkers.


#' Differential analysis within Italy (HEPARIN cohort)
#'
#' @description
#' Compares ALS vs healthy controls in Italian samples only.
#' All samples use HEPARIN tubes, so this is free from tube type confounding.
#'
#' Uses limma with adjustment for age and sex.
#'
#' @param protein_long Long-format protein data with country labels
#' @return limma results table
#' @examples
#' italy_results <- differential_analysis_italy(protein_long)
#' @export
differential_analysis_italy <- function(protein_long) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("DIFFERENTIAL ANALYSIS: Italy (HEPARIN) - ALS vs Healthy")
  message(paste(rep("=", 70), collapse = ""))

  # Filter to Italy, ALS vs Healthy only
  italy_data <- protein_long %>%
    dplyr::filter(
      country == "Italy",
      Diagnosis %in% c("ALS", "Healthy_control")
    ) %>%
    dplyr::mutate(
      Diagnosis = factor(Diagnosis, levels = c("Healthy_control", "ALS"))
    )

  message(sprintf(
    "\nItaly cohort: %s ALS, %s Healthy",
    sum(italy_data$Diagnosis == "ALS" &
      !duplicated(italy_data$SampleID_deidentified[italy_data$Diagnosis == "ALS"])),
    sum(italy_data$Diagnosis == "Healthy_control" &
      !duplicated(italy_data$SampleID_deidentified[italy_data$Diagnosis == "Healthy_control"]))
  ))

  # Pivot to wide format
  italy_wide <- italy_data %>%
    dplyr::distinct(SampleID_deidentified, Assay, .keep_all = TRUE) %>%
    tidyr::pivot_wider(
      id_cols = c(SampleID_deidentified, Diagnosis, Sex, Age_Collection),
      names_from = Assay,
      values_from = NPX
    ) %>%
    # Remove samples with missing covariates
    tidyr::drop_na(Sex, Age_Collection)

  # Prepare design matrix
  design <- model.matrix(~ Diagnosis + Age_Collection + Sex, data = italy_wide)

  # Extract protein matrix
  protein_mat <- italy_wide %>%
    dplyr::select(-SampleID_deidentified, -Diagnosis, -Sex, -Age_Collection) %>%
    as.matrix() %>%
    t() # Proteins as rows, samples as columns

  # Remove proteins with >20% missing
  missing_pct <- rowMeans(is.na(protein_mat))
  protein_mat <- protein_mat[missing_pct < 0.2, ]

  message(sprintf("Analyzing %s proteins with <20%% missing values", nrow(protein_mat)))

  # Fit limma model
  fit <- limma::lmFit(protein_mat, design)
  fit <- limma::eBayes(fit)

  # Extract results for ALS effect
  results <- limma::topTable(
    fit,
    coef = "DiagnosisALS",
    number = Inf,
    adjust.method = "BH"
  ) %>%
    tibble::rownames_to_column("Assay") %>%
    dplyr::mutate(
      cohort = "Italy",
      tube_type = "HEPARIN"
    )

  # Summary
  n_sig <- sum(results$adj.P.Val < 0.05)
  n_sig_fc <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5)

  message(sprintf("\nSignificant proteins (FDR < 0.05): %s", n_sig))
  message(sprintf("Significant with |logFC| > 0.5: %s", n_sig_fc))

  if (n_sig_fc > 0) {
    message("\nTop 10 proteins by adjusted p-value:")
    top10 <- results %>%
      dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 0.5) %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = 10)
    print(as.data.frame(top10 %>% dplyr::select(Assay, logFC, adj.P.Val)))
  }

  return(results)
}


#' Differential analysis within US (EDTA cohort)
#'
#' @description
#' Compares ALS vs controls in US samples only.
#' All samples use EDTA tubes, so this is free from tube type confounding.
#'
#' NOTE: Small sample size (n=40 ALS) → limited power.
#'
#' @param protein_long Long-format protein data with country labels
#' @return limma results table
#' @examples
#' us_results <- differential_analysis_us(protein_long)
#' @export
differential_analysis_us <- function(protein_long) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("DIFFERENTIAL ANALYSIS: US (EDTA) - ALS vs Controls")
  message(paste(rep("=", 70), collapse = ""))
  message("NOTE: Small ALS sample size (n≈40) → limited statistical power")

  # Filter to US, ALS vs Controls
  us_data <- protein_long %>%
    dplyr::filter(
      country == "US",
      Diagnosis %in% c("ALS", "Healthy_control", "Neurological_control")
    ) %>%
    dplyr::mutate(
      Diagnosis_binary = dplyr::if_else(Diagnosis == "ALS", "ALS", "Control"),
      Diagnosis_binary = factor(Diagnosis_binary, levels = c("Control", "ALS"))
    )

  # Count samples
  n_als <- length(unique(us_data$SampleID_deidentified[us_data$Diagnosis == "ALS"]))
  n_healthy <- length(unique(us_data$SampleID_deidentified[us_data$Diagnosis == "Healthy_control"]))
  n_neuro <- length(unique(us_data$SampleID_deidentified[us_data$Diagnosis == "Neurological_control"]))

  message(sprintf(
    "\nUS cohort: %s ALS, %s Healthy, %s Neurological controls",
    n_als, n_healthy, n_neuro
  ))

  # Pivot to wide format
  us_wide <- us_data %>%
    dplyr::distinct(SampleID_deidentified, Assay, .keep_all = TRUE) %>%
    tidyr::pivot_wider(
      id_cols = c(SampleID_deidentified, Diagnosis_binary, Sex, Age_Collection),
      names_from = Assay,
      values_from = NPX
    ) %>%
    # Remove samples with missing covariates
    tidyr::drop_na(Sex, Age_Collection)

  # Design matrix
  design <- model.matrix(~ Diagnosis_binary + Age_Collection + Sex, data = us_wide)

  # Protein matrix
  protein_mat <- us_wide %>%
    dplyr::select(-SampleID_deidentified, -Diagnosis_binary, -Sex, -Age_Collection) %>%
    as.matrix() %>%
    t()

  # Filter proteins
  missing_pct <- rowMeans(is.na(protein_mat))
  protein_mat <- protein_mat[missing_pct < 0.2, ]

  message(sprintf("Analyzing %s proteins with <20%% missing values", nrow(protein_mat)))

  # Fit limma
  fit <- limma::lmFit(protein_mat, design)
  fit <- limma::eBayes(fit)

  # Extract results
  results <- limma::topTable(
    fit,
    coef = "Diagnosis_binaryALS",
    number = Inf,
    adjust.method = "BH"
  ) %>%
    tibble::rownames_to_column("Assay") %>%
    dplyr::mutate(
      cohort = "US",
      tube_type = "EDTA"
    )

  # Summary
  n_sig <- sum(results$adj.P.Val < 0.05)
  n_sig_fc <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5)

  message(sprintf("\nSignificant proteins (FDR < 0.05): %s", n_sig))
  message(sprintf("Significant with |logFC| > 0.5: %s", n_sig_fc))

  if (n_sig_fc > 0) {
    message("\nTop 10 proteins by adjusted p-value:")
    top10 <- results %>%
      dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 0.5) %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = 10)
    print(as.data.frame(top10 %>% dplyr::select(Assay, logFC, adj.P.Val)))
  }

  return(results)
}


#' Pooled differential analysis (original study approach)
#'
#' @description
#' Differential analysis on pooled Italy + US data.
#' Replicates the original study's approach.
#'
#' IMPORTANT: Does NOT adjust for tube type or country.
#' Results confounded by geographic/tube effects.
#'
#' @param protein_long Long-format protein data
#' @return limma results table
#' @examples
#' pooled_results <- differential_analysis_pooled(protein_long)
#' @export
differential_analysis_pooled <- function(protein_long) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("DIFFERENTIAL ANALYSIS: Pooled (Italy + US)")
  message(paste(rep("=", 70), collapse = ""))
  message("WARNING: Results confounded by country/tube type!")

  # Prepare data
  pooled_data <- protein_long %>%
    dplyr::filter(
      country %in% c("Italy", "US"),
      Diagnosis %in% c("ALS", "Healthy_control", "Neurological_control")
    ) %>%
    dplyr::mutate(
      Diagnosis_binary = dplyr::if_else(Diagnosis == "ALS", "ALS", "Control"),
      Diagnosis_binary = factor(Diagnosis_binary, levels = c("Control", "ALS"))
    )

  # Sample counts
  n_als <- length(unique(pooled_data$SampleID_deidentified[pooled_data$Diagnosis == "ALS"]))
  n_control <- length(unique(pooled_data$SampleID_deidentified[pooled_data$Diagnosis != "ALS"]))

  message(sprintf("\nPooled cohort: %s ALS, %s Controls", n_als, n_control))

  # Pivot to wide
  pooled_wide <- pooled_data %>%
    dplyr::distinct(SampleID_deidentified, Assay, .keep_all = TRUE) %>%
    tidyr::pivot_wider(
      id_cols = c(SampleID_deidentified, Diagnosis_binary, Sex, Age_Collection),
      names_from = Assay,
      values_from = NPX
    ) %>%
    # Remove samples with missing covariates
    tidyr::drop_na(Sex, Age_Collection)

  # Design matrix (no tube type adjustment!)
  design <- model.matrix(~ Diagnosis_binary + Age_Collection + Sex, data = pooled_wide)

  # Protein matrix
  protein_mat <- pooled_wide %>%
    dplyr::select(-SampleID_deidentified, -Diagnosis_binary, -Sex, -Age_Collection) %>%
    as.matrix() %>%
    t()

  # Filter
  missing_pct <- rowMeans(is.na(protein_mat))
  protein_mat <- protein_mat[missing_pct < 0.2, ]

  message(sprintf("Analyzing %s proteins", nrow(protein_mat)))

  # Fit limma
  fit <- limma::lmFit(protein_mat, design)
  fit <- limma::eBayes(fit)

  # Results
  results <- limma::topTable(
    fit,
    coef = "Diagnosis_binaryALS",
    number = Inf,
    adjust.method = "BH"
  ) %>%
    tibble::rownames_to_column("Assay") %>%
    dplyr::mutate(
      cohort = "Pooled",
      tube_type = "Mixed"
    )

  # Summary
  n_sig <- sum(results$adj.P.Val < 0.05)
  n_sig_fc <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5)

  message(sprintf("\nSignificant proteins (FDR < 0.05): %s", n_sig))
  message(sprintf("Significant with |logFC| > 0.5: %s", n_sig_fc))

  return(results)
}


#' Meta-analysis: Test concordance between Italy and US
#'
#' @description
#' Identifies proteins that show consistent effects across both countries.
#' These are candidates for tube-robust biomarkers.
#'
#' Tests:
#' - Directional concordance (same sign of effect)
#' - Significance in both cohorts
#' - Heterogeneity (effect size similarity)
#'
#' @param italy_results Results from Italy analysis
#' @param us_results Results from US analysis
#' @return Concordance analysis table
#' @examples
#' concordance <- test_protein_concordance(italy_results, us_results)
#' @export
test_protein_concordance <- function(italy_results, us_results) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("META-ANALYSIS: Concordance between Italy and US")
  message(paste(rep("=", 70), collapse = ""))

  # Merge results
  merged <- dplyr::inner_join(
    italy_results %>%
      dplyr::select(Assay,
        logFC_italy = logFC, P_italy = P.Value,
        adj_P_italy = adj.P.Val, AveExpr_italy = AveExpr
      ),
    us_results %>%
      dplyr::select(Assay,
        logFC_us = logFC, P_us = P.Value,
        adj_P_us = adj.P.Val, AveExpr_us = AveExpr
      ),
    by = "Assay"
  )

  # Concordance metrics
  merged <- merged %>%
    dplyr::mutate(
      # Directional concordance
      same_direction = sign(logFC_italy) == sign(logFC_us),

      # Significance
      sig_italy = adj_P_italy < 0.05 & abs(logFC_italy) > 0.5,
      sig_us = adj_P_us < 0.05 & abs(logFC_us) > 0.5,
      sig_both = sig_italy & sig_us,

      # Effect size correlation
      logFC_diff = abs(logFC_italy - logFC_us),

      # Combined p-value (Fisher's method)
      combined_p = pchisq(
        -2 * (log(P_italy) + log(P_us)),
        df = 4,
        lower.tail = FALSE
      ),

      # Mean effect size
      mean_logFC = (logFC_italy + logFC_us) / 2
    ) %>%
    dplyr::arrange(combined_p)

  # Summary statistics
  n_total <- nrow(merged)
  n_concordant <- sum(merged$same_direction)
  n_sig_both <- sum(merged$sig_both)
  n_sig_both_concordant <- sum(merged$sig_both & merged$same_direction)

  message(sprintf("\nTotal proteins compared: %s", n_total))
  message(sprintf(
    "Directionally concordant: %s (%.1f%%)",
    n_concordant, 100 * n_concordant / n_total
  ))
  message(sprintf("Significant in BOTH cohorts: %s", n_sig_both))
  message(sprintf("Significant in both + concordant: %s", n_sig_both_concordant))

  if (n_sig_both_concordant > 0) {
    message("\n", paste(rep("-", 70), collapse = ""))
    message("TUBE-ROBUST BIOMARKER CANDIDATES:")
    message("(Significant in both Italy and US with same direction)")
    message(paste(rep("-", 70), collapse = ""))

    candidates <- merged %>%
      dplyr::filter(sig_both, same_direction) %>%
      dplyr::select(
        Assay, logFC_italy, logFC_us, mean_logFC,
        adj_P_italy, adj_P_us, combined_p
      ) %>%
      dplyr::arrange(combined_p)

    print(as.data.frame(candidates))
  } else {
    message("\n⚠️  NO proteins are significant in BOTH cohorts with concordant direction.")
    message("This suggests limited tube-robust biomarker signal.")
  }

  # Effect size correlation
  cor_logFC <- cor(merged$logFC_italy, merged$logFC_us, use = "complete.obs")
  message(sprintf("\nEffect size correlation (Italy vs US): r = %.3f", cor_logFC))

  return(merged)
}


#' Compare stratified vs pooled differential analysis
#'
#' @description
#' THE PRIMARY EVIDENCE for differential analysis:
#' Compares protein lists from stratified (within-country) vs pooled analysis.
#'
#' Large discrepancies indicate pooled results are confounded by country/tube.
#'
#' @param italy_results Italy stratified results
#' @param us_results US stratified results
#' @param pooled_results Pooled analysis results
#' @return Comparison summary
#' @examples
#' comparison <- compare_stratified_vs_pooled(italy_results, us_results, pooled_results)
#' @export
compare_stratified_vs_pooled <- function(italy_results, us_results, pooled_results) {
  message("\n\n", paste(rep("=", 80), collapse = ""))
  message("CRITICAL COMPARISON: Stratified vs Pooled Differential Analysis")
  message(paste(rep("=", 80), collapse = ""))

  # Define significance threshold
  sig_threshold_p <- 0.05
  sig_threshold_fc <- 0.5

  # Significant proteins in each analysis
  sig_italy <- italy_results %>%
    dplyr::filter(adj.P.Val < sig_threshold_p, abs(logFC) > sig_threshold_fc) %>%
    dplyr::pull(Assay)

  sig_us <- us_results %>%
    dplyr::filter(adj.P.Val < sig_threshold_p, abs(logFC) > sig_threshold_fc) %>%
    dplyr::pull(Assay)

  sig_pooled <- pooled_results %>%
    dplyr::filter(adj.P.Val < sig_threshold_p, abs(logFC) > sig_threshold_fc) %>%
    dplyr::pull(Assay)

  sig_both_strata <- intersect(sig_italy, sig_us)

  message("\nSignificant proteins (FDR < 0.05, |logFC| > 0.5):")
  message(sprintf("  Italy (HEPARIN):           %s", length(sig_italy)))
  message(sprintf("  US (EDTA):                 %s", length(sig_us)))
  message(sprintf("  Pooled (confounded):       %s", length(sig_pooled)))
  message(sprintf("  Significant in BOTH strata: %s", length(sig_both_strata)))

  # Overlap analysis
  pooled_only <- setdiff(sig_pooled, union(sig_italy, sig_us))
  pooled_overlap_italy <- intersect(sig_pooled, sig_italy)
  pooled_overlap_us <- intersect(sig_pooled, sig_us)
  pooled_overlap_both <- intersect(sig_pooled, sig_both_strata)

  message("\n", paste(rep("-", 80), collapse = ""))
  message("OVERLAP ANALYSIS:")
  message(sprintf(
    "  Pooled ∩ Italy:            %s (%.1f%% of pooled)",
    length(pooled_overlap_italy),
    100 * length(pooled_overlap_italy) / max(1, length(sig_pooled))
  ))
  message(sprintf(
    "  Pooled ∩ US:               %s (%.1f%% of pooled)",
    length(pooled_overlap_us),
    100 * length(pooled_overlap_us) / max(1, length(sig_pooled))
  ))
  message(sprintf(
    "  Pooled ∩ (Italy AND US):   %s (%.1f%% of pooled)",
    length(pooled_overlap_both),
    100 * length(pooled_overlap_both) / max(1, length(sig_pooled))
  ))
  message(sprintf("  Pooled ONLY (not in either stratum): %s", length(pooled_only)))

  # Interpretation
  message("\n", paste(rep("-", 80), collapse = ""))
  message("INTERPRETATION:")

  pct_pooled_in_both <- 100 * length(pooled_overlap_both) / max(1, length(sig_pooled))

  if (pct_pooled_in_both < 30) {
    message("  ⚠️  CRITICAL: <30% of pooled proteins replicate in BOTH strata")
    message("  Pooled results are HEAVILY CONFOUNDED by country/tube effects.")
    message("  Most 'significant' proteins are likely false discoveries.")
  } else if (pct_pooled_in_both < 50) {
    message("  ⚠️  SUBSTANTIAL: 30-50% replication in both strata")
    message("  Major confounding present in pooled analysis.")
  } else if (pct_pooled_in_both < 70) {
    message("  ⚠️  MODERATE: 50-70% replication in both strata")
    message("  Some confounding effects present.")
  } else {
    message("  ✓ Good replication: >70% of pooled proteins in both strata")
    message("  Pooled results are reasonably robust.")
  }

  if (length(pooled_only) > 0.2 * length(sig_pooled)) {
    message(sprintf(
      "\n  ⚠️  WARNING: %s proteins significant in pooled but NOT in either stratum!",
      length(pooled_only)
    ))
    message("  These are likely driven by country/tube confounding.")
  }

  message(paste(rep("=", 80), collapse = ""))

  return(list(
    sig_italy = sig_italy,
    sig_us = sig_us,
    sig_pooled = sig_pooled,
    sig_both_strata = sig_both_strata,
    pooled_only = pooled_only,
    overlap_italy = pooled_overlap_italy,
    overlap_us = pooled_overlap_us,
    overlap_both = pooled_overlap_both,
    pct_replicate_in_both = pct_pooled_in_both
  ))
}
