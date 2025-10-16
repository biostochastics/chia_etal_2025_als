#' Tube Type Effect Analysis with Limma
#'
#' Calculate tube type (EDTA vs HEPARIN) effects within diagnosis groups
#' using limma with age and sex adjustment.
#'
#' @description
#' This function tests for systematic differences between Italy (HEPARIN)
#' and US (EDTA) samples within a specific diagnosis group, adjusting for
#' age and sex confounders using limma.
#'
#' @param protein_long Long-format protein data with country labels
#' @param diagnosis_group Character: "ALS" or "Healthy_control"
#' @return limma results table with tube effect statistics
#' @export
calculate_tube_effects_limma <- function(protein_long, diagnosis_group) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message(sprintf("TUBE EFFECTS ANALYSIS: %s patients", diagnosis_group))
  message(paste(rep("=", 70), collapse = ""))
  message("Testing country/tube differences adjusting for age and sex")

  # Filter to specific diagnosis and relevant countries
  diagnosis_data <- protein_long %>%
    dplyr::filter(
      Diagnosis == diagnosis_group,
      country %in% c("Italy", "US")
    ) %>%
    dplyr::mutate(
      country = factor(country, levels = c("Italy", "US"))  # Italy (HEPARIN) as reference
    )

  # Count samples
  n_italy <- length(unique(diagnosis_data$SampleID_deidentified[diagnosis_data$country == "Italy"]))
  n_us <- length(unique(diagnosis_data$SampleID_deidentified[diagnosis_data$country == "US"]))

  message(sprintf("\n%s cohort: %s Italy (HEPARIN), %s US (EDTA)",
    diagnosis_group, n_italy, n_us
  ))

  # Pivot to wide format
  diagnosis_wide <- diagnosis_data %>%
    dplyr::distinct(SampleID_deidentified, Assay, .keep_all = TRUE) %>%
    tidyr::pivot_wider(
      id_cols = c(SampleID_deidentified, country, Sex, Age_Collection),
      names_from = Assay,
      values_from = NPX
    ) %>%
    # Remove samples with missing covariates
    tidyr::drop_na(Sex, Age_Collection)

  # Prepare design matrix: country effect adjusted for age and sex
  design <- model.matrix(~ country + Age_Collection + Sex, data = diagnosis_wide)

  # Extract protein matrix
  protein_mat <- diagnosis_wide %>%
    dplyr::select(-SampleID_deidentified, -country, -Sex, -Age_Collection) %>%
    as.matrix() %>%
    t() # Proteins as rows, samples as columns

  # Remove proteins with >20% missing
  missing_pct <- rowMeans(is.na(protein_mat))
  protein_mat <- protein_mat[missing_pct < 0.2, ]

  message(sprintf("Analyzing %s proteins with <20%% missing values", nrow(protein_mat)))

  # Fit limma model
  fit <- limma::lmFit(protein_mat, design)
  fit <- limma::eBayes(fit)

  # Extract results for country effect (US vs Italy = EDTA vs HEPARIN)
  results <- limma::topTable(
    fit,
    coef = "countryUS",
    number = Inf,
    adjust.method = "BH"
  ) %>%
    tibble::rownames_to_column("Assay") %>%
    dplyr::mutate(
      diagnosis_group = diagnosis_group,
      tube_comparison = "US/EDTA vs Italy/HEPARIN"
    )

  # Summary
  n_sig <- sum(results$adj.P.Val < 0.05)
  n_sig_fc <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5)

  message(sprintf("\nProteins with significant tube effect (FDR < 0.05): %s", n_sig))
  message(sprintf("Significant with |logFC| > 0.5: %s", n_sig_fc))

  if (n_sig > 0) {
    message("\nTop 10 proteins with strongest tube effects:")
    top10 <- results %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = 10)
    print(as.data.frame(top10 %>% dplyr::select(Assay, logFC, adj.P.Val)))
  }

  message(paste(rep("=", 70), collapse = ""))

  return(results)
}
