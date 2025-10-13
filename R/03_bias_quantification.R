#' Confounding Bias Quantification
#'
#' Functions for quantifying the degree of confounding between tube type,
#' country, and diagnosis in the ALS biomarker study.
#'
#' @description
#' This module implements:
#' - Cramér's V for association strength
#' - Chi-square tests for independence
#' - Contingency tables for visualization
#' - Stratified analyses
#'

#' Calculate Cramér's V for two categorical variables
#'
#' @description
#' Cramér's V measures association strength between categorical variables:
#' - 0 = no association
#' - 0.1 = weak
#' - 0.3 = moderate
#' - 0.5+ = strong
#'
#' @param x First categorical variable (vector or factor)
#' @param y Second categorical variable (vector or factor)
#' @param conf.level Confidence level for bootstrap CI (default 0.95 = 95% CI). Higher values (e.g., 0.99) increase interval width but confidence
#' @param correct Apply continuity correction for 2x2 tables (default TRUE). Reduces bias in small samples per Yates (1934)
#' @return List with V statistic, p-value, and confidence interval
#' @examples
#' cramers_v_result <- calculate_cramers_v(
#'   sample_metadata$Diagnosis,
#'   sample_metadata$Plasma_collection_tube_type
#' )
#' @export
calculate_cramers_v <- function(x, y, conf.level = 0.95, correct = TRUE) {
  # Remove NA values
  complete_cases <- complete.cases(x, y)
  x <- x[complete_cases]
  y <- y[complete_cases]

  if (length(x) == 0) {
    stop("No complete cases available for Cramér's V calculation")
  }

  # Calculate Cramér's V using DescTools
  v_result <- DescTools::CramerV(
    x = x,
    y = y,
    conf.level = conf.level,
    correct = correct
  )

  # If confidence interval is returned, extract components
  if (is.numeric(v_result) && length(v_result) == 3) {
    result <- list(
      statistic = v_result[1],
      conf.low = v_result[2],
      conf.high = v_result[3],
      conf.level = conf.level,
      n = length(x)
    )
  } else {
    result <- list(
      statistic = as.numeric(v_result),
      conf.low = NA,
      conf.high = NA,
      conf.level = conf.level,
      n = length(x)
    )
  }

  # Add interpretation
  result$interpretation <- dplyr::case_when(
    result$statistic < 0.1 ~ "negligible",
    result$statistic < 0.3 ~ "weak",
    result$statistic < 0.5 ~ "moderate",
    TRUE ~ "strong"
  )

  return(result)
}


#' Calculate all pairwise Cramér's V values
#'
#' @description
#' Computes Cramér's V for all pairs of categorical variables of interest:
#' Diagnosis, Tube Type, Country, Sex, Plate Number
#'
#' @param sample_metadata Data frame with sample-level metadata
#' @return Data frame with pairwise Cramér's V values
#' @examples
#' cramers_matrix <- calculate_pairwise_cramers_v(sample_metadata)
#' @export
calculate_pairwise_cramers_v <- function(sample_metadata) {
  # Variables to analyze
  vars <- c(
    "Diagnosis",
    "Plasma_collection_tube_type",
    "country",
    "Sex",
    "Olink_Plate_No"
  )

  # Initialize results
  results <- list()

  # Calculate all pairs
  for (i in 1:(length(vars) - 1)) {
    for (j in (i + 1):length(vars)) {
      var1 <- vars[i]
      var2 <- vars[j]

      tryCatch(
        {
          v_result <- calculate_cramers_v(
            sample_metadata[[var1]],
            sample_metadata[[var2]]
          )

          results[[length(results) + 1]] <- data.frame(
            var1 = var1,
            var2 = var2,
            cramers_v = v_result$statistic,
            conf.low = v_result$conf.low,
            conf.high = v_result$conf.high,
            interpretation = v_result$interpretation,
            n = v_result$n,
            stringsAsFactors = FALSE
          )
        },
        error = function(e) {
          warning(sprintf(
            "Failed to calculate Cramér's V for %s vs %s: %s",
            var1, var2, e$message
          ))
        }
      )
    }
  }

  # Combine results
  results_df <- dplyr::bind_rows(results)

  # Sort by Cramér's V (descending)
  results_df <- results_df %>%
    dplyr::arrange(dplyr::desc(cramers_v))

  message("\n=== Cramér's V Results ===")
  print(as.data.frame(results_df))

  return(results_df)
}


#' Create contingency tables for key variable pairs
#'
#' @description
#' Generates contingency tables showing sample counts for:
#' - Diagnosis × Tube Type
#' - Diagnosis × Country
#' - Tube Type × Country
#' - Diagnosis × Tube Type × Country (3-way)
#'
#' @param sample_metadata Data frame with sample-level metadata
#' @return List of contingency tables
#' @examples
#' tables <- create_contingency_tables(sample_metadata)
#' @export
create_contingency_tables <- function(sample_metadata) {
  tables <- list()

  # Diagnosis × Tube Type
  tables$diagnosis_tube <- sample_metadata %>%
    dplyr::count(Diagnosis, Plasma_collection_tube_type) %>%
    tidyr::pivot_wider(
      names_from = Plasma_collection_tube_type,
      values_from = n,
      values_fill = 0
    )

  # Diagnosis × Country
  tables$diagnosis_country <- sample_metadata %>%
    dplyr::count(Diagnosis, country) %>%
    tidyr::pivot_wider(
      names_from = country,
      values_from = n,
      values_fill = 0
    )

  # Tube Type × Country
  tables$tube_country <- sample_metadata %>%
    dplyr::count(Plasma_collection_tube_type, country) %>%
    tidyr::pivot_wider(
      names_from = country,
      values_from = n,
      values_fill = 0
    )

  # 3-way table: Diagnosis × Tube Type × Country
  tables$diagnosis_tube_country <- sample_metadata %>%
    dplyr::count(Diagnosis, Plasma_collection_tube_type, country)

  # Add row/column totals
  tables$diagnosis_tube <- tables$diagnosis_tube %>%
    dplyr::mutate(Total = rowSums(dplyr::select(., -Diagnosis)))

  message("\n=== Contingency Tables ===\n")
  message("Diagnosis × Tube Type:")
  print(tables$diagnosis_tube)

  message("\nDiagnosis × Country:")
  print(tables$diagnosis_country)

  message("\nTube Type × Country:")
  print(tables$tube_country)

  return(tables)
}


#' Perform chi-square tests for independence
#'
#' @description
#' Tests the null hypothesis that variables are independent using
#' Pearson's chi-square test. IMPORTANT: Given the perfect confounding
#' in this dataset, we EXPECT these tests to be highly significant.
#'
#' @param sample_metadata Data frame with sample-level metadata
#' @return Data frame with test results
#' @examples
#' chi_sq_results <- perform_chi_square_tests(sample_metadata)
#' @export
perform_chi_square_tests <- function(sample_metadata) {
  # Variable pairs to test
  tests <- list(
    list(var1 = "Diagnosis", var2 = "Plasma_collection_tube_type"),
    list(var1 = "Diagnosis", var2 = "country"),
    list(var1 = "Plasma_collection_tube_type", var2 = "country")
  )

  results <- list()

  for (test in tests) {
    var1 <- test$var1
    var2 <- test$var2

    # Create contingency table
    cont_table <- table(
      sample_metadata[[var1]],
      sample_metadata[[var2]]
    )

    # Perform chi-square test
    chi_test <- chisq.test(cont_table)

    # Calculate effect size (Cramér's V)
    cramers <- calculate_cramers_v(
      sample_metadata[[var1]],
      sample_metadata[[var2]]
    )

    results[[length(results) + 1]] <- data.frame(
      var1 = var1,
      var2 = var2,
      chi_square = chi_test$statistic,
      df = chi_test$parameter,
      p_value = chi_test$p.value,
      cramers_v = cramers$statistic,
      interpretation = cramers$interpretation,
      stringsAsFactors = FALSE
    )
  }

  results_df <- dplyr::bind_rows(results)

  message("\n=== Chi-Square Test Results ===")
  print(as.data.frame(results_df))

  # Flag perfect confounding
  perfect_confound <- results_df %>%
    dplyr::filter(cramers_v > 0.9)

  if (nrow(perfect_confound) > 0) {
    message("\n⚠️  WARNING: Near-perfect confounding detected!")
    message("The following variable pairs show Cramér's V > 0.9:")
    print(as.data.frame(perfect_confound %>% dplyr::select(var1, var2, cramers_v)))
  }

  return(results_df)
}


#' Analyze confounding within diagnostic groups
#'
#' @description
#' Stratifies by diagnosis and examines tube type distribution.
#' Critical for understanding the confounding pattern:
#' - Are neurological controls 100% EDTA?
#' - What % of ALS cases are HEPARIN?
#'
#' @param sample_metadata Data frame with sample-level metadata
#' @return List with stratified summaries
#' @examples
#' stratified <- analyze_stratified_confounding(sample_metadata)
#' @export
analyze_stratified_confounding <- function(sample_metadata) {
  # Overall distribution
  overall <- sample_metadata %>%
    dplyr::count(Diagnosis, Plasma_collection_tube_type) %>%
    dplyr::group_by(Diagnosis) %>%
    dplyr::mutate(
      percent = 100 * n / sum(n),
      total = sum(n)
    ) %>%
    dplyr::ungroup()

  # Identify perfect confounding
  perfect_confound_diagnosis <- overall %>%
    dplyr::group_by(Diagnosis) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup()

  # Country distribution by diagnosis
  by_country <- sample_metadata %>%
    dplyr::count(Diagnosis, country) %>%
    dplyr::group_by(Diagnosis) %>%
    dplyr::mutate(
      percent = 100 * n / sum(n),
      total = sum(n)
    ) %>%
    dplyr::ungroup()

  message("\n=== Stratified Confounding Analysis ===\n")
  message("Tube Type Distribution by Diagnosis:")
  print(as.data.frame(overall))

  message("\nCountry Distribution by Diagnosis:")
  print(as.data.frame(by_country))

  if (nrow(perfect_confound_diagnosis) > 0) {
    message("\n⚠️  CRITICAL: Perfect confounding detected!")
    message("The following diagnostic groups have only ONE tube type:")
    print(as.data.frame(perfect_confound_diagnosis %>%
      dplyr::select(Diagnosis, Plasma_collection_tube_type, n, percent)))
  }

  return(list(
    tube_by_diagnosis = overall,
    country_by_diagnosis = by_country,
    perfect_confound = perfect_confound_diagnosis
  ))
}


#' Generate comprehensive bias quantification report
#'
#' @description
#' Master function that runs all bias quantification analyses and
#' returns a comprehensive report object.
#'
#' @param sample_metadata Data frame with sample-level metadata
#' @return List with all bias analysis results
#' @examples
#' bias_report <- quantify_confounding_bias(sample_metadata)
#' @export
quantify_confounding_bias <- function(sample_metadata) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("CONFOUNDING BIAS QUANTIFICATION REPORT")
  message(paste(rep("=", 60), collapse = ""), "\n")

  # 1. Pairwise Cramér's V
  cramers_pairwise <- calculate_pairwise_cramers_v(sample_metadata)

  # 2. Contingency tables
  cont_tables <- create_contingency_tables(sample_metadata)

  # 3. Chi-square tests
  chi_sq_results <- perform_chi_square_tests(sample_metadata)

  # 4. Stratified analysis
  stratified <- analyze_stratified_confounding(sample_metadata)

  # Compile report
  report <- list(
    cramers_v_pairwise = cramers_pairwise,
    contingency_tables = cont_tables,
    chi_square_tests = chi_sq_results,
    stratified_analysis = stratified,
    sample_size = nrow(sample_metadata),
    n_diagnosis_groups = dplyr::n_distinct(sample_metadata$Diagnosis),
    n_tube_types = dplyr::n_distinct(sample_metadata$Plasma_collection_tube_type),
    n_countries = dplyr::n_distinct(sample_metadata$country)
  )

  message("\n", paste(rep("=", 60), collapse = ""))
  message("REPORT COMPLETE")
  message(paste(rep("=", 60), collapse = ""))

  return(report)
}


#' **NEW: Test for Tube Type Effects in Healthy Controls**
#'
#' @description
#' Tests whether "tube-robust" proteins differ between Italy/HEPARIN and US/EDTA
#' WITHIN HEALTHY CONTROLS ONLY.
#'
#' **Rationale**: If a protein differs by tube type/country in disease-free
#' individuals, that difference must be a technical artifact (or population
#' genetics), NOT disease biology.
#'
#' This is a powerful test because:
#'   1. No disease present → isolates tube/country effects
#'   2. Any significant difference = confounding
#'   3. Validates whether "robust" proteins are truly tube-independent
#'
#' **Interpretation**:
#'   - p < 0.05: Protein affected by tube type even without disease
#'   - Effect size: Magnitude of tube artifact in NPX units
#'   - If many robust proteins show effects → severe residual confounding
#'
#' @param protein_wide Wide-format data with proteins as columns
#' @param protein_concordance Concordance analysis results with sig_both/same_direction
#' @param alpha Significance threshold for Type I error rate (default 0.05 = 5% false positive rate). Adjust for multiple testing considerations
#' @return Data frame with test results for each tube-robust protein
#' @examples
#' healthy_tube_effects <- tube_effect_in_healthy(protein_wide, protein_concordance)
#' @export
tube_effect_in_healthy <- function(protein_wide, protein_concordance, alpha = 0.05) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("TESTING TUBE EFFECTS IN HEALTHY CONTROLS")
  message("(Disease-free individuals → isolates technical artifacts)")
  message(paste(rep("=", 70), collapse = ""))

  # Get the tube-robust proteins
  robust_proteins <- protein_concordance %>%
    dplyr::filter(sig_both, same_direction) %>%
    dplyr::pull(Assay)

  n_robust <- length(robust_proteins)
  message(sprintf("\nTesting %d tube-robust proteins in healthy controls only\n", n_robust))

  if (n_robust == 0) {
    stop("No tube-robust proteins found in protein_concordance data")
  }

  # Filter to healthy controls only (no disease)
  healthy_only <- protein_wide %>%
    dplyr::filter(Diagnosis == "Healthy_control", country != "Unknown")

  n_healthy <- nrow(healthy_only)
  message(sprintf("Sample size: %d healthy controls", n_healthy))
  message(sprintf("  - Italy/HEPARIN: %d", sum(healthy_only$country == "Italy")))
  message(sprintf("  - US/EDTA: %d\n", sum(healthy_only$country == "US")))

  if (n_healthy < 20) {
    warning("Very small sample size for healthy controls. Results may be unstable.")
  }

  # Test each protein: Italy/HEPARIN vs US/EDTA
  # Model: protein ~ country + Age_Collection + Sex
  results <- purrr::map_df(robust_proteins, function(prot) {
    # Check if protein column exists
    if (!prot %in% names(healthy_only)) {
      warning(sprintf("Protein %s not found in data, skipping", prot))
      return(NULL)
    }

    # Build model
    formula_str <- sprintf("`%s` ~ country + Age_Collection + Sex", prot)

    tryCatch(
      {
        model <- lm(as.formula(formula_str), data = healthy_only)

        # Extract country coefficient (Italy is reference, so countryUS is the effect)
        model_summary <- broom::tidy(model)
        country_effect <- model_summary %>%
          dplyr::filter(term == "countryUS")

        if (nrow(country_effect) == 0) {
          return(NULL)
        }

        # Calculate effect size in NPX units
        mean_italy <- mean(healthy_only[[prot]][healthy_only$country == "Italy"], na.rm = TRUE)
        mean_us <- mean(healthy_only[[prot]][healthy_only$country == "US"], na.rm = TRUE)

        data.frame(
          protein = prot,
          coefficient = country_effect$estimate,
          std_error = country_effect$std.error,
          t_statistic = country_effect$statistic,
          p_value = country_effect$p.value,
          mean_NPX_italy = mean_italy,
          mean_NPX_us = mean_us,
          raw_difference = mean_us - mean_italy,
          significant = country_effect$p.value < alpha,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        warning(sprintf("Failed to fit model for %s: %s", prot, e$message))
        return(NULL)
      }
    )
  })

  # Remove NULL results
  results <- results %>%
    dplyr::filter(!is.na(protein))

  # Sort by p-value
  results <- results %>%
    dplyr::arrange(p_value)

  # Summary statistics
  n_significant <- sum(results$significant)
  pct_significant <- 100 * n_significant / nrow(results)

  message("\n=== RESULTS ===")
  message(sprintf("Proteins tested: %d", nrow(results)))
  message(sprintf("Significant at α=%.3f: %d (%.1f%%)", alpha, n_significant, pct_significant))

  if (n_significant > 0) {
    message("\n⚠️  TUBE EFFECTS DETECTED IN HEALTHY CONTROLS!")
    message(sprintf("  → %d/%d tube-robust proteins differ by tube type", n_significant, nrow(results)))
    message("  → These differences occur WITHOUT disease present")
    message("  → Evidence of residual technical confounding")

    # Show top 5 most significant
    message("\nTop 5 proteins with strongest tube effects:")
    top5 <- results %>%
      dplyr::filter(significant) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::select(protein, raw_difference, p_value)
    print(as.data.frame(top5))
  } else {
    message("\n✓ No significant tube effects detected in healthy controls")
    message("  → Tube-robust proteins appear truly tube-independent")
  }

  # Add FDR correction
  results <- results %>%
    dplyr::mutate(
      p_value_fdr = p.adjust(p_value, method = "fdr"),
      significant_fdr = p_value_fdr < alpha
    )

  n_significant_fdr <- sum(results$significant_fdr)
  message(sprintf(
    "\nAfter FDR correction: %d significant (%.1f%%)",
    n_significant_fdr, 100 * n_significant_fdr / nrow(results)
  ))

  message(paste(rep("=", 70), collapse = ""), "\n")

  return(results)
}


#' **NEW: Decompose Effect Sizes - Disease vs Tube Type**
#'
#' @description
#' For each tube-robust protein, compares the magnitude of:
#'   1. Disease effect (ALS vs healthy controls)
#'   2. Tube type effect (Italy/HEPARIN vs US/EDTA)
#'
#' **Key question**: Is the tube type effect comparable to or larger than the
#' disease effect? If tube effects are 30-60% as large as disease effects,
#' it's impossible to disentangle them with perfect confounding.
#'
#' **Calculation strategy**:
#'   - Disease effect: ALS vs Healthy, adjusted for country + tube type
#'   - Tube effect: Italy vs US, adjusted for diagnosis
#'   - Ratio: |tube effect| / |disease effect|
#'
#' **Interpretation**:
#'   - Ratio < 0.2: Tube effect negligible compared to disease (good)
#'   - Ratio 0.2-0.5: Tube effect moderate (caution)
#'   - Ratio > 0.5: Tube effect comparable to disease (severe entanglement)
#'   - Ratio > 1.0: Tube effect LARGER than disease effect (artifact?)
#'
#' @param protein_wide Wide-format data with proteins as columns
#' @param protein_concordance Concordance analysis results with sig_both/same_direction
#' @param alpha Significance threshold (default 0.05)
#' @return Data frame with effect sizes and ratios for each tube-robust protein
#' @examples
#' effect_decomp <- effect_size_decomposition(protein_wide, protein_concordance)
#' @export
effect_size_decomposition <- function(protein_wide, protein_concordance, alpha = 0.05) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("EFFECT SIZE DECOMPOSITION: DISEASE vs TUBE TYPE")
  message("(Comparing magnitude of biological vs technical effects)")
  message(paste(rep("=", 70), collapse = ""))

  # Get the tube-robust proteins
  robust_proteins <- protein_concordance %>%
    dplyr::filter(sig_both, same_direction) %>%
    dplyr::pull(Assay)

  n_robust <- length(robust_proteins)
  message(sprintf("\nDecomposing effects for %d tube-robust proteins\n", n_robust))

  if (n_robust == 0) {
    stop("No tube-robust proteins found in protein_concordance data")
  }

  # Prepare data - exclude Unknown country
  data_clean <- protein_wide %>%
    dplyr::filter(country != "Unknown")

  # For disease effect, restrict to ALS and Healthy controls only (exclude Neuro OND)
  # Reason: Neuro OND are 100% EDTA, creating perfect confounding
  data_disease <- data_clean %>%
    dplyr::filter(Diagnosis %in% c("ALS", "Healthy_control"))

  n_samples_disease <- nrow(data_disease)
  n_samples_tube <- nrow(data_clean)

  message(sprintf("Samples for disease effect: %d (ALS + Healthy only)", n_samples_disease))
  message(sprintf("Samples for tube effect: %d (all diagnoses)\n", n_samples_tube))

  # Decompose effects for each protein
  results <- purrr::map_df(robust_proteins, function(prot) {
    # Check if protein exists
    if (!prot %in% names(data_clean)) {
      warning(sprintf("Protein %s not found in data, skipping", prot))
      return(NULL)
    }

    tryCatch(
      {
        # ===== Model 1: Disease effect (ALS vs Healthy, adjusted for country) =====
        # Diagnosis is the predictor, adjusting for country/tube confounding
        formula_disease <- sprintf("`%s` ~ Diagnosis + country + Age_Collection + Sex", prot)
        model_disease <- lm(as.formula(formula_disease), data = data_disease)
        coef_disease <- broom::tidy(model_disease) %>%
          dplyr::filter(stringr::str_detect(term, "Diagnosis"))

        # ALS vs Healthy effect (reference is usually alphabetical, but check)
        disease_effect <- coef_disease$estimate[1]
        disease_pvalue <- coef_disease$p.value[1]
        disease_stderr <- coef_disease$std.error[1]


        # ===== Model 2: Tube effect (Italy vs US, adjusted for diagnosis) =====
        # Country is the predictor, adjusting for disease status
        formula_tube <- sprintf("`%s` ~ country + Diagnosis + Age_Collection + Sex", prot)
        model_tube <- lm(as.formula(formula_tube), data = data_clean)
        coef_tube <- broom::tidy(model_tube) %>%
          dplyr::filter(term == "countryUS")

        if (nrow(coef_tube) == 0) {
          return(NULL)
        }

        tube_effect <- coef_tube$estimate
        tube_pvalue <- coef_tube$p.value
        tube_stderr <- coef_tube$std.error


        # ===== Calculate ratio and classification =====
        abs_disease_effect <- abs(disease_effect)
        abs_tube_effect <- abs(tube_effect)

        # Avoid division by zero
        if (abs_disease_effect < 1e-10) {
          ratio <- NA
          interpretation <- "disease_effect_zero"
        } else {
          ratio <- abs_tube_effect / abs_disease_effect

          interpretation <- dplyr::case_when(
            ratio < 0.2 ~ "negligible_tube_effect",
            ratio < 0.5 ~ "moderate_tube_effect",
            ratio < 1.0 ~ "severe_tube_effect",
            TRUE ~ "tube_dominates_disease"
          )
        }

        # Return results
        data.frame(
          protein = prot,
          disease_effect = disease_effect,
          disease_pvalue = disease_pvalue,
          disease_stderr = disease_stderr,
          disease_significant = disease_pvalue < alpha,
          tube_effect = tube_effect,
          tube_pvalue = tube_pvalue,
          tube_stderr = tube_stderr,
          tube_significant = tube_pvalue < alpha,
          abs_disease_effect = abs_disease_effect,
          abs_tube_effect = abs_tube_effect,
          ratio_tube_to_disease = ratio,
          interpretation = interpretation,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        warning(sprintf("Failed to decompose effects for %s: %s", prot, e$message))
        return(NULL)
      }
    )
  })

  # Remove NULL results
  results <- results %>%
    dplyr::filter(!is.na(protein))

  # Sort by ratio (descending)
  results <- results %>%
    dplyr::arrange(dplyr::desc(ratio_tube_to_disease))

  # Summary statistics
  message("\n=== EFFECT SIZE DECOMPOSITION RESULTS ===\n")
  message(sprintf("Proteins analyzed: %d", nrow(results)))

  # Count by interpretation category
  interpretation_counts <- results %>%
    dplyr::count(interpretation) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nDistribution of tube/disease effect ratios:")
  print(as.data.frame(interpretation_counts))

  # Identify problematic proteins (ratio > 0.5)
  severe_entanglement <- results %>%
    dplyr::filter(ratio_tube_to_disease >= 0.5) %>%
    nrow()

  pct_severe <- 100 * severe_entanglement / nrow(results)

  if (severe_entanglement > 0) {
    message(sprintf(
      "\n⚠️  SEVERE ENTANGLEMENT: %d/%d proteins (%.1f%%)",
      severe_entanglement, nrow(results), pct_severe
    ))
    message("  → Tube effects are ≥50% as large as disease effects")
    message("  → Cannot reliably separate biological from technical signal")

    # Show top 5 worst cases
    message("\nTop 5 proteins with highest tube/disease ratio:")
    top5 <- results %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::select(protein, abs_disease_effect, abs_tube_effect, ratio_tube_to_disease)
    print(as.data.frame(top5))
  } else {
    message("\n✓ No severe entanglement detected")
    message("  → Tube effects are <50% of disease effects for all proteins")
  }

  # Calculate median ratio
  median_ratio <- median(results$ratio_tube_to_disease, na.rm = TRUE)
  message(sprintf("\nMedian tube/disease ratio: %.3f", median_ratio))

  if (median_ratio > 0.5) {
    message("  → CRITICAL: Tube effects typically comparable to disease effects")
  } else if (median_ratio > 0.2) {
    message("  → CAUTION: Tube effects are moderate relative to disease effects")
  } else {
    message("  → Tube effects are relatively small compared to disease effects")
  }

  # Additional context: proteins where tube effect is significant but disease is not
  tube_only_sig <- results %>%
    dplyr::filter(tube_significant, !disease_significant) %>%
    nrow()

  if (tube_only_sig > 0) {
    message(sprintf(
      "\n⚠️  %d proteins: Tube effect significant BUT disease effect NOT significant",
      tube_only_sig
    ))
    message("  → Suggests these may be pure technical artifacts")
  }

  message(paste(rep("=", 70), collapse = ""), "\n")

  return(results)
}
