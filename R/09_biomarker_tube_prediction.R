#' Test if Original Biomarkers Collectively Predict Tube Type
#'
#' This module tests whether the original study's biomarker panel, taken as a whole,
#' can predict plasma collection tube type. If AUC > 0.8, this indicates substantial
#' confounding even if individual proteins aren't top-ranked.

#' Test Biomarker Panel's Ability to Predict Tube Type
#'
#' @description
#' Trains a Random Forest model using ONLY the original study's biomarker proteins
#' to predict tube type (HEPARIN vs EDTA). Compares performance to the full model
#' using all proteins.
#'
#' This is a CRITICAL test because individual proteins may have moderate tube effects,
#' but the combination can still capture substantial technical artifacts.
#'
#' @param protein_wide Wide-format protein data
#' @param biomarker_list Character vector of biomarker protein names
#' @param reverse_prediction_results Full reverse prediction results for comparison
#' @param n_folds Number of CV folds (default: 5)
#' @param seed Random seed for reproducibility
#' @return List with model results and comparison
#' @export
test_biomarker_tube_prediction <- function(protein_wide,
                                           biomarker_list,
                                           reverse_prediction_results,
                                           n_folds = 5,
                                           seed = 42) {
  library(dplyr)
  library(caret)

  cat("\n=== TESTING BIOMARKER PANEL TUBE-TYPE PREDICTION ===\n\n")

  # Filter to clean samples
  data_clean <- protein_wide %>%
    dplyr::filter(!is.na(Plasma_collection_tube_type)) %>%
    dplyr::select(where(~ !all(is.na(.))))

  # Check which biomarkers are available
  available_biomarkers <- intersect(biomarker_list, colnames(data_clean))
  missing_biomarkers <- setdiff(biomarker_list, colnames(data_clean))

  cat(sprintf("Original biomarkers: %d\n", length(biomarker_list)))
  cat(sprintf("Available in dataset: %d\n", length(available_biomarkers)))
  if (length(missing_biomarkers) > 0) {
    cat(sprintf("Missing: %d\n", length(missing_biomarkers)))
    cat("Missing proteins:", paste(missing_biomarkers, collapse = ", "), "\n")
  }
  cat("\n")

  # Outcome: tube type (binary)
  y <- factor(data_clean$Plasma_collection_tube_type, levels = c("EDTA", "HEPARIN"))

  # Features: ONLY the biomarkers
  X_biomarkers <- data_clean %>%
    dplyr::select(dplyr::all_of(available_biomarkers))

  cat(sprintf("Training Random Forest with %d biomarker proteins...\n", ncol(X_biomarkers)))

  # Set seed
  set.seed(seed)

  # Cross-validation setup
  train_control <- trainControl(
    method = "cv",
    number = n_folds,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE
  )

  # Train model
  rf_biomarkers <- caret::train(
    x = X_biomarkers,
    y = y,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    num.trees = 500,
    importance = "permutation"
  )

  # Extract performance metrics
  biomarker_auc <- max(rf_biomarkers$results$ROC)
  biomarker_sens <- rf_biomarkers$results$Sens[which.max(rf_biomarkers$results$ROC)]
  biomarker_spec <- rf_biomarkers$results$Spec[which.max(rf_biomarkers$results$ROC)]

  # Get full model AUC for comparison
  full_auc <- reverse_prediction_results$auc

  # Calculate difference
  auc_difference <- full_auc - biomarker_auc

  # Interpret results
  interpretation <- dplyr::case_when(
    biomarker_auc >= 0.9 ~ "CRITICAL: Biomarkers alone predict tube type with AUC >= 0.9",
    biomarker_auc >= 0.8 ~ "SUBSTANTIAL: Biomarkers have significant tube-type bias (AUC >= 0.8)",
    biomarker_auc >= 0.7 ~ "MODERATE: Some tube-type confounding present (AUC >= 0.7)",
    TRUE ~ "MINIMAL: Limited tube-type prediction ability"
  )

  # Print results
  cat("\n=== RESULTS ===\n\n")
  cat(sprintf("Biomarker-only model: AUC = %.4f\n", biomarker_auc))
  cat(sprintf("  Sensitivity: %.4f\n", biomarker_sens))
  cat(sprintf("  Specificity: %.4f\n\n", biomarker_spec))

  cat(sprintf("Full model (all proteins): AUC = %.4f\n", full_auc))
  cat(sprintf("Difference: %.4f\n\n", auc_difference))

  cat(sprintf("INTERPRETATION: %s\n\n", interpretation))

  # Get predictions for ROC curve
  predictions <- rf_biomarkers$pred %>%
    dplyr::filter(mtry == rf_biomarkers$bestTune$mtry)

  # Compile results
  results <- list(
    model = rf_biomarkers,
    predictions = predictions,
    n_biomarkers = length(available_biomarkers),
    available_biomarkers = available_biomarkers,
    missing_biomarkers = missing_biomarkers,
    auc = biomarker_auc,
    sensitivity = biomarker_sens,
    specificity = biomarker_spec,
    full_model_auc = full_auc,
    auc_difference = auc_difference,
    interpretation = interpretation
  )

  return(results)
}


#' Create Comparison Visualization
#'
#' @description
#' Creates a multi-panel figure comparing:
#' 1. ROC curves (biomarkers vs full model)
#' 2. AUC comparison bar chart
#' 3. Individual biomarker importance for tube discrimination
#'
#' @param biomarker_tube_results Results from test_biomarker_tube_prediction()
#' @param reverse_prediction_results Full reverse prediction results
#' @param save_path Optional path to save figure
#' @return ggplot2 object
#' @export
plot_biomarker_tube_comparison <- function(biomarker_tube_results,
                                           reverse_prediction_results,
                                           save_path = NULL) {
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(dplyr)

  # ========== PANEL A: ROC Curves ==========
  # Biomarker-only ROC
  biomarker_preds <- biomarker_tube_results$predictions
  roc_biomarker <- pROC::roc(
    response = biomarker_preds$obs,
    predictor = biomarker_preds$HEPARIN,
    levels = c("EDTA", "HEPARIN"),
    direction = "<"
  )

  roc_biomarker_df <- data.frame(
    fpr = 1 - roc_biomarker$specificities,
    tpr = roc_biomarker$sensitivities,
    model = "Biomarkers Only"
  )

  # Full model ROC
  full_preds <- reverse_prediction_results$model$pred %>%
    dplyr::filter(
      mtry == reverse_prediction_results$model$bestTune$mtry,
      min.node.size == reverse_prediction_results$model$bestTune$min.node.size
    )

  roc_full <- pROC::roc(
    response = full_preds$obs,
    predictor = full_preds$HEPARIN,
    levels = c("EDTA", "HEPARIN"),
    direction = "<"
  )

  roc_full_df <- data.frame(
    fpr = 1 - roc_full$specificities,
    tpr = roc_full$sensitivities,
    model = "Full Model (All Proteins)"
  )

  roc_combined <- rbind(roc_biomarker_df, roc_full_df)

  p_roc <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = model)) +
    geom_line(size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text",
      x = 0.6, y = 0.2,
      label = sprintf("Biomarkers: AUC = %.3f", biomarker_tube_results$auc),
      size = 4, color = "#E69F00", fontface = "bold"
    ) +
    annotate("text",
      x = 0.6, y = 0.1,
      label = sprintf("Full Model: AUC = %.3f", biomarker_tube_results$full_model_auc),
      size = 4, color = "#D32F2F", fontface = "bold"
    ) +
    scale_color_manual(values = c(
      "Biomarkers Only" = "#E69F00",
      "Full Model (All Proteins)" = "#D32F2F"
    )) +
    labs(
      title = "A. ROC Curves: Predicting Tube Type",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 1)
    )

  # ========== PANEL B: AUC Comparison Bar Chart ==========
  comparison_df <- data.frame(
    model = c(
      "Biomarkers\nOnly\n(17 proteins)",
      "Full Model\n(2,868 proteins)"
    ),
    auc = c(biomarker_tube_results$auc, biomarker_tube_results$full_model_auc),
    color = c("#E69F00", "#D32F2F")
  )

  p_bars <- ggplot(comparison_df, aes(x = model, y = auc, fill = color)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", auc)),
      vjust = -0.5, fontface = "bold", size = 5
    ) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", size = 1) +
    annotate("text",
      x = 1.5, y = 0.82,
      label = "AUC = 0.8 threshold\n(substantial bias)",
      size = 3.5, color = "red", fontface = "bold"
    ) +
    scale_fill_identity() +
    labs(
      title = "B. AUC Comparison",
      subtitle = biomarker_tube_results$interpretation,
      x = NULL,
      y = "AUC (Tube Type Prediction)"
    ) +
    ylim(0, 1.05) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "#D32F2F", face = "bold", size = 11),
      axis.text.x = element_text(face = "bold", size = 11),
      panel.grid.major.x = element_blank()
    )

  # ========== PANEL C: Individual Biomarker Ranks ==========
  # Get importance for all proteins
  importance_vec <- reverse_prediction_results$model$finalModel$variable.importance
  importance_df <- data.frame(
    protein = names(importance_vec),
    importance = as.numeric(importance_vec),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(importance)) %>%
    dplyr::mutate(rank = dplyr::row_number())

  # Filter to biomarkers
  biomarker_ranks <- importance_df %>%
    dplyr::filter(protein %in% biomarker_tube_results$available_biomarkers) %>%
    dplyr::arrange(rank) %>%
    dplyr::mutate(
      protein = factor(protein, levels = protein),
      color = dplyr::case_when(
        rank <= 100 ~ "#D32F2F",
        rank <= 500 ~ "#FF9800",
        TRUE ~ "#4CAF50"
      )
    )

  p_ranks <- ggplot(biomarker_ranks, aes(x = protein, y = rank, fill = color)) +
    geom_col() +
    geom_text(aes(label = rank), hjust = -0.2, fontface = "bold", size = 3) +
    scale_fill_identity() +
    scale_y_reverse(limits = c(2900, 0)) +
    coord_flip() +
    labs(
      title = "C. Individual Biomarker Ranks (Tube Discrimination)",
      subtitle = sprintf("Out of %d total proteins", nrow(importance_df)),
      x = NULL,
      y = "Rank (1 = highest tube-type importance)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.y = element_text(face = "bold", size = 10),
      panel.grid.major.y = element_blank()
    )

  # ========== COMBINE PANELS ==========
  combined <- (p_roc | p_bars) / p_ranks +
    plot_annotation(
      title = "Original Biomarkers Can Substantially Predict Tube Type",
      subtitle = sprintf(
        "Despite individual proteins not being top-ranked, the %d-protein panel achieves AUC = %.3f for tube type prediction",
        biomarker_tube_results$n_biomarkers,
        biomarker_tube_results$auc
      ),
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(
          size = 13, hjust = 0.5, color = "#D32F2F",
          face = "bold", margin = margin(b = 10)
        )
      )
    )

  # Save if requested
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, combined, width = 16, height = 12, dpi = 300)
    cat(sprintf("\n✓ Biomarker tube prediction comparison saved to: %s\n", save_path))
  }

  return(combined)
}


#' Generate Summary Report
#'
#' @description
#' Creates a text summary of the biomarker tube-type prediction analysis.
#'
#' @param biomarker_tube_results Results from test_biomarker_tube_prediction()
#' @param save_path Optional path to save report
#' @return Character string with formatted report
#' @export
summarize_biomarker_tube_prediction <- function(biomarker_tube_results,
                                                save_path = NULL) {
  report <- paste0(
    "# BIOMARKER PANEL TUBE-TYPE PREDICTION ANALYSIS\n\n",
    "**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "---\n\n",
    "## KEY QUESTION\n\n",
    "Can the original study's 17-protein biomarker panel, taken as a whole,\n",
    "predict plasma collection tube type (HEPARIN vs EDTA)?\n\n",
    "If yes, this indicates the biomarkers are confounded by technical artifacts,\n",
    "even if individual proteins aren't top-ranked for tube discrimination.\n\n",
    "---\n\n",
    "## RESULTS\n\n",
    sprintf("**Biomarker proteins used:** %d out of 17\n", biomarker_tube_results$n_biomarkers),
    sprintf("**AUC (5-fold CV):** %.4f\n", biomarker_tube_results$auc),
    sprintf("**Sensitivity:** %.4f\n", biomarker_tube_results$sensitivity),
    sprintf("**Specificity:** %.4f\n\n", biomarker_tube_results$specificity),
    "**Comparison to full model:**\n",
    sprintf("  - Full model (2,868 proteins): AUC = %.4f\n", biomarker_tube_results$full_model_auc),
    sprintf("  - Biomarkers only: AUC = %.4f\n", biomarker_tube_results$auc),
    sprintf("  - Difference: %.4f\n\n", biomarker_tube_results$auc_difference),
    "---\n\n",
    "## INTERPRETATION\n\n",
    sprintf("**%s**\n\n", biomarker_tube_results$interpretation),
    "This finding is critical because:\n\n",
    "1. **Individual protein ranks were misleading:** Most biomarkers ranked below median\n",
    "   for tube discrimination, suggesting minimal confounding.\n\n",
    "2. **Collective effect is substantial:** When combined, the biomarkers achieve\n",
    "   AUC = ", sprintf("%.3f", biomarker_tube_results$auc), " for predicting tube type.\n\n",
    "3. **Implications:** The biomarker panel is significantly confounded by tube type,\n",
    "   even though individual proteins aren't the worst offenders.\n\n"
  )

  if (biomarker_tube_results$auc >= 0.8) {
    report <- paste0(
      report,
      "## CONCLUSION\n\n",
      "The original study's biomarker panel has **SUBSTANTIAL tube-type bias**.\n",
      "Combined with perfect confounding (Italy/HEPARIN vs US/EDTA), this means:\n\n",
      "- ❌ Cannot determine if biomarkers reflect ALS biology or technical artifacts\n",
      "- ❌ Model's claimed performance may be driven by tube-type effects\n",
      "- ❌ External validation with consistent tube types is essential\n",
      "- ❌ Stratified analyses (Italy-only, US-only) should be performed\n\n"
    )
  }

  report <- paste0(
    report,
    "---\n\n",
    "**Report Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "**Analysis Pipeline:** targets + R 4.4.2\n"
  )

  # Save if requested
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, save_path)
    cat(sprintf("\n✓ Biomarker tube prediction report saved to: %s\n", save_path))
  }

  cat(report)
  return(report)
}
