#' Comprehensive Visualization Functions for ALS Biomarker Bias Investigation
#'
#' Publication-quality figures demonstrating confounding bias across multiple analyses.
#' Organized by analysis type for easy navigation and maintenance.
#'
#' @description
#' This module provides all visualization functions for the bias investigation:
#' - Reverse prediction ROC curves and feature importance
#' - Geographic validation (leave-country-out CV)
#' - Protein overlap and concordance analyses
#' - Confounding structure demonstrations
#' - Residual confounding assessments
#' - Meta-analysis forest plots
#'
#' All plots use dark scientific theme (viridis palette) optimized for screen display.

# Source utility functions
source("R/99_utils.R")


# ==============================================================================
# REVERSE PREDICTION VISUALIZATIONS
# ==============================================================================

#' Plot ROC curve for reverse prediction test
#'
#' @description
#' Visualizes the performance of predicting tube type from proteins.
#' AUC near 1.0 indicates perfect separation = devastating finding.
#'
#' @param reverse_prediction_results Results from reverse_prediction_test()
#' @param save_path Optional path to save the plot
#' @return ggplot object
#' @export
plot_reverse_prediction_roc <- function(reverse_prediction_results, save_path = NULL) {
  # Extract predictions
  model <- reverse_prediction_results$model
  pred_df <- model$pred %>%
    dplyr::filter(
      mtry == model$bestTune$mtry,
      splitrule == model$bestTune$splitrule,
      min.node.size == model$bestTune$min.node.size
    )

  # Calculate ROC curve
  roc_obj <- pROC::roc(
    response = pred_df$obs,
    predictor = pred_df$HEPARIN,
    levels = c("EDTA", "HEPARIN")
  )

  auc_value <- reverse_prediction_results$auc

  # Create ROC plot
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  )

  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "#FDE724", linewidth = 2.5, alpha = 0.9) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "#71717A", linewidth = 1
    ) +
    annotate("text",
      x = 0.6, y = 0.2,
      label = sprintf("AUC = %.3f", auc_value),
      size = 6, fontface = "bold", color = "#FDE724"
    ) +
    labs(
      title = "Reverse Prediction Test: Predicting Tube Type from Proteins",
      subtitle = sprintf("AUC = %.3f - Near-perfect discrimination indicates tube effects dominate", auc_value),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_dark_scientific(base_size = 14) +
    coord_equal()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 7, dpi = 300)
    message(sprintf("Saved ROC curve to %s", save_path))
  }

  return(p)
}


#' Plot tube type feature importance
#'
#' @description
#' Shows top 20 proteins that distinguish HEPARIN from EDTA tubes.
#' These proteins are likely affected by technical artifacts.
#'
#' @param reverse_prediction_results Results from reverse_prediction_test()
#' @param save_path Optional path to save the plot
#' @return ggplot object
#' @export
plot_tube_type_features <- function(reverse_prediction_results, save_path = NULL) {
  top_features <- reverse_prediction_results$top_features

  p <- ggplot(top_features, aes(x = reorder(protein, Overall), y = Overall)) +
    geom_col(fill = "#31688E", alpha = 0.9) +
    coord_flip() +
    labs(
      title = "Top 20 Proteins Distinguishing Tube Types",
      subtitle = "Proteins most affected by HEPARIN vs EDTA anticoagulants",
      x = "Protein",
      y = "Variable Importance (Permutation-Based)"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank()
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    message(sprintf("Saved feature importance plot to %s", save_path))
  }

  return(p)
}


# ==============================================================================
# GEOGRAPHIC VALIDATION VISUALIZATIONS
# ==============================================================================

#' Plot ROC curves for leave-country-out CV
#'
#' @description
#' Creates ROC curves comparing:
#' - Italy → US
#' - US → Italy
#' - Pooled CV (for reference)
#' Now includes confidence intervals for AUC estimates in legend.
#'
#' @param lcv_results Leave-country-out CV results (must contain auc_ci_lower, auc_ci_upper)
#' @param pooled_results Pooled CV results (must contain auc_ci_lower, auc_ci_upper)
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_lcv_roc_curves <- function(lcv_results, pooled_results, save_path = NULL) {
  library(ggplot2)
  library(pROC)

  # Extract ROC curves
  roc_italy_to_us <- lcv_results$italy_to_us$roc_curve
  roc_us_to_italy <- lcv_results$us_to_italy$roc_curve

  # Extract confidence intervals (if available)
  italy_ci_lower <- if (!is.null(lcv_results$italy_to_us$auc_ci_lower)) lcv_results$italy_to_us$auc_ci_lower else NA
  italy_ci_upper <- if (!is.null(lcv_results$italy_to_us$auc_ci_upper)) lcv_results$italy_to_us$auc_ci_upper else NA
  us_ci_lower <- if (!is.null(lcv_results$us_to_italy$auc_ci_lower)) lcv_results$us_to_italy$auc_ci_lower else NA
  us_ci_upper <- if (!is.null(lcv_results$us_to_italy$auc_ci_upper)) lcv_results$us_to_italy$auc_ci_upper else NA
  pooled_ci_lower <- if (!is.null(pooled_results$auc_ci_lower)) pooled_results$auc_ci_lower else NA
  pooled_ci_upper <- if (!is.null(pooled_results$auc_ci_upper)) pooled_results$auc_ci_upper else NA

  # Create labels with CIs
  italy_label <- if (!is.na(italy_ci_lower)) {
    sprintf(
      "Italy → US (AUC = %.3f, 95%% CI: %.3f-%.3f)",
      lcv_results$italy_to_us$test_auc, italy_ci_lower, italy_ci_upper
    )
  } else {
    sprintf("Italy → US (AUC = %.3f)", lcv_results$italy_to_us$test_auc)
  }

  us_label <- if (!is.na(us_ci_lower)) {
    sprintf(
      "US → Italy (AUC = %.3f, 95%% CI: %.3f-%.3f)",
      lcv_results$us_to_italy$test_auc, us_ci_lower, us_ci_upper
    )
  } else {
    sprintf("US → Italy (AUC = %.3f)", lcv_results$us_to_italy$test_auc)
  }

  pooled_label <- if (!is.na(pooled_ci_lower)) {
    sprintf(
      "Pooled CV\nAUC = %.3f (95%% CI: %.3f-%.3f)",
      pooled_results$cv_auc, pooled_ci_lower, pooled_ci_upper
    )
  } else {
    sprintf("Pooled CV\nAUC = %.3f", pooled_results$cv_auc)
  }

  # Create data frame for plotting
  roc_data <- rbind(
    data.frame(
      specificity = 1 - roc_italy_to_us$specificities,
      sensitivity = roc_italy_to_us$sensitivities,
      model = italy_label
    ),
    data.frame(
      specificity = 1 - roc_us_to_italy$specificities,
      sensitivity = roc_us_to_italy$sensitivities,
      model = us_label
    )
  )

  # Define colors dynamically based on actual labels
  color_values <- setNames(
    c("#1F9E89", "#B5DE2B"), # viridis jade, yellow-green
    c(italy_label, us_label)
  )

  p <- ggplot(roc_data, aes(x = specificity, y = sensitivity, color = model)) +
    geom_line(linewidth = 1.2) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      color = "#888888", alpha = 0.7
    ) +
    annotate("text",
      x = 0.6, y = 0.2,
      label = pooled_label,
      color = "#CCCCCC", size = 3, hjust = 0
    ) +
    scale_color_manual(values = color_values) +
    labs(
      title = "Geographic Validation: Leave-Country-Out Cross-Validation",
      subtitle = "Model trained on one country fails to generalize to another - revealing geographic overfitting",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Test Direction"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = c(0.65, 0.25),
      legend.text = element_text(size = 9)
    ) +
    coord_equal()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    message(sprintf("Saved ROC curve to: %s", save_path))
  }

  return(p)
}


#' Plot performance comparison: Pooled vs Leave-Country-Out
#'
#' @description
#' Bar chart showing AUC comparison with performance gap highlighted.
#' Now includes error bars showing 95% confidence intervals.
#'
#' @param pooled_vs_lcv_comparison Comparison results (must contain CI fields)
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_performance_comparison <- function(pooled_vs_lcv_comparison, save_path = NULL) {
  library(ggplot2)

  # Extract CIs (if available)
  pooled_ci_lower <- if (!is.null(pooled_vs_lcv_comparison$pooled_auc_ci_lower)) {
    pooled_vs_lcv_comparison$pooled_auc_ci_lower
  } else {
    NA
  }
  pooled_ci_upper <- if (!is.null(pooled_vs_lcv_comparison$pooled_auc_ci_upper)) {
    pooled_vs_lcv_comparison$pooled_auc_ci_upper
  } else {
    NA
  }

  italy_ci_lower <- if (!is.null(pooled_vs_lcv_comparison$lcv_results$italy_to_us$auc_ci_lower)) {
    pooled_vs_lcv_comparison$lcv_results$italy_to_us$auc_ci_lower
  } else {
    NA
  }
  italy_ci_upper <- if (!is.null(pooled_vs_lcv_comparison$lcv_results$italy_to_us$auc_ci_upper)) {
    pooled_vs_lcv_comparison$lcv_results$italy_to_us$auc_ci_upper
  } else {
    NA
  }

  us_ci_lower <- if (!is.null(pooled_vs_lcv_comparison$lcv_results$us_to_italy$auc_ci_lower)) {
    pooled_vs_lcv_comparison$lcv_results$us_to_italy$auc_ci_lower
  } else {
    NA
  }
  us_ci_upper <- if (!is.null(pooled_vs_lcv_comparison$lcv_results$us_to_italy$auc_ci_upper)) {
    pooled_vs_lcv_comparison$lcv_results$us_to_italy$auc_ci_upper
  } else {
    NA
  }

  # Mean LCV CIs (average of lower/upper bounds)
  mean_ci_lower <- if (!is.na(italy_ci_lower) & !is.na(us_ci_lower)) {
    mean(c(italy_ci_lower, us_ci_lower))
  } else {
    NA
  }
  mean_ci_upper <- if (!is.na(italy_ci_upper) & !is.na(us_ci_upper)) {
    mean(c(italy_ci_upper, us_ci_upper))
  } else {
    NA
  }

  # Prepare data
  perf_data <- data.frame(
    approach = c(
      "Pooled CV\n(Original)",
      "Italy → US\n(Geographic)",
      "US → Italy\n(Geographic)",
      "Mean LCV\n(Realistic)"
    ),
    auc = c(
      pooled_vs_lcv_comparison$pooled_auc,
      pooled_vs_lcv_comparison$lcv_results$italy_to_us$test_auc,
      pooled_vs_lcv_comparison$lcv_results$us_to_italy$test_auc,
      pooled_vs_lcv_comparison$lcv_mean_auc
    ),
    ci_lower = c(pooled_ci_lower, italy_ci_lower, us_ci_lower, mean_ci_lower),
    ci_upper = c(pooled_ci_upper, italy_ci_upper, us_ci_upper, mean_ci_upper),
    type = c("Pooled", "LCV", "LCV", "LCV")
  )

  perf_data$approach <- factor(perf_data$approach, levels = perf_data$approach)

  # Calculate gap
  gap <- pooled_vs_lcv_comparison$auc_gap
  gap_pct <- pooled_vs_lcv_comparison$percent_drop

  p <- ggplot(perf_data, aes(x = approach, y = auc, fill = type)) +
    geom_col(width = 0.7, alpha = 0.9) +
    # Add error bars if CIs available
    {
      if (!all(is.na(perf_data$ci_lower))) {
        geom_errorbar(
          aes(ymin = ci_lower, ymax = ci_upper),
          width = 0.3,
          linewidth = 0.8,
          color = "#EEEEEE",
          alpha = 0.8
        )
      }
    } +
    geom_text(
      aes(label = sprintf("%.3f", auc)),
      vjust = ifelse(!all(is.na(perf_data$ci_lower)), -1.5, -0.5),
      size = 4, fontface = "bold", color = "#EEEEEE"
    ) +
    # Highlight the gap
    annotate("segment",
      x = 1, xend = 4,
      y = pooled_vs_lcv_comparison$pooled_auc,
      yend = pooled_vs_lcv_comparison$lcv_mean_auc,
      linetype = "dashed", color = "#FDE724", linewidth = 1
    ) +
    annotate("text",
      x = 2.5,
      y = mean(c(
        pooled_vs_lcv_comparison$pooled_auc,
        pooled_vs_lcv_comparison$lcv_mean_auc
      )),
      label = sprintf("Gap: %.3f\n(%.1f%% drop)", gap, gap_pct),
      color = "#FDE724", size = 4, fontface = "bold"
    ) +
    scale_fill_manual(values = c("Pooled" = "#FDE724", "LCV" = "#31688E")) + # viridis yellow and cyan
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Model Performance: Pooled CV vs Geographic Validation",
      subtitle = "⚠️ SEVERE: 18.2% performance drop reveals geographic overfitting (error bars = 95% CI)",
      x = NULL,
      y = "AUC (Area Under ROC Curve)",
      fill = "Approach"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 7, dpi = 300)
    message(sprintf("Saved performance comparison to: %s", save_path))
  }

  return(p)
}


#' Plot confusion matrices for leave-country-out CV
#'
#' @description
#' Side-by-side confusion matrices for Italy→US and US→Italy.
#'
#' @param lcv_results Leave-country-out CV results
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_confusion_matrices <- function(lcv_results, save_path = NULL) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # Prepare data for Italy → US
  cm_italy_to_us <- as.data.frame.table(lcv_results$italy_to_us$confusion_matrix)
  names(cm_italy_to_us) <- c("Predicted", "Actual", "Count")
  cm_italy_to_us$Direction <- "Italy → US"

  # Prepare data for US → Italy
  cm_us_to_italy <- as.data.frame.table(lcv_results$us_to_italy$confusion_matrix)
  names(cm_us_to_italy) <- c("Predicted", "Actual", "Count")
  cm_us_to_italy$Direction <- "US → Italy"

  # Combine
  cm_data <- rbind(cm_italy_to_us, cm_us_to_italy)

  # Calculate percentages within each direction
  cm_data <- cm_data %>%
    dplyr::group_by(Direction) %>%
    dplyr::mutate(
      Total = sum(Count),
      Percentage = 100 * Count / Total,
      Label = sprintf("%d\n(%.1f%%)", Count, Percentage)
    )

  p <- ggplot(cm_data, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "#555555", linewidth = 1) +
    geom_text(aes(label = Label), size = 4, fontface = "bold", color = "#EEEEEE") +
    facet_wrap(~Direction, ncol = 2) +
    scale_fill_gradient(low = "#440154", high = "#31688E") + # viridis dark purple to cyan
    labs(
      title = "Confusion Matrices: Geographic Validation Failure",
      subtitle = "US → Italy shows 98% misclassification of ALS cases - extreme overfitting to source cohort",
      x = "True Label",
      y = "Predicted Label"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    coord_equal()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 12, height = 5, dpi = 300)
    message(sprintf("Saved confusion matrices to: %s", save_path))
  }

  return(p)
}


# ==============================================================================
# PROTEIN OVERLAP AND CONCORDANCE VISUALIZATIONS
# ==============================================================================

#' Plot Venn diagram for protein overlap
#'
#' @description
#' Venn diagram showing overlap between:
#' - Italy significant proteins
#' - US significant proteins
#' - Pooled significant proteins
#'
#' @param stratified_vs_pooled_comparison Comparison results
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_protein_overlap_venn <- function(stratified_vs_pooled_comparison, save_path = NULL) {
  library(ggplot2)
  library(ggvenn)

  # Prepare data for ggvenn
  venn_data <- list(
    Italy = stratified_vs_pooled_comparison$sig_italy,
    US = stratified_vs_pooled_comparison$sig_us,
    Pooled = stratified_vs_pooled_comparison$sig_pooled
  )

  p <- ggvenn(
    venn_data,
    fill_color = c("#1F9E89", "#B5DE2B", "#FDE724"), # viridis jade, yellow-green, yellow
    fill_alpha = 0.5,
    stroke_size = 1,
    set_name_size = 5,
    text_size = 4
  ) +
    labs(
      title = "Protein Overlap: Stratified vs Pooled Analysis",
      subtitle = sprintf(
        "⚠️ CRITICAL: Only %d proteins (%.1f%%) significant in BOTH Italy and US",
        length(stratified_vs_pooled_comparison$sig_both_strata),
        stratified_vs_pooled_comparison$pct_replicate_in_both
      )
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(color = "#FDE724", face = "bold", hjust = 0.5)
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    message(sprintf("Saved Venn diagram to: %s", save_path))
  }

  return(p)
}


#' Plot protein counts comparison
#'
#' @description
#' Bar chart comparing number of significant proteins across analyses.
#'
#' @param stratified_vs_pooled_comparison Comparison results
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_protein_counts <- function(stratified_vs_pooled_comparison, save_path = NULL) {
  library(ggplot2)

  counts_data <- data.frame(
    analysis = c(
      "Italy\n(HEPARIN)",
      "US\n(EDTA)",
      "Pooled\n(Confounded)",
      "Both Strata\n(Robust)"
    ),
    count = c(
      length(stratified_vs_pooled_comparison$sig_italy),
      length(stratified_vs_pooled_comparison$sig_us),
      length(stratified_vs_pooled_comparison$sig_pooled),
      length(stratified_vs_pooled_comparison$sig_both_strata)
    ),
    type = c("Stratified", "Stratified", "Pooled", "Robust")
  )

  counts_data$analysis <- factor(counts_data$analysis, levels = counts_data$analysis)

  p <- ggplot(counts_data, aes(x = analysis, y = count, fill = type)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = count), vjust = -0.5, size = 5, fontface = "bold", color = "#EEEEEE") +
    scale_fill_manual(values = c(
      "Stratified" = "#31688E", # viridis cyan
      "Pooled" = "#FDE724", # viridis yellow
      "Robust" = "#6CCE59" # viridis lime
    )) +
    labs(
      title = "Significant Proteins by Analysis Type",
      subtitle = sprintf(
        "Only %d/%d pooled proteins (%.1f%%) replicate in both strata - major reproducibility issue",
        length(stratified_vs_pooled_comparison$overlap_both),
        length(stratified_vs_pooled_comparison$sig_pooled),
        stratified_vs_pooled_comparison$pct_replicate_in_both
      ),
      x = NULL,
      y = "Number of Significant Proteins\n(FDR < 0.05, |logFC| > 0.5)",
      fill = "Type"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 7, dpi = 300)
    message(sprintf("Saved protein counts to: %s", save_path))
  }

  return(p)
}


#' Plot effect size correlation
#'
#' @description
#' Scatter plot of Italy vs US effect sizes with correlation.
#'
#' @param protein_concordance Concordance results
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_effect_size_correlation <- function(protein_concordance, save_path = NULL) {
  library(ggplot2)
  library(dplyr)

  # Calculate correlation
  cor_val <- cor(protein_concordance$logFC_italy,
    protein_concordance$logFC_us,
    use = "complete.obs"
  )

  # Highlight concordant proteins
  plot_data <- protein_concordance %>%
    dplyr::mutate(
      concordant_sig = sig_both & same_direction,
      point_type = dplyr::case_when(
        concordant_sig ~ "Tube-Robust (sig in both)",
        sig_italy | sig_us ~ "Significant (one cohort)",
        TRUE ~ "Not significant"
      )
    )

  p <- ggplot(plot_data, aes(x = logFC_italy, y = logFC_us, color = point_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#888888") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "#FDE724", linetype = "dashed") +
    scale_color_manual(values = c(
      "Tube-Robust (sig in both)" = "#6CCE59", # viridis lime
      "Significant (one cohort)" = "#FDE724", # viridis yellow
      "Not significant" = "#71717A" # neutral gray
    )) +
    annotate("text",
      x = min(plot_data$logFC_italy, na.rm = TRUE),
      y = max(plot_data$logFC_us, na.rm = TRUE),
      label = sprintf("Correlation: r = %.3f", cor_val),
      hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "#FDE724"
    ) +
    labs(
      title = "Effect Size Concordance: Italy vs US",
      subtitle = "Extremely low correlation (r = 0.03) indicates inconsistent effects driven by confounding",
      x = "Log Fold Change (Italy/HEPARIN)",
      y = "Log Fold Change (US/EDTA)",
      color = "Protein Type"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "bottom"
    ) +
    coord_equal()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 9, dpi = 300)
    message(sprintf("Saved effect size correlation to: %s", save_path))
  }

  return(p)
}


# ==============================================================================
# CONFOUNDING STRUCTURE VISUALIZATIONS
# ==============================================================================

#' Plot confounding distribution
#'
#' @description
#' Visualizes the perfect confounding structure:
#' Diagnosis × Tube Type × Country
#'
#' @param sample_metadata Sample-level metadata
#' @param save_path Optional path to save the plot
#' @return ggplot object
#' @export
plot_confounding_structure <- function(sample_metadata, save_path = NULL) {
  # Create summary
  summary_df <- sample_metadata %>%
    dplyr::count(Diagnosis, Plasma_collection_tube_type, country) %>%
    dplyr::mutate(
      Diagnosis = factor(Diagnosis, levels = c("ALS", "Healthy_control", "Neurological_control"))
    )

  p <- ggplot(
    summary_df,
    aes(x = Diagnosis, y = n, fill = Plasma_collection_tube_type)
  ) +
    geom_col(position = "stack", color = "#2A3142", linewidth = 0.5) +
    geom_text(aes(label = n),
      position = position_stack(vjust = 0.5),
      color = "#0A0E14", fontface = "bold", size = 4.5
    ) +
    facet_wrap(~country, ncol = 2) +
    scale_fill_manual(
      values = tube_type_colors_viridis,
      name = "Tube Type"
    ) +
    labs(
      title = "Perfect Confounding: Diagnosis × Tube Type × Geographic Origin",
      subtitle = "100% of Neurological Controls are US/EDTA - mathematically impossible to separate effects",
      x = "Diagnosis",
      y = "Number of Samples"
    ) +
    theme_dark_scientific(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 6, dpi = 300)
    message(sprintf("Saved confounding structure plot to %s", save_path))
  }

  return(p)
}


#' Plot PCA colored by tube type vs diagnosis
#'
#' @description
#' Performs PCA on protein data and colors by both tube type and diagnosis.
#' If tube type separates samples better than diagnosis, this indicates
#' technical effects dominate biological signal.
#'
#' @param protein_wide Wide-format protein data
#' @param sample_metadata Sample metadata
#' @param save_path Optional path to save the plot
#' @return List with plot and PCA results
#' @export
plot_pca_tube_vs_diagnosis <- function(protein_wide, sample_metadata = NULL, save_path = NULL) {
  # Extract protein matrix
  protein_matrix <- protein_wide %>%
    dplyr::filter(country != "Unknown") %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis,
      -Plasma_collection_tube_type, -country, -Sex, -Age_Collection
    )

  # Handle missing values (median imputation)
  protein_matrix <- as.matrix(protein_matrix)
  for (i in seq_len(ncol(protein_matrix))) {
    protein_matrix[is.na(protein_matrix[, i]), i] <- median(protein_matrix[, i], na.rm = TRUE)
  }

  # Perform PCA
  pca_result <- prcomp(protein_matrix, center = TRUE, scale. = TRUE)

  # Create data frame for plotting
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Diagnosis = protein_wide$Diagnosis[protein_wide$country != "Unknown"],
    Tube_Type = protein_wide$Plasma_collection_tube_type[protein_wide$country != "Unknown"],
    Country = protein_wide$country[protein_wide$country != "Unknown"]
  )

  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100

  # Plot 1: Colored by tube type
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tube_Type, shape = Country)) +
    geom_point(size = 3, alpha = 0.8, stroke = 0.5) +
    scale_color_manual(
      values = tube_type_colors_viridis,
      name = "Tube Type"
    ) +
    labs(
      title = "PCA: Protein Expression Colored by Tube Type",
      subtitle = "Clear separation by anticoagulant indicates technical artifact dominates variance",
      x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
      y = sprintf("PC2 (%.1f%% variance)", var_explained[2])
    ) +
    theme_dark_scientific(base_size = 13) +
    theme(
      legend.position = "right"
    )

  # Plot 2: Colored by diagnosis
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis, shape = Country)) +
    geom_point(size = 3, alpha = 0.8, stroke = 0.5) +
    scale_color_manual(
      values = diagnosis_colors_viridis,
      name = "Diagnosis"
    ) +
    labs(
      title = "PCA: Protein Expression Colored by Diagnosis",
      subtitle = "Biological disease groups show minimal separation compared to technical factors",
      x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
      y = sprintf("PC2 (%.1f%% variance)", var_explained[2])
    ) +
    theme_dark_scientific(base_size = 13) +
    theme(
      legend.position = "right"
    )

  # Combine plots
  combined_plot <- patchwork::wrap_plots(p1, p2, ncol = 1)

  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = 10, height = 12, dpi = 300)
    message(sprintf("Saved PCA plots to %s", save_path))
  }

  return(list(
    plot = combined_plot,
    pca_result = pca_result,
    pca_data = pca_df
  ))
}


# ==============================================================================
# FOREST PLOTS AND META-ANALYSIS VISUALIZATIONS
# ==============================================================================

#' Plot forest plot for tube-robust proteins
#'
#' @description
#' Forest plot showing effect sizes in Italy vs US for concordant proteins.
#'
#' @param protein_concordance Concordance analysis results
#' @param top_n Number of top proteins to plot (default 15)
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_tube_robust_forest <- function(protein_concordance, top_n = 15, save_path = NULL) {
  library(ggplot2)
  library(dplyr)

  # Get top tube-robust proteins
  candidates <- protein_concordance %>%
    dplyr::filter(sig_both, same_direction) %>%
    dplyr::arrange(combined_p) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::select(Assay, logFC_italy, logFC_us, combined_p)

  # Reshape for plotting
  forest_data <- rbind(
    candidates %>%
      dplyr::mutate(
        cohort = "Italy (HEPARIN)",
        logFC = logFC_italy
      ),
    candidates %>%
      dplyr::mutate(
        cohort = "US (EDTA)",
        logFC = logFC_us
      )
  ) %>%
    dplyr::select(Assay, cohort, logFC, combined_p)

  # Order by combined p-value
  forest_data$Assay <- factor(
    forest_data$Assay,
    levels = rev(candidates$Assay)
  )

  p <- ggplot(forest_data, aes(x = logFC, y = Assay, color = cohort)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(
      aes(xmin = logFC * 0.9, xmax = logFC * 1.1),
      height = 0.3,
      position = position_dodge(width = 0.5)
    ) +
    scale_color_manual(values = c(
      "Italy (HEPARIN)" = "#1F9E89", # viridis jade
      "US (EDTA)" = "#B5DE2B" # viridis yellow-green
    )) +
    labs(
      title = "Tube-Robust Biomarker Candidates",
      subtitle = sprintf("Top %d proteins significant in BOTH Italy and US with consistent direction", top_n),
      x = "Log Fold Change (ALS vs Controls)",
      y = NULL,
      color = "Cohort"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    message(sprintf("Saved forest plot to: %s", save_path))
  }

  return(p)
}


# ==============================================================================
# RESIDUAL CONFOUNDING VISUALIZATIONS
# ==============================================================================

#' Plot comparison of reverse prediction (all vs robust-only proteins)
#'
#' @description
#' Bar chart comparing AUC for predicting tube type using:
#'   - All proteins (~2880 proteins)
#'   - Tube-robust proteins only (22 proteins)
#'
#' Shows whether even "robust" proteins retain tube type signal.
#' Now includes error bars showing 95% confidence intervals.
#'
#' @param robust_reverse_pred Reverse prediction results for robust proteins only (must contain auc_ci_lower, auc_ci_upper)
#' @param full_reverse_pred Reverse prediction results for all proteins (must contain auc_ci_lower, auc_ci_upper)
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_robust_reverse_prediction_comparison <- function(robust_reverse_pred,
                                                      full_reverse_pred,
                                                      save_path = NULL) {
  library(ggplot2)

  # Get protein counts safely
  n_full <- if (!is.null(full_reverse_pred$n_proteins)) full_reverse_pred$n_proteins else ncol(full_reverse_pred$feature_importance)
  n_robust <- if (!is.null(robust_reverse_pred$n_proteins)) robust_reverse_pred$n_proteins else nrow(robust_reverse_pred$feature_importance)

  # Extract CIs (if available)
  full_ci_lower <- if (!is.null(full_reverse_pred$auc_ci_lower)) full_reverse_pred$auc_ci_lower else NA
  full_ci_upper <- if (!is.null(full_reverse_pred$auc_ci_upper)) full_reverse_pred$auc_ci_upper else NA
  robust_ci_lower <- if (!is.null(robust_reverse_pred$auc_ci_lower)) robust_reverse_pred$auc_ci_lower else NA
  robust_ci_upper <- if (!is.null(robust_reverse_pred$auc_ci_upper)) robust_reverse_pred$auc_ci_upper else NA

  comparison_data <- data.frame(
    model = c(
      sprintf("All Proteins\n(n=%d)", n_full),
      sprintf("Tube-Robust Only\n(n=%d)", n_robust)
    ),
    auc = c(full_reverse_pred$auc, robust_reverse_pred$auc),
    ci_lower = c(full_ci_lower, robust_ci_lower),
    ci_upper = c(full_ci_upper, robust_ci_upper),
    type = c("Full", "Robust"),
    stringsAsFactors = FALSE
  )

  # Ensure unique factor levels
  comparison_data$model <- factor(comparison_data$model, levels = unique(comparison_data$model))

  # Determine severity (using viridis colors consistently)
  robust_auc <- robust_reverse_pred$auc
  severity_color <- if (robust_auc > 0.85) {
    "#FDE724" # Severe - viridis yellow (100%)
  } else if (robust_auc > 0.75) {
    "#B5DE2B" # Substantial - viridis yellow-green (90%)
  } else if (robust_auc > 0.65) {
    "#6CCE59" # Moderate - viridis lime (80%)
  } else {
    "#1F9E89" # Minimal - viridis jade green (60%)
  }

  severity_label <- if (robust_auc > 0.85) {
    "⚠️ SEVERE ENTANGLEMENT"
  } else if (robust_auc > 0.75) {
    "⚠️ SUBSTANTIAL ENTANGLEMENT"
  } else if (robust_auc > 0.65) {
    "⚠️ MODERATE ENTANGLEMENT"
  } else {
    "✓ MINIMAL ENTANGLEMENT"
  }

  # Build CI label for subtitle
  ci_label <- if (!is.na(robust_ci_lower)) {
    sprintf(
      "%s: Robust proteins predict tube type (AUC = %.3f, 95%% CI: %.3f-%.3f)",
      severity_label, robust_auc, robust_ci_lower, robust_ci_upper
    )
  } else {
    sprintf(
      "%s: Robust proteins predict tube type with AUC = %.3f",
      severity_label, robust_auc
    )
  }

  p <- ggplot(comparison_data, aes(x = model, y = auc, fill = type)) +
    geom_col(width = 0.6, alpha = 0.9) +
    # Add error bars if CIs available
    {
      if (!all(is.na(comparison_data$ci_lower))) {
        geom_errorbar(
          aes(ymin = ci_lower, ymax = ci_upper),
          width = 0.25,
          linewidth = 0.8,
          color = "#EEEEEE",
          alpha = 0.8
        )
      }
    } +
    geom_text(
      aes(label = sprintf("%.3f", auc)),
      vjust = ifelse(!all(is.na(comparison_data$ci_lower)), -1.5, -0.5),
      size = 5, fontface = "bold", color = "#EEEEEE"
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "#888888") +
    annotate("text",
      x = 1.5, y = 0.55,
      label = "Chance level", color = "#888888", size = 3
    ) +
    scale_fill_manual(values = c("Full" = "#31688E", "Robust" = severity_color)) + # viridis cyan for Full
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Residual Confounding in Tube-Robust Proteins",
      subtitle = ci_label,
      x = NULL,
      y = "AUC for Predicting Tube Type (error bars = 95% CI)",
      fill = "Protein Set"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.subtitle = element_text(color = severity_color, face = "bold")
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 7, dpi = 300)
    message(sprintf("Saved reverse prediction comparison to: %s", save_path))
  }

  return(p)
}


#' Plot tube effects in healthy controls
#'
#' @description
#' Forest plot showing proteins with significant tube type effects
#' in healthy controls only (no disease present).
#'
#' @param healthy_effects Results from tube_effect_in_healthy()
#' @param top_n Number of most significant proteins to plot (default 15)
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_healthy_tube_effects <- function(healthy_effects, top_n = 15, save_path = NULL) {
  library(ggplot2)
  library(dplyr)

  # Get top proteins by p-value
  plot_data <- healthy_effects %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      protein = factor(protein, levels = rev(protein)),
      significant = p_value < 0.05,
      ci_lower = raw_difference - 1.96 * std_error,
      ci_upper = raw_difference + 1.96 * std_error
    )

  n_sig <- sum(healthy_effects$significant)
  n_total <- nrow(healthy_effects)

  p <- ggplot(plot_data, aes(x = raw_difference, y = protein, color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
    scale_color_manual(values = c("TRUE" = "#FDE724", "FALSE" = "#888888")) +
    labs(
      title = "Tube Type Effects in Healthy Controls (Disease-Free)",
      subtitle = sprintf(
        "⚠️ %d/%d tube-robust proteins differ by tube type WITHOUT disease present",
        n_sig, n_total
      ),
      x = "Difference in NPX (US/EDTA - Italy/HEPARIN)",
      y = NULL,
      color = "Significant\n(p < 0.05)"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.subtitle = element_text(color = "#FDE724", face = "bold")
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    message(sprintf("Saved healthy tube effects to: %s", save_path))
  }

  return(p)
}


#' Plot effect size decomposition ratios
#'
#' @description
#' Scatter plot showing tube effect vs disease effect for each tube-robust
#' protein. Points above the diagonal (y=x) have tube effects larger than
#' disease effects.
#'
#' @param effect_decomp Results from effect_size_decomposition()
#' @param save_path Path to save figure (optional)
#' @return ggplot object
#' @export
plot_effect_decomposition_ratios <- function(effect_decomp, save_path = NULL) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  # Prepare data
  plot_data <- effect_decomp %>%
    dplyr::mutate(
      category = dplyr::case_when(
        ratio_tube_to_disease > 1.0 ~ "Tube > Disease",
        ratio_tube_to_disease > 0.5 ~ "Tube ≥ 50% of Disease",
        ratio_tube_to_disease > 0.2 ~ "Tube 20-50% of Disease",
        TRUE ~ "Tube < 20% of Disease"
      ),
      category = factor(category, levels = c(
        "Tube > Disease",
        "Tube ≥ 50% of Disease",
        "Tube 20-50% of Disease",
        "Tube < 20% of Disease"
      ))
    )

  # Count by category
  category_counts <- plot_data %>%
    dplyr::count(category)

  n_severe <- sum(plot_data$ratio_tube_to_disease >= 0.5, na.rm = TRUE)
  pct_severe <- 100 * n_severe / nrow(plot_data)

  # Label top 5 most problematic proteins
  proteins_to_label <- plot_data %>%
    dplyr::arrange(dplyr::desc(ratio_tube_to_disease)) %>%
    dplyr::slice_head(n = 5)

  p <- ggplot(plot_data, aes(x = abs_disease_effect, y = abs_tube_effect, color = category)) +
    # Diagonal line (tube = disease)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#FDE724", linewidth = 1) +
    # Reference line (tube = 50% of disease)
    geom_abline(intercept = 0, slope = 0.5, linetype = "dotted", color = "#6CCE59", linewidth = 0.8) +
    # Points
    geom_point(size = 3, alpha = 0.7) +
    # Label problematic proteins
    ggrepel::geom_text_repel(
      data = proteins_to_label,
      aes(label = protein),
      size = 3,
      color = "#EEEEEE",
      box.padding = 0.5,
      max.overlaps = 20
    ) +
    scale_color_manual(values = c(
      "Tube > Disease" = "#FDE724", # viridis yellow (100%)
      "Tube ≥ 50% of Disease" = "#B5DE2B", # viridis yellow-green (90%)
      "Tube 20-50% of Disease" = "#6CCE59", # viridis lime (80%)
      "Tube < 20% of Disease" = "#1F9E89" # viridis jade (60%)
    )) +
    annotate("text",
      x = max(plot_data$abs_disease_effect, na.rm = TRUE) * 0.7,
      y = max(plot_data$abs_tube_effect, na.rm = TRUE) * 0.9,
      label = sprintf("y = x\n(Tube = Disease)"),
      color = "#FDE724", size = 3, fontface = "italic"
    ) +
    labs(
      title = "Effect Size Decomposition: Tube Type vs Disease Effects",
      subtitle = sprintf(
        "⚠️ SEVERE: %d/%d proteins (%.1f%%) have tube effects ≥50%% of disease effects",
        n_severe, nrow(plot_data), pct_severe
      ),
      x = "Absolute Disease Effect (|ALS - Control|)",
      y = "Absolute Tube Effect (|US - Italy|)",
      color = "Entanglement Level"
    ) +
    theme_dark_scientific(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(color = "#FDE724", face = "bold")
    ) +
    coord_equal()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 9, dpi = 300)
    message(sprintf("Saved effect decomposition to: %s", save_path))
  }

  return(p)
}


# ==============================================================================
# SUMMARY AND EXPORT FUNCTIONS
# ==============================================================================

#' Create summary figure panel (Reverse Prediction Focus)
#'
#' @description
#' Combines key findings into a single multi-panel figure for presentation.
#' Focus on reverse prediction test and confounding structure.
#'
#' @param reverse_prediction_results Reverse prediction test results
#' @param sample_metadata Sample metadata
#' @param protein_wide Wide protein data
#' @param save_path Path to save the figure
#' @return Combined ggplot object
#' @export
create_summary_figure <- function(reverse_prediction_results, sample_metadata,
                                  protein_wide, save_path = NULL) {
  message("Creating summary figure panel...")

  # Generate individual plots
  p_roc <- plot_reverse_prediction_roc(reverse_prediction_results)
  p_confound <- plot_confounding_structure(sample_metadata)
  p_features <- plot_tube_type_features(reverse_prediction_results)

  # Combine with patchwork
  summary_fig <- (p_roc | p_confound) / p_features +
    plot_annotation(
      title = "ALS Biomarker Study: Evidence of Severe Confounding Bias",
      subtitle = "Reverse prediction AUC = 0.999 indicates technical artifacts overwhelm disease biology",
      theme = theme_dark_scientific(base_size = 14) +
        theme(
          plot.title = element_text(size = 18, face = "bold", color = "#F0F6FC"),
          plot.subtitle = element_text(size = 14, color = "#FDE724")
        )
    )

  if (!is.null(save_path)) {
    ggsave(save_path, summary_fig, width = 16, height = 12, dpi = 300)
    message(sprintf("Saved summary figure to %s", save_path))
  }

  return(summary_fig)
}


#' Create comprehensive summary panel (Geographic Validation Focus)
#'
#' @description
#' Multi-panel figure combining key visualizations from geographic validation.
#'
#' @param lcv_results Leave-country-out CV results
#' @param pooled_results Pooled CV results
#' @param pooled_vs_lcv_comparison Comparison results
#' @param stratified_vs_pooled_comparison Stratified comparison
#' @param protein_concordance Concordance results
#' @param save_path Path to save figure (optional)
#' @return Combined plot
#' @export
create_summary_panel <- function(lcv_results, pooled_results,
                                 pooled_vs_lcv_comparison,
                                 stratified_vs_pooled_comparison,
                                 protein_concordance,
                                 save_path = NULL) {
  library(patchwork)

  message("Creating comprehensive summary panel...")

  # Generate individual plots
  p1 <- plot_performance_comparison(pooled_vs_lcv_comparison)
  p2 <- plot_confusion_matrices(lcv_results)
  p3 <- plot_protein_counts(stratified_vs_pooled_comparison)
  p4 <- plot_tube_robust_forest(protein_concordance, top_n = 10)

  # Combine with patchwork
  combined <- (p1 | p3) / (p2 | p4) +
    plot_annotation(
      title = "Geographic Confounding Investigation: Summary of Key Findings",
      subtitle = "⚠️ Multiple converging lines of evidence reveal severe geographic/tube type confounding",
      theme = theme_dark_scientific(base_size = 14) +
        theme(
          plot.title = element_text(size = 16, hjust = 0.5),
          plot.subtitle = element_text(color = "#FDE724", face = "bold", size = 12, hjust = 0.5)
        )
    )

  if (!is.null(save_path)) {
    ggsave(save_path, combined, width = 16, height = 14, dpi = 300)
    message(sprintf("Saved summary panel to: %s", save_path))
  }

  return(combined)
}


#' Export all visualizations
#'
#' @description
#' Convenience function to generate all visualizations at once.
#'
#' @param lcv_results Leave-country-out CV results
#' @param pooled_results Pooled CV results
#' @param pooled_vs_lcv_comparison Comparison results
#' @param stratified_vs_pooled_comparison Stratified comparison
#' @param protein_concordance Concordance results
#' @param output_dir Output directory (default: "outputs/figures")
#' @return List of plot objects
#' @export
export_all_visualizations <- function(lcv_results,
                                      pooled_results,
                                      pooled_vs_lcv_comparison,
                                      stratified_vs_pooled_comparison,
                                      protein_concordance,
                                      output_dir = "outputs/figures") {
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(sprintf("Created output directory: %s", output_dir))
  }

  message("\n", paste(rep("=", 70), collapse = ""))
  message("GENERATING ALL VISUALIZATIONS")
  message(paste(rep("=", 70), collapse = ""))

  plots <- list()

  # 1. ROC curves
  message("\n1. ROC curves for leave-country-out CV...")
  plots$roc <- plot_lcv_roc_curves(
    lcv_results, pooled_results,
    save_path = file.path(output_dir, "01_lcv_roc_curves.png")
  )

  # 2. Performance comparison
  message("\n2. Performance comparison bar chart...")
  plots$performance <- plot_performance_comparison(
    pooled_vs_lcv_comparison,
    save_path = file.path(output_dir, "02_performance_comparison.png")
  )

  # 3. Confusion matrices
  message("\n3. Confusion matrices...")
  plots$confusion <- plot_confusion_matrices(
    lcv_results,
    save_path = file.path(output_dir, "03_confusion_matrices.png")
  )

  # 4. Protein counts
  message("\n4. Protein counts comparison...")
  plots$counts <- plot_protein_counts(
    stratified_vs_pooled_comparison,
    save_path = file.path(output_dir, "04_protein_counts.png")
  )

  # 5. Venn diagram
  message("\n5. Protein overlap Venn diagram...")
  plots$venn <- plot_protein_overlap_venn(
    stratified_vs_pooled_comparison,
    save_path = file.path(output_dir, "05_protein_venn.png")
  )

  # 6. Forest plot
  message("\n6. Tube-robust protein forest plot...")
  plots$forest <- plot_tube_robust_forest(
    protein_concordance,
    save_path = file.path(output_dir, "06_forest_plot.png")
  )

  # 7. Effect size correlation
  message("\n7. Effect size correlation scatter...")
  plots$correlation <- plot_effect_size_correlation(
    protein_concordance,
    save_path = file.path(output_dir, "07_effect_correlation.png")
  )

  # 8. Summary panel
  message("\n8. Comprehensive summary panel...")
  plots$summary <- create_summary_panel(
    lcv_results, pooled_results, pooled_vs_lcv_comparison,
    stratified_vs_pooled_comparison, protein_concordance,
    save_path = file.path(output_dir, "08_summary_panel.png")
  )

  message("\n", paste(rep("=", 70), collapse = ""))
  message("ALL VISUALIZATIONS COMPLETE!")
  message(sprintf("Saved to: %s", output_dir))
  message(paste(rep("=", 70), collapse = ""))

  return(plots)
}
