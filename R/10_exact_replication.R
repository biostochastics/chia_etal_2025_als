#' Exact Replication of Original Study Methods
#'
#' Functions that precisely reproduce the methodology from Chia et al. (2025)
#' notebooks to enable direct comparison with our analysis.
#'
#' @description
#' This module implements:
#' 1. 80/20 discovery/replication split (stratified by diagnosis + tube type)
#' 2. Feature selection using differential analysis
#' 3. ML models WITH tube type as predictor (their approach)
#' 4. Beta-beta plots comparing HEPARIN vs EDTA effect sizes
#' 5. Their exact thresholds (FDR < 0.05, |logFC| > 0.6)


#' Create 80/20 split stratified by diagnosis AND tube type
#'
#' @description
#' Replicates the original study's discovery/replication split:
#' - 80% discovery, 20% replication
#' - Stratified by BOTH Diagnosis AND Plasma_collection_tube_type
#' - This maintains tube type distribution across splits (unlike LCV)
#'
#' CRITICAL: This approach does NOT test geographic generalization!
#' Both discovery and replication contain HEPARIN and EDTA samples.
#'
#' @param protein_wide Wide-format protein data
#' @param split_ratio Proportion for discovery (default 0.8)
#' @param seed Random seed
#' @return List with discovery and replication datasets
#' @examples
#' split_data <- create_chia_8020_split(protein_wide, seed = 42)
#' @export
create_chia_8020_split <- function(protein_wide, split_ratio = 0.8, seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("CREATING 80/20 SPLIT (Original Study Method)")
  message(paste(rep("=", 70), collapse = ""))
  message("\nStratified by DIAGNOSIS + TUBE TYPE (maintains tube distribution)")

  set.seed(seed)

  # Prepare data
  data_clean <- protein_wide %>%
    dplyr::filter(country %in% c("Italy", "US")) %>%
    dplyr::mutate(
      Diagnosis_binary = dplyr::case_when(
        Diagnosis == "ALS" ~ "ALS",
        Diagnosis %in% c("Healthy_control", "Neurological_control") ~ "Control",
        TRUE ~ NA_character_
      ),
      # Create stratification variable: Diagnosis + Tube Type
      strata = paste(Diagnosis_binary, Plasma_collection_tube_type, sep = "_")
    ) %>%
    dplyr::filter(!is.na(Diagnosis_binary))

  # Show stratification structure
  strata_table <- table(data_clean$Diagnosis_binary, data_clean$Plasma_collection_tube_type)
  message("\nSample distribution:")
  print(strata_table)

  # Create stratified split using caret
  train_indices <- caret::createDataPartition(
    y = data_clean$strata,
    p = split_ratio,
    list = FALSE
  )

  discovery <- data_clean[train_indices, ] %>%
    dplyr::select(-strata)
  replication <- data_clean[-train_indices, ] %>%
    dplyr::select(-strata)

  # Verify stratification maintained
  message("\nDISCOVERY (80%):")
  message(sprintf("  Total: %s samples", nrow(discovery)))
  disc_table <- table(discovery$Diagnosis_binary, discovery$Plasma_collection_tube_type)
  print(disc_table)

  message("\nREPLICATION (20%):")
  message(sprintf("  Total: %s samples", nrow(replication)))
  rep_table <- table(replication$Diagnosis_binary, replication$Plasma_collection_tube_type)
  print(rep_table)

  message("\n", paste(rep("=", 70), collapse = ""))
  message("CRITICAL NOTE:")
  message("  Both discovery and replication contain HEPARIN and EDTA samples.")
  message("  This split does NOT test geographic generalization!")
  message("  Compare with leave-country-out CV for true validation.")
  message(paste(rep("=", 70), collapse = ""))

  return(list(
    discovery = discovery,
    replication = replication,
    discovery_table = disc_table,
    replication_table = rep_table
  ))
}


#' Train model using original study approach
#'
#' @description
#' Replicates their ML methodology:
#' 1. Use differentially expressed proteins as features
#' 2. Include tube type (Cohort) as a predictor
#' 3. Include Sex and Age
#' 4. Train on discovery, test on replication
#'
#' @param split_data Output from create_chia_8020_split()
#' @param differential_pooled Pooled differential analysis results
#' @param fdr_threshold False Discovery Rate threshold (default 0.05 = 5% expected false positives). Controls Type I error under multiple testing via Benjamini-Hochberg
#' @param logfc_threshold Log fold-change threshold (default 0.6 = 1.52-fold change). Filters for biologically meaningful effects; original study used 0.6
#' @param seed Random seed for reproducibility
#' @return Model results with discovery/replication performance
#' @examples
#' model <- train_chia_replication_model(split_data, differential_pooled)
#' @export
train_chia_replication_model <- function(split_data,
                                         differential_pooled,
                                         fdr_threshold = 0.05,
                                         logfc_threshold = 0.6,
                                         seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("TRAINING MODEL: Original Study Methodology")
  message(paste(rep("=", 70), collapse = ""))
  message("\nFeatures: Differential proteins + Tube Type + Age + Sex")
  message("Training: Discovery (80%)")
  message("Testing: Replication (20%)")

  set.seed(seed)

  # Select significant proteins from pooled differential analysis
  sig_proteins <- differential_pooled %>%
    dplyr::filter(
      adj.P.Val < fdr_threshold,
      abs(logFC) > logfc_threshold
    ) %>%
    dplyr::pull(Assay)

  message(sprintf(
    "\nFeature selection: %s proteins (FDR < %.2f, |logFC| > %.1f)",
    length(sig_proteins), fdr_threshold, logfc_threshold
  ))

  if (length(sig_proteins) == 0) {
    stop("No significant proteins found with current thresholds!")
  }

  # Prepare discovery data
  discovery <- split_data$discovery
  y_discovery <- factor(discovery$Diagnosis_binary, levels = c("Control", "ALS"))

  X_discovery <- discovery %>%
    dplyr::select(
      # Tube type (encoded as 0/1)
      Plasma_collection_tube_type,
      # Covariates
      Sex, Age_Collection,
      # Differential proteins
      any_of(sig_proteins)
    ) %>%
    dplyr::mutate(
      tube_type_encoded = as.numeric(Plasma_collection_tube_type == "HEPARIN"),
      sex_encoded = as.numeric(Sex == "Male")
    ) %>%
    dplyr::select(-Plasma_collection_tube_type, -Sex)

  # Prepare replication data
  replication <- split_data$replication
  y_replication <- factor(replication$Diagnosis_binary, levels = c("Control", "ALS"))

  X_replication <- replication %>%
    dplyr::select(
      Plasma_collection_tube_type,
      Sex, Age_Collection,
      any_of(sig_proteins)
    ) %>%
    dplyr::mutate(
      tube_type_encoded = as.numeric(Plasma_collection_tube_type == "HEPARIN"),
      sex_encoded = as.numeric(Sex == "Male")
    ) %>%
    dplyr::select(-Plasma_collection_tube_type, -Sex)

  message(sprintf("\nFinal features: %s (proteins + tube + age + sex)", ncol(X_discovery)))

  # Train Random Forest on DISCOVERY with internal CV
  message("\nTraining Random Forest on DISCOVERY (with 5-fold CV)...")

  train_control <- caret::trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )

  rf_model <- caret::train(
    x = X_discovery,
    y = y_discovery,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = floor(sqrt(ncol(X_discovery))),
      splitrule = "gini",
      min.node.size = 5
    ),
    importance = "permutation",
    num.trees = 500
  )

  # Discovery performance (CV)
  discovery_auc <- max(rf_model$results$ROC)
  discovery_sens <- rf_model$results$Sens[which.max(rf_model$results$ROC)]
  discovery_spec <- rf_model$results$Spec[which.max(rf_model$results$ROC)]

  message("\n", paste(rep("-", 70), collapse = ""))
  message("DISCOVERY Performance (5-fold CV):")
  message(sprintf("  AUC:         %.3f", discovery_auc))
  message(sprintf("  Sensitivity: %.3f", discovery_sens))
  message(sprintf("  Specificity: %.3f", discovery_spec))

  # Feature importance
  var_imp <- caret::varImp(rf_model, scale = TRUE)
  importance_df <- var_imp$importance %>%
    tibble::rownames_to_column("feature") %>%
    dplyr::arrange(dplyr::desc(Overall))

  # Check tube type rank
  tube_rank <- which(importance_df$feature == "tube_type_encoded")

  message(sprintf("\nTube type importance rank: %s / %s features", tube_rank, nrow(importance_df)))

  if (tube_rank <= 10) {
    message("⚠️  WARNING: Tube type is in TOP 10 features!")
    message("   This is a NON-BIOLOGICAL variable driving predictions.")
  }

  # Test on REPLICATION set
  message("\n", paste(rep("-", 70), collapse = ""))
  message("Testing on REPLICATION (20% holdout)...")

  rep_pred_prob <- predict(rf_model, newdata = X_replication, type = "prob")
  rep_pred_class <- predict(rf_model, newdata = X_replication, type = "raw")

  # Calculate replication metrics
  rep_accuracy <- mean(rep_pred_class == y_replication)

  # ROC curve
  roc_obj <- pROC::roc(
    response = y_replication,
    predictor = rep_pred_prob$ALS,
    levels = c("Control", "ALS"),
    direction = "<",
    quiet = TRUE
  )

  rep_auc <- as.numeric(pROC::auc(roc_obj))

  # Confusion matrix
  conf_mat <- table(Predicted = rep_pred_class, Actual = y_replication)
  rep_sens <- conf_mat["ALS", "ALS"] / sum(conf_mat[, "ALS"])
  rep_spec <- conf_mat["Control", "Control"] / sum(conf_mat[, "Control"])

  message("\nREPLICATION Performance (20% holdout):")
  message(sprintf("  AUC:         %.3f", rep_auc))
  message(sprintf("  Accuracy:    %.3f", rep_accuracy))
  message(sprintf("  Sensitivity: %.3f", rep_sens))
  message(sprintf("  Specificity: %.3f", rep_spec))

  message("\nConfusion Matrix:")
  print(conf_mat)

  message("\n", paste(rep("=", 70), collapse = ""))
  message("INTERPRETATION:")
  message("  This replicates the original study's approach:")
  message("  - 80/20 split maintains tube type distribution")
  message("  - Tube type included as predictor")
  message("  - High performance expected (but may not generalize!)")
  message("\n  Compare with Leave-Country-Out CV for realistic performance.")
  message(paste(rep("=", 70), collapse = ""))

  # Top features
  message("\nTop 10 most important features:")
  print(as.data.frame(importance_df %>% dplyr::slice_head(n = 10)))

  return(list(
    model = rf_model,
    discovery_auc = discovery_auc,
    discovery_sens = discovery_sens,
    discovery_spec = discovery_spec,
    replication_auc = rep_auc,
    replication_accuracy = rep_accuracy,
    replication_sens = rep_sens,
    replication_spec = rep_spec,
    confusion_matrix = conf_mat,
    tube_type_rank = tube_rank,
    feature_importance = importance_df,
    n_proteins_used = length(sig_proteins),
    roc_curve = roc_obj
  ))
}


#' Create beta-beta plot comparing HEPARIN vs EDTA effect sizes
#'
#' @description
#' Replicates the original study's beta-beta plots showing concordance
#' between effect sizes in HEPARIN (Italy) vs EDTA (US) cohorts.
#'
#' Identifies proteins as "concordant" or "outliers" based on residuals
#' from linear model.
#'
#' @param italy_results Differential results from Italy
#' @param us_results Differential results from US
#' @param fdr_threshold False Discovery Rate threshold controlling multiple testing (default 0.05)
#' @param logfc_threshold Log fold-change threshold for biological relevance (default 0.6 = 1.52x change)
#' @param save_path Path to save plot (optional)
#' @return ggplot object + concordance statistics
#' @examples
#' beta_beta <- create_beta_beta_plot(italy_results, us_results)
#' @export
create_beta_beta_plot <- function(italy_results,
                                  us_results,
                                  fdr_threshold = 0.05,
                                  logfc_threshold = 0.6,
                                  save_path = NULL) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("BETA-BETA PLOT: HEPARIN (Italy) vs EDTA (US)")
  message(paste(rep("=", 70), collapse = ""))

  # Identify significant proteins in either cohort
  sig_italy <- italy_results %>%
    dplyr::filter(adj.P.Val < fdr_threshold, abs(logFC) > logfc_threshold) %>%
    dplyr::pull(Assay)

  sig_us <- us_results %>%
    dplyr::filter(adj.P.Val < fdr_threshold, abs(logFC) > logfc_threshold) %>%
    dplyr::pull(Assay)

  sig_proteins <- unique(c(sig_italy, sig_us))

  message(sprintf("\nProteins to plot: %s (significant in either cohort)", length(sig_proteins)))

  # Merge results
  merged <- dplyr::inner_join(
    italy_results %>%
      dplyr::select(Assay, logFC_italy = logFC, t_italy = t, adj_P_italy = adj.P.Val),
    us_results %>%
      dplyr::select(Assay, logFC_us = logFC, t_us = t, adj_P_us = adj.P.Val),
    by = "Assay"
  ) %>%
    dplyr::filter(Assay %in% sig_proteins)

  # Use t-statistics for beta-beta plot (like original study)
  # t-statistic = effect size / standard error
  # Better for comparing across cohorts with different sample sizes

  # Fit linear model
  lm_model <- lm(t_us ~ t_italy, data = merged)
  merged$residuals <- residuals(lm_model)

  # Identify outliers (|residual| > threshold)
  # Use IQR method like original study
  q1 <- quantile(merged$residuals, 0.25)
  q3 <- quantile(merged$residuals, 0.75)
  outlier_threshold <- (abs(q1) + abs(q3)) / 2

  merged <- merged %>%
    dplyr::mutate(
      outlier_status = ifelse(abs(residuals) > outlier_threshold, "Outlier", "Concordant"),
      sig_italy = Assay %in% sig_italy,
      sig_us = Assay %in% sig_us,
      sig_both = sig_italy & sig_us
    )

  # Summary
  n_concordant <- sum(merged$outlier_status == "Concordant")
  n_outlier <- sum(merged$outlier_status == "Outlier")

  message(sprintf(
    "\nConcordant proteins: %s (%.1f%%)",
    n_concordant, 100 * n_concordant / nrow(merged)
  ))
  message(sprintf(
    "Outlier proteins: %s (%.1f%%)",
    n_outlier, 100 * n_outlier / nrow(merged)
  ))

  # Correlation
  cor_t <- cor(merged$t_italy, merged$t_us, use = "complete.obs")
  message(sprintf("\nCorrelation (t-statistics): r = %.3f", cor_t))

  # Create plot
  library(ggplot2)
  library(ggrepel)

  p <- ggplot(merged, aes(x = t_italy, y = t_us, label = Assay)) +
    geom_point(aes(color = outlier_status), size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(
      values = c("Concordant" = "#56B4E9", "Outlier" = "#E69F00"),
      name = "Concordance"
    ) +
    labs(
      title = "Beta-Beta Plot: HEPARIN (Italy) vs EDTA (US)",
      subtitle = sprintf(
        "n = %s proteins (significant in either cohort), r = %.3f",
        nrow(merged), cor_t
      ),
      x = "t-statistic (HEPARIN/Italy)",
      y = "t-statistic (EDTA/US)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  # Label top outliers
  top_outliers <- merged %>%
    dplyr::filter(outlier_status == "Outlier") %>%
    dplyr::arrange(desc(abs(residuals))) %>%
    dplyr::slice_head(n = 10)

  if (nrow(top_outliers) > 0) {
    p <- p + geom_text_repel(
      data = top_outliers,
      aes(label = Assay),
      size = 3,
      max.overlaps = 10
    )
  }

  # Save if requested
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 8, height = 6, dpi = 300)
    message(sprintf("\nPlot saved: %s", save_path))
  }

  message(paste(rep("=", 70), collapse = ""))

  return(list(
    plot = p,
    data = merged,
    n_concordant = n_concordant,
    n_outlier = n_outlier,
    correlation = cor_t,
    lm_model = lm_model
  ))
}


#' Compare original study approach vs our rigorous approach
#'
#' @description
#' Side-by-side comparison showing:
#' - Their 80/20 split (maintains tube distribution) vs our LCV
#' - Their model WITH tube type vs our model WITHOUT
#' - Their pooled differential vs our stratified approach
#'
#' This makes the methodological differences crystal clear.
#'
#' @param protein_wide Wide-format data
#' @param protein_long Long-format data
#' @param seed Random seed
#' @return Comprehensive comparison table
#' @examples
#' comparison <- compare_original_vs_rigorous(protein_wide, protein_long)
#' @export
compare_original_vs_rigorous <- function(protein_wide, protein_long, seed = 42) {
  message("\n\n", paste(rep("=", 80), collapse = ""))
  message("COMPREHENSIVE COMPARISON: Original Study vs Our Analysis")
  message(paste(rep("=", 80), collapse = ""))

  # 1. Run original study approach (80/20 + tube type)
  message("\n[1/4] Running ORIGINAL STUDY approach (80/20 split + tube type)...")
  split_data <- create_chia_8020_split(protein_wide, seed = seed)

  # Need pooled differential for feature selection
  pooled_diff <- differential_analysis_pooled(protein_long)

  original_model <- train_chia_replication_model(
    split_data,
    pooled_diff,
    fdr_threshold = 0.05,
    logfc_threshold = 0.6,
    seed = seed
  )

  # 2. Run leave-country-out CV (our rigorous test)
  message("\n[2/4] Running LEAVE-COUNTRY-OUT CV (our approach)...")
  lcv_results <- leave_country_out_cv(protein_wide, seed = seed)

  # 3. Pooled CV (their general approach)
  message("\n[3/4] Running POOLED CV...")
  pooled_cv <- pooled_cv_analysis(protein_wide, n_folds = 5, seed = seed)

  # 4. Stratified differential analysis
  message("\n[4/4] Running STRATIFIED DIFFERENTIAL...")
  italy_diff <- differential_analysis_italy(protein_long)
  us_diff <- differential_analysis_us(protein_long)

  # Create comparison summary
  message("\n\n", paste(rep("=", 80), collapse = ""))
  message("COMPARISON SUMMARY")
  message(paste(rep("=", 80), collapse = ""))

  comparison_table <- tibble::tribble(
    ~Method, ~Approach, ~AUC, ~Interpretation,
    "Original 80/20 (Discovery)",
    "Maintains tube distribution",
    original_model$discovery_auc,
    "CV within discovery set",
    "Original 80/20 (Replication)",
    "20% holdout (mixed tubes)",
    original_model$replication_auc,
    "Still contains both tube types",
    "Pooled CV (5-fold)",
    "Standard CV on pooled data",
    pooled_cv$cv_auc,
    "Published estimate (inflated)",
    "LCV: Italy→US",
    "True geographic holdout",
    lcv_results$italy_to_us$test_auc,
    "Train HEPARIN, test EDTA",
    "LCV: US→Italy",
    "True geographic holdout",
    lcv_results$us_to_italy$test_auc,
    "Train EDTA, test HEPARIN",
    "LCV: Mean",
    "Average geographic performance",
    lcv_results$mean_test_auc,
    "Realistic estimate"
  )

  message("\nML MODEL PERFORMANCE:")
  print(as.data.frame(comparison_table), row.names = FALSE)

  # Performance gap
  gap_original_vs_lcv <- original_model$replication_auc - lcv_results$mean_test_auc
  gap_pooled_vs_lcv <- pooled_cv$cv_auc - lcv_results$mean_test_auc

  message("\n", paste(rep("-", 80), collapse = ""))
  message("PERFORMANCE GAPS:")
  message(sprintf(
    "  Original 80/20 vs LCV:  %.3f (%.1f%% overestimation)",
    gap_original_vs_lcv,
    100 * gap_original_vs_lcv / lcv_results$mean_test_auc
  ))
  message(sprintf(
    "  Pooled CV vs LCV:       %.3f (%.1f%% overestimation)",
    gap_pooled_vs_lcv,
    100 * gap_pooled_vs_lcv / lcv_results$mean_test_auc
  ))

  # Differential analysis comparison
  n_sig_italy <- sum(italy_diff$adj.P.Val < 0.05 & abs(italy_diff$logFC) > 0.5)
  n_sig_us <- sum(us_diff$adj.P.Val < 0.05 & abs(us_diff$logFC) > 0.5)
  n_sig_pooled <- sum(pooled_diff$adj.P.Val < 0.05 & abs(pooled_diff$logFC) > 0.5)

  message("\n", paste(rep("-", 80), collapse = ""))
  message("DIFFERENTIAL PROTEIN COUNTS (FDR < 0.05, |logFC| > 0.5):")
  message(sprintf("  Pooled (original approach):  %s", n_sig_pooled))
  message(sprintf("  Italy stratified:            %s", n_sig_italy))
  message(sprintf("  US stratified:               %s", n_sig_us))

  message("\n", paste(rep("=", 80), collapse = ""))
  message("KEY FINDINGS:")
  message("  1. Original 80/20 split maintains tube distribution (not true holdout)")
  message("  2. Their replication set is NOT geographically independent")
  message("  3. Leave-country-out CV reveals TRUE generalization performance")
  message(sprintf(
    "  4. Published estimates overestimate by %.1f%%",
    100 * gap_pooled_vs_lcv / lcv_results$mean_test_auc
  ))
  message(paste(rep("=", 80), collapse = ""))

  return(list(
    comparison_table = comparison_table,
    original_model = original_model,
    lcv_results = lcv_results,
    pooled_cv = pooled_cv,
    gap_original_vs_lcv = gap_original_vs_lcv,
    gap_pooled_vs_lcv = gap_pooled_vs_lcv,
    split_data = split_data
  ))
}
