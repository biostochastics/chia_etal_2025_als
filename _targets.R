# _targets.R
# ALS Biomarker Confounding Bias Investigation Pipeline
#
# This targets pipeline implements Week 1 analysis:
# 1. Data loading and country label reconstruction
# 2. Confounding bias quantification (Cramér's V)
# 3. REVERSE PREDICTION TEST (predicting tube type from proteins)
#
# Run with: targets::tar_make()
# Visualize with: targets::tar_visnetwork()

# Load packages
library(targets)
library(tarchetypes)

# Set pipeline options
tar_option_set(
  packages = c(
    # Data manipulation
    "tidyverse",
    "data.table",
    "dtplyr",
    "here",

    # Statistical analysis
    "broom",
    "DescTools",
    "vcd",
    "limma",

    # Machine learning
    "caret",
    "ranger",
    "glmnet",
    "DALEX",
    "pROC",

    # Visualization
    "ggplot2",
    "patchwork",
    "ggpubr",
    "ggrepel",
    "ggforce",
    "corrplot"
  ),
  format = "rds", # Default format for targets
  error = "continue" # Continue pipeline even if one target fails
)

# Preserve source code for computational transparency appendix
options(keep.source = TRUE)

# Source all R functions
tar_source("R/")

# Define pipeline
list(
  # ============================================================================
  # DATA LOADING
  # ============================================================================

  # Track the data file
  tar_target(
    name = raw_data_file,
    command = "original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt",
    format = "file"
  ),

  # Load raw OLINK data
  tar_target(
    name = raw_data,
    command = load_raw_olink_data(raw_data_file)
  ),

  # Add country labels based on plate + tube type
  tar_target(
    name = data_with_country,
    command = add_country_labels(raw_data)
  ),

  # Validate data integrity
  tar_target(
    name = data_validation,
    command = validate_data_integrity(data_with_country)
  ),

  # Create sample-level metadata
  tar_target(
    name = sample_metadata,
    command = create_sample_metadata(data_with_country)
  ),

  # ============================================================================
  # CONFOUNDING BIAS QUANTIFICATION
  # ============================================================================

  # Calculate pairwise Cramér's V values
  tar_target(
    name = cramers_v_pairwise,
    command = calculate_pairwise_cramers_v(sample_metadata)
  ),

  # Create contingency tables
  tar_target(
    name = contingency_tables,
    command = create_contingency_tables(sample_metadata)
  ),

  # Perform chi-square tests
  tar_target(
    name = chi_square_tests,
    command = perform_chi_square_tests(sample_metadata)
  ),

  # Stratified confounding analysis
  tar_target(
    name = stratified_analysis,
    command = analyze_stratified_confounding(sample_metadata)
  ),

  # Comprehensive bias quantification report
  tar_target(
    name = bias_report,
    command = quantify_confounding_bias(sample_metadata)
  ),

  # ============================================================================
  # MACHINE LEARNING EVALUATION
  # ============================================================================

  # Pivot protein data to wide format
  tar_target(
    name = protein_wide,
    command = pivot_protein_wide(data_with_country)
  ),

  # **CRITICAL: Reverse Prediction Test**
  # Can we predict tube type from protein expression?
  # If AUC > 0.9 → tube type effects dominate biological signal
  tar_target(
    name = reverse_prediction_results,
    command = reverse_prediction_test(
      protein_wide,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Train model WITH tube type as feature
  tar_target(
    name = model_with_tube,
    command = train_model_with_tube_type(
      protein_wide,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Train model WITHOUT tube type
  tar_target(
    name = model_without_tube,
    command = train_model_without_tube_type(
      protein_wide,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Compare models to quantify tube type contribution
  tar_target(
    name = tube_type_contribution,
    command = compare_tube_type_contribution(
      protein_wide,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # ============================================================================
  # LEAVE-COUNTRY-OUT CROSS-VALIDATION (CRITICAL TEST)
  # ============================================================================

  # Leave-country-out CV: Train Italy → Test US, Train US → Test Italy
  tar_target(
    name = lcv_results,
    command = leave_country_out_cv(
      protein_wide,
      seed = 25667777
    )
  ),

  # Pooled CV (original study approach)
  tar_target(
    name = pooled_cv_results,
    command = pooled_cv_analysis(
      protein_wide,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # **PRIMARY EVIDENCE: Compare pooled vs leave-country-out**
  tar_target(
    name = pooled_vs_lcv_comparison,
    command = compare_pooled_vs_lcv(
      protein_wide,
      seed = 25667777
    )
  ),

  # Within-country CV (Italy and US separately)
  # Tests performance in homogeneous anticoagulant contexts
  tar_target(
    name = within_country_cv_results,
    command = within_country_cross_validation(
      protein_wide,
      n_folds = 5,
      seed = 42
    )
  ),

  # ============================================================================
  # DIFFERENTIAL PROTEIN EXPRESSION ANALYSIS
  # ============================================================================

  # Differential analysis: Italy only (HEPARIN)
  tar_target(
    name = differential_italy,
    command = differential_analysis_italy(data_with_country)
  ),

  # Differential analysis: US only (EDTA)
  tar_target(
    name = differential_us,
    command = differential_analysis_us(data_with_country)
  ),

  # Differential analysis: Pooled (confounded)
  tar_target(
    name = differential_pooled,
    command = differential_analysis_pooled(data_with_country)
  ),

  # Meta-analysis: Test concordance between Italy and US
  tar_target(
    name = protein_concordance,
    command = test_protein_concordance(
      differential_italy,
      differential_us
    )
  ),

  # **PRIMARY EVIDENCE: Compare stratified vs pooled differential analysis**
  tar_target(
    name = stratified_vs_pooled_comparison,
    command = compare_stratified_vs_pooled(
      differential_italy,
      differential_us,
      differential_pooled
    )
  ),

  # ============================================================================
  # TUBE TYPE EFFECTS WITH COVARIATE ADJUSTMENT (LIMMA)
  # ============================================================================
  # Test for tube type (EDTA vs HEPARIN) effects within diagnosis groups
  # adjusting for age and sex using limma

  # Tube effects in ALS patients
  tar_target(
    name = tube_effects_als,
    command = calculate_tube_effects_limma(
      data_with_country,
      diagnosis_group = "ALS"
    )
  ),

  # Tube effects in Healthy controls
  tar_target(
    name = tube_effects_healthy,
    command = calculate_tube_effects_limma(
      data_with_country,
      diagnosis_group = "Healthy_control"
    )
  ),

  # ============================================================================
  # RESIDUAL CONFOUNDING IN TUBE-ROBUST PROTEINS (NEW ANALYSIS)
  # ============================================================================
  # Even proteins that are significant in BOTH Italy and US (with concordant
  # direction) may still be confounded with tube type. This section quantifies
  # the extent of residual confounding in "tube-robust" proteins.

  # Test 1: Can tube-robust proteins predict tube type?
  # (Reverse prediction using ONLY the 22 tube-robust proteins)
  tar_target(
    name = robust_proteins_reverse_prediction,
    command = reverse_prediction_test_robust_only(
      protein_wide,
      protein_concordance,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Test 2: Do tube-robust proteins differ by tube type in healthy controls?
  # (If yes, the difference is technical artifact, not disease biology)
  tar_target(
    name = healthy_tube_effects,
    command = tube_effect_in_healthy(
      protein_wide,
      protein_concordance,
      alpha = 0.05
    )
  ),

  # Test 3: Effect size decomposition - compare disease vs tube effects
  # (Calculate ratio: |tube effect| / |disease effect| for each protein)
  tar_target(
    name = effect_decomposition,
    command = effect_size_decomposition(
      protein_wide,
      protein_concordance,
      alpha = 0.05
    )
  ),

  # ============================================================================
  # VISUALIZATIONS
  # ============================================================================

  # ROC curve for reverse prediction test
  tar_target(
    name = plot_roc,
    command = plot_reverse_prediction_roc(
      reverse_prediction_results,
      save_path = "outputs/figures/reverse_prediction_roc.png"
    )
  ),

  # Tube type feature importance
  tar_target(
    name = plot_features,
    command = plot_tube_type_features(
      reverse_prediction_results,
      save_path = "outputs/figures/tube_type_features.png"
    )
  ),

  # Confounding structure
  tar_target(
    name = plot_confounding,
    command = plot_confounding_structure(
      sample_metadata,
      save_path = "outputs/figures/confounding_structure.png"
    )
  ),

  # PCA analysis
  tar_target(
    name = plot_pca,
    command = plot_pca_tube_vs_diagnosis(
      protein_wide,
      save_path = "outputs/figures/pca_tube_vs_diagnosis.png"
    )
  ),

  # Combined summary figure
  tar_target(
    name = summary_figure,
    command = create_summary_figure(
      reverse_prediction_results,
      sample_metadata,
      protein_wide,
      save_path = "outputs/figures/summary_panel.png"
    )
  ),

  # ============================================================================
  # EXACT REPLICATION OF ORIGINAL STUDY
  # ============================================================================

  # Create 80/20 split (their approach: maintains tube distribution)
  tar_target(
    name = chia_8020_split,
    command = create_chia_8020_split(
      protein_wide,
      split_ratio = 0.8,
      seed = 25667777
    )
  ),

  # Train model using their exact approach
  # (80/20 split + tube type as feature + diff proteins)
  tar_target(
    name = chia_replication_model,
    command = train_chia_replication_model(
      chia_8020_split,
      differential_pooled,
      fdr_threshold = 0.05,
      logfc_threshold = 0.6,
      seed = 25667777
    )
  ),

  # Beta-beta plot (HEPARIN vs EDTA effect concordance)
  # NOTE: Commented out - not used in final report, uses hardcoded colors
  # tar_target(
  #   name = beta_beta_plot,
  #   command = create_beta_beta_plot(
  #     differential_italy,
  #     differential_us,
  #     fdr_threshold = 0.05,
  #     logfc_threshold = 0.6,
  #     save_path = "outputs/figures/beta_beta_plot.png"
  #   )
  # ),

  # **CRITICAL COMPARISON: Their approach vs ours**
  tar_target(
    name = original_vs_rigorous_comparison,
    command = compare_original_vs_rigorous(
      protein_wide,
      data_with_country,
      seed = 25667777
    )
  ),

  # ============================================================================
  # PROTEIN OVERLAP ANALYSIS (CRITICAL COMPARISON)
  # ============================================================================

  # Get original biomarkers from Chia et al. (2025)
  tar_target(
    name = original_biomarkers,
    command = get_original_biomarkers()
  ),

  # Get all 33 differential proteins from original study
  tar_target(
    name = original_differential_proteins,
    command = get_original_differential_proteins()
  ),

  # Get tube-affected proteins from our reverse prediction
  tar_target(
    name = tube_affected_proteins,
    command = get_tube_affected_proteins(
      reverse_prediction_results,
      top_n = 20
    )
  ),

  # Calculate overlap statistics
  tar_target(
    name = overlap_stats,
    command = calculate_overlap_statistics(
      original_biomarkers,
      tube_affected_proteins
    )
  ),

  # Venn diagram showing overlap
  # NOTE: Commented out - not used in final report, uses hardcoded colors
  # tar_target(
  #   name = plot_overlap_venn,
  #   command = plot_overlap_venn(
  #     overlap_stats,
  #     save_path = "outputs/figures/overlap_venn.png"
  #   )
  # ),

  # Importance heatmap for overlapping proteins
  tar_target(
    name = plot_overlap_heatmap,
    command = plot_overlap_importance_heatmap(
      reverse_prediction_results,
      overlap_stats,
      save_path = "outputs/figures/overlap_importance_heatmap.png"
    )
  ),

  # Concordance scatter plot
  tar_target(
    name = plot_overlap_concordance,
    command = plot_overlap_concordance(
      reverse_prediction_results,
      overlap_stats,
      save_path = "outputs/figures/overlap_concordance.png"
    )
  ),

  # Comprehensive overlap report
  # NOTE: Commented out - not used in final report
  # tar_target(
  #   name = overlap_report,
  #   command = generate_overlap_report(
  #     overlap_stats,
  #     reverse_prediction_results,
  #     save_path = "outputs/reports/protein_overlap_analysis.md"
  #   )
  # ),

  # ============================================================================
  # BIOMARKER COLLECTIVE TUBE-TYPE PREDICTION (CRITICAL TEST)
  # ============================================================================

  # Extract proteins significant in both Italy and US (concordant direction)
  tar_target(
    name = concordant_biomarkers,
    command = protein_concordance %>%
      dplyr::filter(sig_both, same_direction) %>%
      dplyr::pull(Assay)
  ),

  # Test if original 17-protein biomarker panel predicts tube type
  tar_target(
    name = biomarker_tube_prediction,
    command = test_biomarker_tube_prediction(
      protein_wide,
      original_biomarkers,
      reverse_prediction_results,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Test if concordant biomarkers predict tube type
  tar_target(
    name = concordant_tube_prediction,
    command = test_biomarker_tube_prediction(
      protein_wide,
      concordant_biomarkers,
      reverse_prediction_results,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Test top N proteins for tube type prediction
  tar_target(
    name = top_proteins_tube_prediction,
    command = test_top_proteins_tube_prediction(
      protein_wide,
      reverse_prediction_results,
      n_top = 10,
      n_folds = 5,
      seed = 25667777
    )
  ),

  # Visualization comparing biomarker-only vs full model
  # NOTE: Commented out - not used in final report, uses hardcoded colors
  # tar_target(
  #   name = plot_biomarker_tube_comparison,
  #   command = plot_biomarker_tube_comparison(
  #     biomarker_tube_prediction,
  #     reverse_prediction_results,
  #     save_path = "outputs/figures/biomarker_tube_prediction.png"
  #   )
  # ),

  # Summary report
  # NOTE: Commented out - not used in final report
  # tar_target(
  #   name = biomarker_tube_summary,
  #   command = summarize_biomarker_tube_prediction(
  #     biomarker_tube_prediction,
  #     save_path = "outputs/reports/biomarker_tube_prediction.md"
  #   )
  # ),

  # ============================================================================
  # VISUALIZATIONS
  # ============================================================================

  # 1. ROC curves for leave-country-out CV
  tar_target(
    name = viz_roc_curves,
    command = plot_lcv_roc_curves(
      lcv_results,
      pooled_cv_results,
      save_path = "outputs/figures/01_lcv_roc_curves.png"
    )
  ),

  # 2. Performance comparison bar chart
  tar_target(
    name = viz_performance,
    command = plot_performance_comparison(
      pooled_vs_lcv_comparison,
      save_path = "outputs/figures/02_performance_comparison.png"
    )
  ),

  # 3. Confusion matrices
  tar_target(
    name = viz_confusion,
    command = plot_confusion_matrices(
      lcv_results,
      save_path = "outputs/figures/03_confusion_matrices.png"
    )
  ),

  # 4. Protein counts comparison
  tar_target(
    name = viz_protein_counts,
    command = plot_protein_counts(
      stratified_vs_pooled_comparison,
      save_path = "outputs/figures/04_protein_counts.png"
    )
  ),

  # 5. Protein overlap Venn diagram
  tar_target(
    name = viz_venn,
    command = plot_protein_overlap_venn(
      stratified_vs_pooled_comparison,
      save_path = "outputs/figures/05_protein_venn.png"
    )
  ),

  # 6. Tube-robust protein forest plot
  tar_target(
    name = viz_forest,
    command = plot_tube_robust_forest(
      protein_concordance,
      top_n = 15,
      save_path = "outputs/figures/06_forest_plot.png"
    )
  ),

  # 7. Effect size correlation
  tar_target(
    name = viz_correlation,
    command = plot_effect_size_correlation(
      protein_concordance,
      save_path = "outputs/figures/07_effect_correlation.png"
    )
  ),

  # 7b. Effect size correlation by protein set (All vs Pooled-significant)
  tar_target(
    name = viz_correlation_by_set,
    command = plot_effect_correlation_by_set(
      protein_concordance,
      differential_pooled,
      save_path = "outputs/figures/07b_effect_correlation_by_set.png"
    )
  ),

  # 8. Comprehensive summary panel
  tar_target(
    name = viz_summary_panel,
    command = create_summary_panel(
      lcv_results,
      pooled_cv_results,
      pooled_vs_lcv_comparison,
      stratified_vs_pooled_comparison,
      protein_concordance,
      save_path = "outputs/figures/08_summary_panel.png"
    )
  ),

  # Export all visualizations at once (master function)
  tar_target(
    name = all_visualizations,
    command = export_all_visualizations(
      lcv_results,
      pooled_cv_results,
      pooled_vs_lcv_comparison,
      stratified_vs_pooled_comparison,
      protein_concordance,
      output_dir = "outputs/figures"
    )
  ),

  # ==== NEW: Residual Confounding Visualizations ====

  # 9. Reverse prediction comparison (all vs robust-only proteins)
  tar_target(
    name = viz_robust_reverse_comparison,
    command = plot_robust_reverse_prediction_comparison(
      robust_proteins_reverse_prediction,
      reverse_prediction_results,
      save_path = "outputs/figures/09_robust_reverse_comparison.png"
    )
  ),

  # 10. Tube effects in healthy controls
  tar_target(
    name = viz_healthy_tube_effects,
    command = plot_healthy_tube_effects(
      healthy_tube_effects,
      top_n = 15,
      save_path = "outputs/figures/10_healthy_tube_effects.png"
    )
  ),

  # 11. Effect size decomposition ratios
  tar_target(
    name = viz_effect_decomposition,
    command = plot_effect_decomposition_ratios(
      effect_decomposition,
      save_path = "outputs/figures/11_effect_decomposition.png"
    )
  ),

  # 12a. Tube effect size distribution in healthy controls (full range, labeled outliers)
  tar_target(
    name = viz_tube_effect_distribution_full,
    command = plot_tube_effect_size_distribution(
      tube_effects_healthy,
      save_path = "outputs/figures/12a_tube_effect_distribution_full.png",
      xlim_range = NULL,
      label_outliers = FALSE
    )
  ),

  # 12b. Tube effect size distribution in healthy controls (zoomed to -1000 to 1000)
  tar_target(
    name = viz_tube_effect_distribution_zoomed,
    command = plot_tube_effect_size_distribution(
      tube_effects_healthy,
      save_path = "outputs/figures/12b_tube_effect_distribution_zoomed.png",
      xlim_range = c(-1000, 1000),
      label_outliers = TRUE
    )
  ),

  # ============================================================================
  # SUMMARY REPORT
  # ============================================================================

  # Generate comprehensive summary of ALL findings
  tar_target(
    name = investigation_summary,
    command = {
      list(
        # ============ DATASET OVERVIEW ============
        n_samples = nrow(sample_metadata),
        n_proteins = ncol(protein_wide) - 6,
        diagnosis_counts = table(sample_metadata$Diagnosis),
        tube_type_counts = table(sample_metadata$Plasma_collection_tube_type),
        country_counts = table(sample_metadata$country),

        # ============ CONFOUNDING STRUCTURE ============
        neuro_controls_all_edta = all(
          sample_metadata$Plasma_collection_tube_type[
            sample_metadata$Diagnosis == "Neurological_control"
          ] == "EDTA"
        ),
        als_pct_heparin = 100 * sum(
          sample_metadata$Diagnosis == "ALS" &
            sample_metadata$Plasma_collection_tube_type == "HEPARIN"
        ) / sum(sample_metadata$Diagnosis == "ALS"),

        # ============ REVERSE PREDICTION TEST ============
        reverse_prediction_auc = reverse_prediction_results$auc,
        reverse_prediction_sens = reverse_prediction_results$sensitivity,
        reverse_prediction_spec = reverse_prediction_results$specificity,

        # ============ MODEL COMPARISON (WITH/WITHOUT TUBE) ============
        auc_with_tube = model_with_tube$auc,
        auc_without_tube = model_without_tube$auc,
        tube_type_rank = model_with_tube$tube_type_rank,
        tube_contribution_auc_diff = tube_type_contribution$auc_difference,
        tube_contribution_pct = tube_type_contribution$percent_change,

        # ============ LEAVE-COUNTRY-OUT CV ============
        lcv_italy_to_us_auc = lcv_results$italy_to_us$test_auc,
        lcv_us_to_italy_auc = lcv_results$us_to_italy$test_auc,
        lcv_mean_test_auc = lcv_results$mean_test_auc,
        lcv_performance_drop = lcv_results$performance_drop,

        # ============ POOLED VS LCV COMPARISON ============
        pooled_cv_auc = pooled_vs_lcv_comparison$pooled_auc,
        lcv_auc = pooled_vs_lcv_comparison$lcv_mean_auc,
        pooled_vs_lcv_gap = pooled_vs_lcv_comparison$auc_gap,
        pooled_vs_lcv_pct_drop = pooled_vs_lcv_comparison$percent_drop,

        # ============ DIFFERENTIAL ANALYSIS ============
        # Significant proteins (FDR < 0.05, |logFC| > 0.5)
        n_sig_italy = sum(differential_italy$adj.P.Val < 0.05 & abs(differential_italy$logFC) > 0.5),
        n_sig_us = sum(differential_us$adj.P.Val < 0.05 & abs(differential_us$logFC) > 0.5),
        n_sig_pooled = sum(differential_pooled$adj.P.Val < 0.05 & abs(differential_pooled$logFC) > 0.5),

        # Concordance
        n_sig_both_countries = length(stratified_vs_pooled_comparison$sig_both_strata),
        n_pooled_only = length(stratified_vs_pooled_comparison$pooled_only),
        pct_pooled_replicate_both = stratified_vs_pooled_comparison$pct_replicate_in_both,

        # ============ EXACT REPLICATION OF ORIGINAL STUDY ============
        # Their 80/20 split + tube type approach
        chia_discovery_auc = chia_replication_model$discovery_auc,
        chia_replication_auc = chia_replication_model$replication_auc,
        chia_tube_rank = chia_replication_model$tube_type_rank,

        # Beta-beta plot concordance (commented out - not generated)
        # beta_beta_correlation = beta_beta_plot$correlation,
        # beta_beta_n_concordant = beta_beta_plot$n_concordant,
        # beta_beta_n_outlier = beta_beta_plot$n_outlier,

        # Performance gaps (original vs rigorous)
        gap_chia_vs_lcv = original_vs_rigorous_comparison$gap_original_vs_lcv,
        gap_pooled_vs_lcv = original_vs_rigorous_comparison$gap_pooled_vs_lcv,

        # ============ INTERPRETATION ============
        analysis_date = Sys.time(),

        # Overall verdict
        overall_interpretation = paste(
          # Reverse prediction
          if (reverse_prediction_results$auc > 0.95) {
            "CRITICAL: Reverse prediction AUC > 0.95 - tube type effects dominate."
          } else if (reverse_prediction_results$auc > 0.9) {
            "SEVERE: Reverse prediction AUC > 0.9 - major tube confounding."
          } else {
            "MODERATE: Some tube type effects present."
          },
          "\n",
          # Leave-country-out
          if (pooled_vs_lcv_comparison$auc_gap > 0.2) {
            "CRITICAL: Pooled CV overestimates performance by >0.2 AUC. Model does NOT generalize."
          } else if (pooled_vs_lcv_comparison$auc_gap > 0.1) {
            "SUBSTANTIAL: Pooled CV overestimates by 0.1-0.2 AUC. Limited generalizability."
          } else {
            "ACCEPTABLE: Model generalizes reasonably across countries."
          },
          "\n",
          # Differential analysis
          if (stratified_vs_pooled_comparison$pct_replicate_in_both < 30) {
            "CRITICAL: <30% of pooled proteins replicate in both strata. Heavy confounding."
          } else if (stratified_vs_pooled_comparison$pct_replicate_in_both < 50) {
            "SUBSTANTIAL: 30-50% replication. Major confounding in differential analysis."
          } else {
            "MODERATE: Reasonable replication across strata."
          }
        )
      )
    }
  ),

  # ============================================================================
  # COMPUTATIONAL TRANSPARENCY APPENDIX
  # ============================================================================
  # NOTE: Appendix metadata (tar_meta(), tar_manifest(), tar_network()) is
  # called directly in the Quarto document appendix chunks, not as a target.
  # This avoids the limitation that tar_meta() cannot be called during pipeline.

  # ============================================================================
  # REPORT GENERATION
  # ============================================================================

  # Integrate main report into pipeline
  # Report automatically rebuilds when upstream targets change
  tarchetypes::tar_quarto(
    name = main_report,
    path = "reports/confounding_investigation.qmd",
    extra_files = c("reports/references.bib", "reports/nature.csl", "reports/custom-styles-dark.css"),
    quiet = FALSE
  )
)
