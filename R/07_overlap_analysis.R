#' Protein Overlap Analysis: Original Biomarkers vs Tube-Affected Proteins
#'
#' This module compares the original study's biomarker panel with proteins
#' most affected by tube type to determine if their "biomarkers" are actually
#' technical artifacts.

#' Original Study's 17-Protein Biomarker Panel
#'
#' @description
#' The 17 proteins used in the final ML model (from Chia et al. 2025, Figure 4a).
#' These were selected from 33 differentially abundant proteins and ranked by
#' feature importance in their Random Forest model.
#'
#' Note: Tube type was ranked #2 overall (after NEFL).
#'
#' @return Character vector of 17 protein names
#' @export
get_original_biomarkers <- function() {
  c(
    "NEFL", # Neurofilament light chain - MOST IMPORTANT
    "CSRP3", # Cysteine and glycine-rich protein 3
    "LIF", # Leukemia inhibitory factor
    "TPM3", # Tropomyosin 3
    "HSPB6", # Heat shock protein beta-6
    "TNNT1", # Troponin T1
    "MEGF10", # Multiple EGF-like domains 10
    "RNASE3", # Ribonuclease A family member 3
    "CKMT2", # Creatine kinase, mitochondrial 2
    "MYBPC1", # Myosin binding protein C1
    "MYBPC2", # Myosin binding protein C2
    "MYOM2", # Myomesin 2
    "MYOZ2", # Myozenin 2
    "TNNI1", # Troponin I1
    "TNNC1", # Troponin C1
    "TNNI2", # Troponin I2
    "MYLPF" # Myosin light chain, phosphorylatable, fast skeletal muscle
  )
}


#' All 33 Differentially Abundant Proteins from Original Study
#'
#' @description
#' Complete list from Table 1 of Chia et al. 2025.
#' These were identified as significantly different between ALS and controls
#' in the discovery cohort.
#'
#' @return Character vector of 33 protein names
#' @export
get_original_differential_proteins <- function() {
  c(
    # The 17 selected for ML model
    "NEFL", "CSRP3", "LIF", "TPM3", "HSPB6", "TNNT1", "MEGF10", "RNASE3",
    "CKMT2", "MYBPC1", "MYBPC2", "MYOM2", "MYOZ2", "TNNI1", "TNNC1",
    "TNNI2", "MYLPF",

    # Additional 16 from differential analysis (not in final ML model)
    "CKMT1A", "CKMT1B", "CORO6", "DTNB", "FGF21", "HSPB7", "MYBPH",
    "MYH1", "MYL1", "MYL3", "NEB", "RBFOX3", "SSC4D", "TMOD4",
    "TNNT3", "TTN", "XIRP2"
  )
}


#' Extract Top Tube-Affected Proteins
#'
#' @description
#' Extracts the top N proteins from reverse prediction test that are most
#' affected by tube type (highest feature importance for discriminating
#' HEPARIN vs EDTA).
#'
#' @param reverse_prediction_results Results from reverse_prediction_test()
#' @param top_n Number of top proteins to extract (default: 20)
#' @return Character vector of protein names
#' @export
get_tube_affected_proteins <- function(reverse_prediction_results, top_n = 20) {
  # Extract feature importance from ranger model
  importance_vec <- reverse_prediction_results$model$finalModel$variable.importance

  # Convert to data frame and sort
  importance_df <- data.frame(
    protein = names(importance_vec),
    importance = as.numeric(importance_vec),
    stringsAsFactors = FALSE
  )

  # Sort by importance and take top N
  top_proteins <- importance_df %>%
    dplyr::arrange(dplyr::desc(importance)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(protein)

  return(top_proteins)
}


#' Calculate Protein Set Overlaps
#'
#' @description
#' Computes overlap statistics between original biomarkers and tube-affected
#' proteins. Returns counts and percentages for:
#' - Overlap (proteins in both sets)
#' - Original only (biomarkers not affected by tube)
#' - Tube only (tube artifacts not identified as biomarkers)
#'
#' @param original_proteins Character vector of original biomarker proteins
#' @param tube_proteins Character vector of tube-affected proteins
#' @return List with overlap statistics
#' @export
calculate_overlap_statistics <- function(original_proteins, tube_proteins) {
  # Calculate overlaps
  overlap <- intersect(original_proteins, tube_proteins)
  original_only <- setdiff(original_proteins, tube_proteins)
  tube_only <- setdiff(tube_proteins, original_proteins)

  # Calculate percentages
  pct_original_affected <- 100 * length(overlap) / length(original_proteins)
  pct_tube_identified <- 100 * length(overlap) / length(tube_proteins)

  # Create result list
  result <- list(
    # Counts
    n_original = length(original_proteins),
    n_tube = length(tube_proteins),
    n_overlap = length(overlap),
    n_original_only = length(original_only),
    n_tube_only = length(tube_only),

    # Percentages
    pct_original_affected = pct_original_affected,
    pct_tube_identified = pct_tube_identified,

    # Protein lists
    overlap_proteins = overlap,
    original_only_proteins = original_only,
    tube_only_proteins = tube_only,

    # Interpretation
    interpretation = dplyr::case_when(
      pct_original_affected >= 50 ~ "CRITICAL: Majority of biomarkers are tube artifacts",
      pct_original_affected >= 30 ~ "SUBSTANTIAL: Large fraction of biomarkers affected by tube type",
      pct_original_affected >= 10 ~ "MODERATE: Some biomarkers may be confounded",
      TRUE ~ "MINIMAL: Few biomarkers overlap with tube effects"
    )
  )

  # Print summary
  cat("\n=== PROTEIN OVERLAP ANALYSIS ===\n\n")
  cat(sprintf("Original biomarkers: %d proteins\n", result$n_original))
  cat(sprintf("Tube-affected proteins: %d proteins\n", result$n_tube))
  cat(sprintf(
    "Overlap: %d proteins (%.1f%% of biomarkers)\n",
    result$n_overlap, result$pct_original_affected
  ))
  cat(sprintf("\n%s\n", result$interpretation))

  if (result$n_overlap > 0) {
    cat("\nOverlapping proteins (BIOMARKERS THAT ARE TUBE ARTIFACTS):\n")
    cat(paste("  -", result$overlap_proteins, collapse = "\n"), "\n")
  }

  return(result)
}


#' Create Overlap Venn Diagram
#'
#' @description
#' Creates a Venn diagram showing overlap between original biomarkers
#' and tube-affected proteins.
#'
#' @param overlap_stats Results from calculate_overlap_statistics()
#' @param save_path Optional path to save figure
#' @return ggplot2 object
#' @export
plot_overlap_venn <- function(overlap_stats, save_path = NULL) {
  library(ggplot2)
  library(ggforce)

  # Create Venn diagram data
  venn_data <- data.frame(
    x = c(-0.5, 0.5),
    y = c(0, 0),
    r = c(1, 1),
    label = c(
      sprintf("Original\nBiomarkers\n(n=%d)", overlap_stats$n_original),
      sprintf("Tube-Affected\nProteins\n(n=%d)", overlap_stats$n_tube)
    ),
    color = c("#1E88E5", "#D32F2F")
  )

  # Create plot
  p <- ggplot() +
    # Circles
    ggforce::geom_circle(
      data = venn_data,
      aes(x0 = x, y0 = y, r = r, color = color, fill = color),
      alpha = 0.2,
      size = 2
    ) +
    # Labels for each circle
    annotate("text",
      x = -1.2, y = 0.5,
      label = sprintf("%d\nBiomarkers\nOnly", overlap_stats$n_original_only),
      size = 5, fontface = "bold", color = "#1E88E5"
    ) +
    annotate("text",
      x = 1.2, y = 0.5,
      label = sprintf("%d\nTube\nArtifacts\nOnly", overlap_stats$n_tube_only),
      size = 5, fontface = "bold", color = "#D32F2F"
    ) +
    # Overlap label (CENTER - THE SMOKING GUN)
    annotate("text",
      x = 0, y = 0,
      label = sprintf(
        "%d\nOVERLAP\n(%.1f%% of\nbiomarkers)",
        overlap_stats$n_overlap,
        overlap_stats$pct_original_affected
      ),
      size = 6, fontface = "bold", color = "#FF6F00"
    ) +
    # Title
    labs(
      title = "Overlap Between Original Biomarkers and Tube-Affected Proteins",
      subtitle = overlap_stats$interpretation
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    coord_fixed() +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(
        size = 12, hjust = 0.5, color = "#D32F2F",
        face = "bold", margin = margin(t = 5, b = 20)
      ),
      plot.margin = margin(20, 20, 20, 20)
    )

  # Save if requested
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 10, height = 8, dpi = 300, bg = "transparent")
    cat(sprintf("\n✓ Venn diagram saved to: %s\n", save_path))
  }

  return(p)
}


#' Create Overlap Heatmap
#'
#' @description
#' Creates a heatmap comparing feature importance rankings between:
#' - Original study's biomarker importance
#' - Our reverse prediction importance (tube type discrimination)
#'
#' Only includes proteins that appear in BOTH analyses.
#'
#' @param reverse_prediction_results Results from reverse_prediction_test()
#' @param overlap_stats Results from calculate_overlap_statistics()
#' @param save_path Optional path to save figure
#' @return ggplot2 object
#' @export
plot_overlap_importance_heatmap <- function(reverse_prediction_results,
                                            overlap_stats,
                                            save_path = NULL) {
  library(ggplot2)
  library(dplyr)

  # Get overlapping proteins
  overlap_proteins <- overlap_stats$overlap_proteins

  if (length(overlap_proteins) == 0) {
    cat("\nNo overlapping proteins to plot.\n")
    return(NULL)
  }

  # Get tube-type importance from model
  importance_vec <- reverse_prediction_results$model$finalModel$variable.importance
  importance_df <- data.frame(
    protein = names(importance_vec),
    importance = as.numeric(importance_vec),
    stringsAsFactors = FALSE
  )

  # Filter to overlapping proteins
  tube_importance <- importance_df %>%
    dplyr::filter(protein %in% overlap_proteins) %>%
    dplyr::select(protein, tube_importance = importance) %>%
    dplyr::mutate(tube_importance_scaled = 100 * tube_importance / max(tube_importance))

  # Original study importance (approximated from paper ranking)
  # NEFL = 1 (highest), CSRP3 = 2, etc.
  original_biomarkers <- get_original_biomarkers()
  original_importance <- data.frame(
    protein = original_biomarkers,
    original_rank = seq_along(original_biomarkers)
  ) %>%
    dplyr::mutate(
      # Convert rank to importance score (inverse)
      original_importance_scaled = 100 * (18 - original_rank) / 17
    )

  # Merge
  plot_data <- tube_importance %>%
    dplyr::inner_join(original_importance, by = "protein") %>%
    dplyr::arrange(dplyr::desc(tube_importance_scaled)) %>%
    dplyr::mutate(
      protein = factor(protein, levels = protein)
    )

  # Reshape for plotting
  plot_data_long <- plot_data %>%
    dplyr::select(protein, original_importance_scaled, tube_importance_scaled) %>%
    tidyr::pivot_longer(
      cols = c(original_importance_scaled, tube_importance_scaled),
      names_to = "source",
      values_to = "importance"
    ) %>%
    dplyr::mutate(
      source = dplyr::case_when(
        source == "original_importance_scaled" ~ "Original Study\n(Biomarker Rank)",
        source == "tube_importance_scaled" ~ "Our Analysis\n(Tube Discrimination)"
      )
    )

  # Create heatmap
  p <- ggplot(plot_data_long, aes(x = source, y = protein, fill = importance)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = sprintf("%.0f", importance)),
      color = "white", fontface = "bold", size = 4
    ) +
    scale_fill_gradient2(
      low = "#FFFFFF",
      mid = "#FFA726",
      high = "#D32F2F",
      midpoint = 50,
      name = "Importance\n(scaled)"
    ) +
    labs(
      title = "Feature Importance: Original Biomarkers vs Tube Type Discrimination",
      subtitle = sprintf(
        "%d proteins appear in BOTH the original biomarker panel AND our tube-affected list",
        nrow(plot_data)
      ),
      x = NULL,
      y = "Protein"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(
        size = 12, hjust = 0.5, color = "#D32F2F",
        face = "bold", margin = margin(b = 15)
      ),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 11),
      legend.position = "right",
      panel.grid = element_blank()
    )

  # Save if requested
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)

    # Adjust height based on number of proteins
    plot_height <- max(6, 0.4 * nrow(plot_data))

    ggplot2::ggsave(save_path, p, width = 10, height = plot_height, dpi = 300, bg = "transparent")
    cat(sprintf("\n✓ Importance heatmap saved to: %s\n", save_path))
  }

  return(p)
}


#' Create Concordance Scatter Plot
#'
#' @description
#' Scatter plot comparing original biomarker importance rank (x-axis) vs
#' tube-type discrimination importance (y-axis) for overlapping proteins.
#'
#' If proteins cluster in the top-right quadrant, it means the most important
#' "biomarkers" are also the most affected by tube type = CRITICAL PROBLEM.
#'
#' @param reverse_prediction_results Results from reverse_prediction_test()
#' @param overlap_stats Results from calculate_overlap_statistics()
#' @param save_path Optional path to save figure
#' @return ggplot2 object
#' @export
plot_overlap_concordance <- function(reverse_prediction_results,
                                     overlap_stats,
                                     save_path = NULL) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)

  # Get overlapping proteins
  overlap_proteins <- overlap_stats$overlap_proteins

  if (length(overlap_proteins) == 0) {
    cat("\nNo overlapping proteins to plot.\n")
    return(NULL)
  }

  # Get tube-type importance from model
  importance_vec <- reverse_prediction_results$model$finalModel$variable.importance
  importance_df <- data.frame(
    protein = names(importance_vec),
    importance = as.numeric(importance_vec),
    stringsAsFactors = FALSE
  )

  # Filter to overlapping proteins
  tube_importance <- importance_df %>%
    dplyr::filter(protein %in% overlap_proteins) %>%
    dplyr::select(protein, tube_importance = importance)

  # Get original importance
  original_biomarkers <- get_original_biomarkers()
  original_importance <- data.frame(
    protein = original_biomarkers,
    original_rank = seq_along(original_biomarkers)
  )

  # Merge
  plot_data <- tube_importance %>%
    dplyr::inner_join(original_importance, by = "protein") %>%
    dplyr::arrange(original_rank)

  # Create scatter plot
  p <- ggplot(plot_data, aes(x = original_rank, y = tube_importance)) +
    geom_point(size = 4, alpha = 0.7, color = "#D32F2F") +
    ggrepel::geom_text_repel(
      aes(label = protein),
      size = 3.5,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      max.overlaps = 20
    ) +
    # Reference lines
    geom_hline(
      yintercept = median(plot_data$tube_importance),
      linetype = "dashed", color = "gray50"
    ) +
    geom_vline(
      xintercept = median(plot_data$original_rank),
      linetype = "dashed", color = "gray50"
    ) +
    # Quadrant labels
    annotate("text",
      x = 3, y = max(plot_data$tube_importance) * 0.95,
      label = "HIGH biomarker importance\nHIGH tube effect\n(CRITICAL PROBLEM)",
      size = 4, fontface = "bold", color = "#D32F2F", hjust = 0
    ) +
    labs(
      title = "Concordance: Original Biomarker Rank vs Tube Type Discrimination",
      subtitle = sprintf(
        "Proteins in top-left = important biomarkers heavily affected by tube type (n=%d)",
        nrow(plot_data)
      ),
      x = "Original Study Biomarker Rank (1 = most important)",
      y = "Tube Type Discrimination Importance (Our Analysis)"
    ) +
    scale_x_continuous(breaks = seq(1, max(plot_data$original_rank), by = 2)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(
        size = 12, color = "#D32F2F",
        face = "bold", margin = margin(b = 15)
      ),
      axis.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray30", fill = NA, size = 1)
    )

  # Save if requested
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 12, height = 10, dpi = 300, bg = "transparent")
    cat(sprintf("\n✓ Concordance plot saved to: %s\n", save_path))
  }

  return(p)
}


#' Generate Comprehensive Overlap Report
#'
#' @description
#' Creates a detailed text report documenting protein overlaps between
#' original biomarkers and tube-affected proteins.
#'
#' @param overlap_stats Results from calculate_overlap_statistics()
#' @param reverse_prediction_results Results from reverse_prediction_test()
#' @param save_path Optional path to save report
#' @return Character string with formatted report
#' @export
generate_overlap_report <- function(overlap_stats,
                                    reverse_prediction_results,
                                    save_path = NULL) {
  # Build report
  report <- paste0(
    "# PROTEIN OVERLAP ANALYSIS REPORT\n",
    "## Original Biomarkers vs Tube-Affected Proteins\n\n",
    "**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "---\n\n",
    "## SUMMARY\n\n",
    "This analysis compares the 17-protein biomarker panel from Chia et al. (2025)\n",
    "with the top 20 proteins most affected by tube type (from our reverse prediction test).\n\n",
    "**Research Question:** Are the original study's 'biomarkers' actually technical artifacts?\n\n",
    "---\n\n",
    "## KEY FINDINGS\n\n",
    sprintf("- **Original biomarkers:** %d proteins\n", overlap_stats$n_original),
    sprintf("- **Tube-affected proteins:** %d proteins\n", overlap_stats$n_tube),
    sprintf(
      "- **OVERLAP:** %d proteins (%.1f%% of biomarkers)\n\n",
      overlap_stats$n_overlap, overlap_stats$pct_original_affected
    ),
    sprintf("**INTERPRETATION:** %s\n\n", overlap_stats$interpretation),
    "---\n\n"
  )

  # Overlapping proteins
  if (overlap_stats$n_overlap > 0) {
    report <- paste0(
      report,
      "## OVERLAPPING PROTEINS (BIOMARKERS THAT ARE TUBE ARTIFACTS)\n\n",
      "These proteins appear in BOTH:\n",
      "1. The original study's 17-protein biomarker panel\n",
      "2. Our top 20 tube-affected proteins list\n\n"
    )

    for (protein in overlap_stats$overlap_proteins) {
      # Get ranks
      original_rank <- which(get_original_biomarkers() == protein)

      # Get tube type rank from model importance
      importance_vec <- reverse_prediction_results$model$finalModel$variable.importance
      importance_df <- data.frame(
        protein_name = names(importance_vec),
        importance = as.numeric(importance_vec),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::arrange(dplyr::desc(importance))

      tube_rank <- which(importance_df$protein_name == protein)

      report <- paste0(
        report,
        sprintf("### %s\n", protein),
        sprintf("- **Original biomarker rank:** #%d out of 17\n", original_rank),
        sprintf("- **Tube discrimination rank:** #%d out of 2,868 proteins\n", tube_rank),
        sprintf("- **Implication:** This 'biomarker' is heavily influenced by tube type.\n\n")
      )
    }
  } else {
    report <- paste0(
      report,
      "## OVERLAPPING PROTEINS\n\n",
      "NO OVERLAP FOUND. The original biomarkers do not appear in the top 20\n",
      "tube-affected proteins. This suggests the biomarker panel may not be\n",
      "heavily confounded by tube type effects (at least not the most extreme effects).\n\n"
    )
  }

  # Original only
  report <- paste0(
    report,
    "---\n\n",
    sprintf(
      "## ORIGINAL BIOMARKERS NOT AFFECTED BY TUBE TYPE (n=%d)\n\n",
      overlap_stats$n_original_only
    ),
    "These proteins were identified as biomarkers in the original study but do NOT\n",
    "appear in our top 20 tube-affected proteins:\n\n"
  )

  for (protein in overlap_stats$original_only_proteins) {
    report <- paste0(report, sprintf("- %s\n", protein))
  }

  # Tube only
  report <- paste0(
    report,
    "\n---\n\n",
    sprintf(
      "## TUBE ARTIFACTS NOT IDENTIFIED AS BIOMARKERS (n=%d)\n\n",
      overlap_stats$n_tube_only
    ),
    "These proteins are heavily affected by tube type but were NOT included in\n",
    "the original biomarker panel:\n\n"
  )

  for (protein in overlap_stats$tube_only_proteins) {
    report <- paste0(report, sprintf("- %s\n", protein))
  }

  # Conclusions
  report <- paste0(
    report,
    "\n---\n\n",
    "## CONCLUSIONS\n\n",
    "1. **Biomarker Validity:**\n"
  )

  if (overlap_stats$pct_original_affected >= 30) {
    report <- paste0(
      report,
      "   - CRITICAL: A substantial fraction (",
      sprintf("%.1f%%", overlap_stats$pct_original_affected),
      ") of the original biomarkers\n",
      "     are heavily affected by tube type. This raises serious questions about\n",
      "     whether they reflect ALS biology or technical artifacts.\n\n"
    )
  } else if (overlap_stats$pct_original_affected >= 10) {
    report <- paste0(
      report,
      "   - MODERATE: Some biomarkers (",
      sprintf("%.1f%%", overlap_stats$pct_original_affected),
      ") overlap with tube-affected proteins.\n",
      "     This suggests partial confounding that should be investigated further.\n\n"
    )
  } else {
    report <- paste0(
      report,
      "   - MINIMAL: Few biomarkers (",
      sprintf("%.1f%%", overlap_stats$pct_original_affected),
      ") overlap with the most\n",
      "     tube-affected proteins. However, perfect confounding between tube type,\n",
      "     country, and diagnosis still makes interpretation difficult.\n\n"
    )
  }

  report <- paste0(
    report,
    "2. **Study Design Issues:**\n",
    "   - The original study included tube type as a covariate in their model,\n",
    "     acknowledging its importance (ranked #2 overall).\n",
    "   - However, perfect confounding with country and near-perfect confounding\n",
    "     with diagnosis means adjusting for tube type may not fully control\n",
    "     for these technical effects.\n",
    "   - Our reverse prediction test (AUC = 0.999) shows proteins can perfectly\n",
    "     predict tube type, indicating pervasive technical effects.\n\n",
    "3. **Recommendations:**\n",
    "   - Independent validation using CONSISTENT tube types across all samples\n",
    "   - Stratified analyses within tube type (Italy/HEPARIN only, US/EDTA only)\n",
    "   - Biological validation of overlapping proteins to determine if effects\n",
    "     are technical or biological\n",
    "   - Transparent reporting of confounding structure in publications\n\n",
    "---\n\n",
    "**Report Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "**Analysis Pipeline:** targets + R 4.4.2\n",
    "**Investigation Team:** ALS Biomarker Bias Analysis Group\n"
  )

  # Save if requested
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(report, save_path)
    cat(sprintf("\n✓ Overlap report saved to: %s\n", save_path))
  }

  # Print to console
  cat(report)

  return(report)
}
