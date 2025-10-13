# Test Differential Analysis Functions
# Tests for limma implementation and meta-analysis correctness

test_that("limma results have expected structure and valid values", {
  skip("Requires full dataset - run manually")

  # This test validates that differential analysis returns correct structure
  # italy_results <- differential_analysis_italy(protein_long)

  # Required columns
  # required_cols <- c("Assay", "logFC", "P.Value", "adj.P.Val", "AveExpr", "t", "B")
  # expect_true(all(required_cols %in% names(italy_results)))

  # FDR control: adjusted p-values should be >= raw p-values
  # expect_true(all(italy_results$adj.P.Val >= italy_results$P.Value))

  # All p-values should be in [0, 1]
  # expect_true(all(italy_results$P.Value >= 0 & italy_results$P.Value <= 1))
  # expect_true(all(italy_results$adj.P.Val >= 0 & italy_results$adj.P.Val <= 1))
})

test_that("meta-analysis combined p-values are valid", {
  # Test Fisher's method for combining p-values

  # Create mock results
  mock_italy <- tibble::tibble(
    Assay = c("Protein1", "Protein2", "Protein3"),
    logFC = c(1.5, -1.2, 0.8),
    P.Value = c(0.001, 0.05, 0.10),
    adj.P.Val = c(0.01, 0.15, 0.30),
    AveExpr = c(5, 6, 7)
  )

  mock_us <- tibble::tibble(
    Assay = c("Protein1", "Protein2", "Protein3"),
    logFC = c(1.4, -1.1, -0.3),  # Protein3 has opposite direction
    P.Value = c(0.002, 0.06, 0.12),
    adj.P.Val = c(0.02, 0.18, 0.35),
    AveExpr = c(5.5, 6.2, 6.8)
  )

  # Test directional concordance
  merged <- dplyr::inner_join(mock_italy, mock_us, by = "Assay", suffix = c("_italy", "_us"))

  concordant <- sign(merged$logFC_italy) == sign(merged$logFC_us)

  expect_equal(concordant, c(TRUE, TRUE, FALSE))

  # Test Fisher's method calculation
  # Combined p-value should be smaller than individual p-values for concordant proteins
  for (i in 1:nrow(merged)) {
    chi_stat <- -2 * (log(merged$P.Value_italy[i]) + log(merged$P.Value_us[i]))
    combined_p <- pchisq(chi_stat, df = 4, lower.tail = FALSE)

    if (concordant[i]) {
      expect_lte(combined_p, min(merged$P.Value_italy[i], merged$P.Value_us[i]),
                 label = "Combined p-value should be <= min individual p-value")
    }
  }
})

test_that("stratified vs pooled comparison calculates overlaps correctly", {
  # Test the overlap calculation logic with mock data

  mock_sig_italy <- c("P1", "P2", "P3", "P4", "P5")
  mock_sig_us <- c("P2", "P3", "P6", "P7")
  mock_sig_pooled <- c("P1", "P2", "P3", "P8", "P9", "P10")

  # Both strata
  sig_both <- intersect(mock_sig_italy, mock_sig_us)
  expect_equal(sig_both, c("P2", "P3"))

  # Pooled overlap with both
  pooled_overlap_both <- intersect(mock_sig_pooled, sig_both)
  expect_equal(pooled_overlap_both, c("P2", "P3"))

  # Pooled only (not in either stratum)
  pooled_only <- setdiff(mock_sig_pooled, union(mock_sig_italy, mock_sig_us))
  expect_equal(pooled_only, c("P8", "P9", "P10"))

  # Percentage replication
  pct_replicate <- 100 * length(pooled_overlap_both) / length(mock_sig_pooled)
  expect_equal(pct_replicate, 100 * 2 / 6, tolerance = 0.1)
})

test_that("significance thresholds are applied correctly", {
  # Test that FDR and logFC thresholds work as expected

  mock_results <- tibble::tibble(
    Assay = paste0("P", 1:10),
    logFC = c(0.8, -0.7, 0.3, -0.9, 1.2, 0.4, -1.5, 0.2, 0.6, -0.45),
    adj.P.Val = c(0.01, 0.02, 0.08, 0.001, 0.03, 0.10, 0.005, 0.20, 0.04, 0.06)
  )

  # Apply thresholds
  sig <- mock_results %>%
    dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 0.5)

  # Should get proteins with adj.P.Val < 0.05 AND |logFC| > 0.5
  # P1: logFC=0.8 (>0.5), adj.P=0.01 (<0.05) ✓
  # P2: logFC=-0.7 (abs>0.5), adj.P=0.02 (<0.05) ✓
  # P4: logFC=-0.9 (abs>0.5), adj.P=0.001 (<0.05) ✓
  # P5: logFC=1.2 (>0.5), adj.P=0.03 (<0.05) ✓
  # P7: logFC=-1.5 (abs>0.5), adj.P=0.005 (<0.05) ✓
  # P9: logFC=0.6 (>0.5), adj.P=0.04 (<0.05) ✓
  expected <- c("P1", "P2", "P4", "P5", "P7", "P9")
  expect_equal(sig$Assay, expected)
})
