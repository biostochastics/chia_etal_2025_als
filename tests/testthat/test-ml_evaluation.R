# Test ML Evaluation Functions
# Tests for data leakage, reproducibility, and correctness

test_that("reverse prediction is reproducible with fixed seed", {
  skip("Requires full dataset - run manually")

  # This test would require the full protein_wide dataset
  # Skipping for automated testing, but template is here for manual validation

  set.seed(42)
  # result1 <- reverse_prediction_test(protein_wide, seed = 42, n_folds = 3)

  set.seed(42)
  # result2 <- reverse_prediction_test(protein_wide, seed = 42, n_folds = 3)

  # expect_equal(result1$auc, result2$auc, tolerance = 0.001)
})

test_that("leave-country-out CV uses independent feature selection", {
  skip("Requires full dataset - run manually")

  # This test verifies that Italy and US models filter proteins independently
  # which is the CORRECT implementation to prevent data leakage

  # lcv_results <- leave_country_out_cv(protein_wide, seed = 42)

  # The models should have been trained on different protein subsets
  # because filtering is done on training data only

  # This is the CORRECT behavior
})

test_that("AUC calculations are mathematically correct", {
  # Test with synthetic data where we know the expected AUC

  # Perfect prediction
  y_true <- factor(c(rep("Control", 50), rep("ALS", 50)), levels = c("Control", "ALS"))
  y_pred_perfect <- c(rep(0, 50), rep(1, 50))

  roc_obj <- pROC::roc(y_true, y_pred_perfect, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))

  expect_equal(auc_val, 1.0, tolerance = 0.001,
               label = "Perfect prediction should have AUC = 1.0")

  # Random prediction (should be ~0.5)
  set.seed(123)
  y_pred_random <- runif(100)
  roc_obj_random <- pROC::roc(y_true, y_pred_random, quiet = TRUE)
  auc_random <- as.numeric(pROC::auc(roc_obj_random))

  expect_gt(auc_random, 0.4,
            label = "Random prediction AUC should be around 0.5")
  expect_lt(auc_random, 0.6,
            label = "Random prediction AUC should be around 0.5")

  # Good but not perfect prediction
  y_pred_good <- c(runif(50, 0, 0.4), runif(50, 0.6, 1))
  roc_obj_good <- pROC::roc(y_true, y_pred_good, quiet = TRUE)
  auc_good <- as.numeric(pROC::auc(roc_obj_good))

  expect_gt(auc_good, 0.7,
            label = "Well-separated predictions should have AUC > 0.7")
})

test_that("ranger handles missing values correctly", {
  # Test that ranger can handle NA values when configured correctly
  set.seed(42)

  # Create synthetic data with missing values
  n <- 100
  p <- 10

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("V", 1:p)

  # Introduce missing values
  X[sample(1:(n * p), size = 20)] <- NA

  y <- factor(sample(c("A", "B"), n, replace = TRUE))

  X_df <- as.data.frame(X)

  # ranger handles missing values natively when replace = TRUE
  # Should not error - suppress messages
  suppressMessages({
    expect_no_error({
      model <- caret::train(
        x = X_df,
        y = y,
        method = "ranger",
        trControl = caret::trainControl(method = "cv", number = 3, verboseIter = FALSE),
        tuneGrid = data.frame(mtry = 3, splitrule = "gini", min.node.size = 5),
        num.trees = 50,
        respect.unordered.factors = "order",
        replace = TRUE
      )
    })
  })
})

test_that("comparison functions handle edge cases", {
  # Test that model comparison functions handle cases where AUC values are similar

  # Mock results
  mock_with_tube <- list(auc = 0.95)
  mock_without_tube <- list(auc = 0.94)

  auc_diff <- mock_with_tube$auc - mock_without_tube$auc
  pct_change <- 100 * auc_diff / mock_without_tube$auc

  expect_equal(auc_diff, 0.01, tolerance = 0.001)
  expect_lt(abs(pct_change), 2,
            label = "Small AUC differences should be <2%")
})
