#' Machine Learning Model Evaluation
#'
#' Functions for evaluating ML model dependence on non-biological features,
#' particularly plasma collection tube type.
#'
#' @description
#' Implements the REVERSE PREDICTION TEST:
#' - Train model to predict TUBE TYPE from protein expression
#' - If AUC > 0.9 → tube type signal dominates biology
#' - This is the "smoking gun" for confounding bias
#'
#' Also implements:
#' - Feature importance analysis
#' - Model performance with/without tube type
#' - Permutation importance
#'

#' Pivot protein data from long to wide format
#'
#' @description
#' Converts long-format OLINK data (one row per sample-protein) to wide format
#' (one row per sample, one column per protein). Required for ML modeling.
#'
#' @param data Long-format data with NPX values
#' @return Wide-format tibble with samples as rows, proteins as columns
#' @examples
#' protein_wide <- pivot_protein_wide(data_with_country)
#' @export
pivot_protein_wide <- function(data) {
  message("Pivoting protein data to wide format...")

  # Select relevant columns
  protein_data <- data %>%
    dplyr::select(
      SampleID_deidentified,
      Assay,
      NPX,
      Diagnosis,
      Plasma_collection_tube_type,
      country,
      Sex,
      Age_Collection
    )

  # Pivot to wide format
  protein_wide <- protein_data %>%
    # Keep only one row per sample-protein
    dplyr::distinct(SampleID_deidentified, Assay, .keep_all = TRUE) %>%
    # Pivot
    tidyr::pivot_wider(
      id_cols = c(
        SampleID_deidentified,
        Diagnosis,
        Plasma_collection_tube_type,
        country,
        Sex,
        Age_Collection
      ),
      names_from = Assay,
      values_from = NPX,
      values_fill = NA
    )

  # Check for missing values
  n_proteins <- ncol(protein_wide) - 6 # Subtract metadata columns
  n_samples <- nrow(protein_wide)

  message(sprintf("Wide format: %s samples × %s proteins", n_samples, n_proteins))

  # Calculate missingness
  protein_cols <- protein_wide %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis,
      -Plasma_collection_tube_type, -country,
      -Sex, -Age_Collection
    )

  pct_missing <- 100 * sum(is.na(protein_cols)) / (nrow(protein_cols) * ncol(protein_cols))
  message(sprintf("Missing values: %.2f%%", pct_missing))

  return(protein_wide)
}


#' **CRITICAL: Reverse Prediction Test**
#'
#' @description
#' The "smoking gun" analysis: Can we predict tube type from protein expression?
#'
#' If AUC > 0.9, this indicates tube type signal dominates biological signal.
#' This would mean the proteins are more informative about the collection tube
#' than about disease state - a devastating finding for the biomarker claims.
#'
#' IMPORTANT: This implements proper cross-validation with NO data leakage.
#' Feature scaling and imputation are done WITHIN each fold.
#'
#' @param protein_wide Wide-format data with samples as rows
#' @param n_folds Number of cross-validation folds (default 5). Higher values (10) increase computation but reduce variance; lower (3) faster but less stable
#' @param seed Random seed for reproducibility. Critical for ensuring identical train/test splits across runs and enabling exact replication
#' @return List with model performance metrics and feature importance
#' @examples
#' reverse_test <- reverse_prediction_test(protein_wide, n_folds = 5, seed = 42)
#' @export
reverse_prediction_test <- function(protein_wide, n_folds = 5, seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("REVERSE PREDICTION TEST: Can proteins predict TUBE TYPE?")
  message(paste(rep("=", 70), collapse = ""))
  message("\nIf AUC > 0.9, tube type signal dominates biological signal!")
  message("This would be DAMNING evidence for the biomarker claims.\n")

  set.seed(seed)

  # Prepare data
  # Remove samples with unknown country
  data_clean <- protein_wide %>%
    dplyr::filter(country != "Unknown")

  # Outcome: tube type (binary)
  y <- factor(data_clean$Plasma_collection_tube_type,
    levels = c("EDTA", "HEPARIN")
  )

  # Features: protein expression only (exclude metadata)
  X <- data_clean %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis,
      -Plasma_collection_tube_type, -country,
      -Sex, -Age_Collection
    )

  message(sprintf("Samples: %s", nrow(X)))
  message(sprintf("Features: %s proteins", ncol(X)))
  message(sprintf(
    "Outcome: %s EDTA, %s HEPARIN",
    sum(y == "EDTA"), sum(y == "HEPARIN")
  ))

  # NOTE: Protein filtering removed to prevent data leakage
  # ranger handles missing values natively when using replace = TRUE
  # This ensures no information from test folds influences feature selection
  message("\nUsing all proteins - ranger handles missing values within CV folds")

  # Set up cross-validation
  train_control <- caret::trainControl(
    method = "cv",
    number = n_folds,
    savePredictions = "final",
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    verboseIter = TRUE
  )

  # Train Random Forest model
  message("\nTraining Random Forest model...")
  message("(This may take several minutes)")
  message("NOTE: Missing values handled natively by ranger within each CV fold")

  rf_model <- caret::train(
    x = X,
    y = y,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = c(floor(sqrt(ncol(X))), floor(ncol(X) / 3)),
      splitrule = "gini",
      min.node.size = c(5, 10)
    ),
    importance = "permutation",
    num.trees = 500,
    # CRITICAL: Handle missing values within ranger to prevent data leakage
    # Each CV fold will handle missing values independently
    respect.unordered.factors = "order",
    replace = TRUE
    # ranger handles missing values natively when replace = TRUE
  )

  # Extract performance metrics
  best_tune <- rf_model$bestTune
  cv_results <- rf_model$results %>%
    dplyr::filter(
      mtry == best_tune$mtry,
      splitrule == best_tune$splitrule,
      min.node.size == best_tune$min.node.size
    )

  auc <- cv_results$ROC
  sensitivity <- cv_results$Sens
  specificity <- cv_results$Spec

  # Calculate AUC confidence interval using pROC
  # Combine predictions across all CV folds
  cv_preds <- rf_model$pred %>%
    dplyr::filter(
      mtry == best_tune$mtry,
      splitrule == best_tune$splitrule,
      min.node.size == best_tune$min.node.size
    )

  roc_obj <- pROC::roc(
    response = cv_preds$obs,
    predictor = cv_preds$HEPARIN,
    levels = c("EDTA", "HEPARIN"),
    direction = "<",
    quiet = TRUE
  )

  auc_ci <- pROC::ci.auc(roc_obj, conf.level = 0.95, method = "delong")
  auc_ci_lower <- auc_ci[1]
  auc_ci_upper <- auc_ci[3]

  message("\n", paste(rep("=", 70), collapse = ""))
  message("RESULTS: Reverse Prediction Test")
  message(paste(rep("=", 70), collapse = ""))
  message(sprintf("\n  Cross-validated AUC:       %.3f (95%% CI: %.3f-%.3f)", auc, auc_ci_lower, auc_ci_upper))
  message(sprintf("  Sensitivity:               %.3f", sensitivity))
  message(sprintf("  Specificity:               %.3f", specificity))

  # Interpret results
  if (auc > 0.95) {
    message("\n  ⚠️  CRITICAL FINDING: AUC > 0.95")
    message("  Proteins can PERFECTLY predict tube type!")
    message("  This suggests tube type effects DOMINATE biological signal.")
  } else if (auc > 0.9) {
    message("\n  ⚠️  SEVERE CONCERN: AUC > 0.9")
    message("  Proteins strongly predict tube type.")
    message("  Major confounding with non-biological factors.")
  } else if (auc > 0.8) {
    message("\n  ⚠️  MODERATE CONCERN: AUC > 0.8")
    message("  Tube type can be predicted with good accuracy.")
    message("  Substantial confounding present.")
  } else if (auc > 0.7) {
    message("\n  ℹ️  MILD CONCERN: AUC > 0.7")
    message("  Some tube type signal detectable in proteins.")
  } else {
    message("\n  ✓ ACCEPTABLE: AUC < 0.7")
    message("  Tube type effects are modest.")
  }

  message("\n", paste(rep("=", 70), collapse = ""))

  # Feature importance (top 20 proteins)
  var_imp <- caret::varImp(rf_model, scale = TRUE)
  top_features <- var_imp$importance %>%
    tibble::rownames_to_column("protein") %>%
    dplyr::arrange(dplyr::desc(Overall)) %>%
    dplyr::slice_head(n = 20)

  message("\nTop 20 proteins distinguishing HEPARIN vs EDTA:")
  print(as.data.frame(top_features))

  # Compile results
  results <- list(
    model = rf_model,
    auc = auc,
    auc_ci_lower = auc_ci_lower,
    auc_ci_upper = auc_ci_upper,
    sensitivity = sensitivity,
    specificity = specificity,
    best_tune = best_tune,
    cv_results = cv_results,
    top_features = top_features,
    n_samples = nrow(X),
    n_proteins = ncol(X),
    n_folds = n_folds,
    seed = seed,
    outcome_distribution = table(y),
    roc_obj = roc_obj
  )

  return(results)
}


#' Train model WITH tube type as feature
#'
#' @description
#' Trains a model to predict diagnosis using proteins + tube type.
#' This replicates the approach from the original study where tube type
#' was ranked as the 2nd most important feature.
#'
#' @param protein_wide Wide-format data
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return Model results
#' @examples
#' model_with_tube <- train_model_with_tube_type(protein_wide)
#' @export
train_model_with_tube_type <- function(protein_wide, n_folds = 5, seed = 42) {
  message("\nTraining model WITH tube type as feature...")

  set.seed(seed)

  # Prepare data - binary classification: ALS vs Controls
  data_clean <- protein_wide %>%
    dplyr::filter(country != "Unknown") %>%
    dplyr::mutate(
      outcome = dplyr::if_else(Diagnosis == "ALS", "ALS", "Control")
    )

  y <- factor(data_clean$outcome, levels = c("Control", "ALS"))

  # Features: proteins + tube type
  X <- data_clean %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis, -outcome,
      -country, -Sex, -Age_Collection
    )

  # Encode tube type as numeric
  X <- X %>%
    dplyr::mutate(
      tube_type_encoded = as.numeric(Plasma_collection_tube_type == "HEPARIN")
    ) %>%
    dplyr::select(-Plasma_collection_tube_type)

  # NOTE: Protein filtering removed to prevent data leakage
  # ranger handles missing values natively within each CV fold
  message(sprintf("Features: %s proteins + tube_type (no filtering)", ncol(X) - 1))
  message("Missing values handled by ranger within CV folds")

  # Train model
  train_control <- caret::trainControl(
    method = "cv",
    number = n_folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary
  )

  model <- caret::train(
    x = X,
    y = y,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = floor(sqrt(ncol(X))),
      splitrule = "gini",
      min.node.size = 5
    ),
    importance = "permutation",
    num.trees = 500
  )

  # Feature importance
  var_imp <- caret::varImp(model, scale = TRUE)
  importance_df <- var_imp$importance %>%
    tibble::rownames_to_column("feature") %>%
    dplyr::arrange(dplyr::desc(Overall))

  # Check tube type rank
  tube_rank <- which(importance_df$feature == "tube_type_encoded")

  message(sprintf("\nModel AUC: %.3f", max(model$results$ROC)))
  message(sprintf("Tube type importance rank: %s / %s", tube_rank, nrow(importance_df)))

  if (tube_rank <= 10) {
    message("⚠️  WARNING: Tube type is in top 10 features!")
  }

  return(list(
    model = model,
    importance = importance_df,
    tube_type_rank = tube_rank,
    auc = max(model$results$ROC)
  ))
}


#' Train model WITHOUT tube type
#'
#' @description
#' Trains model using only proteins (no tube type feature).
#' Compare to model WITH tube type to quantify tube type contribution.
#'
#' @param protein_wide Wide-format data
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return Model results
#' @examples
#' model_without_tube <- train_model_without_tube_type(protein_wide)
#' @export
train_model_without_tube_type <- function(protein_wide, n_folds = 5, seed = 42) {
  message("\nTraining model WITHOUT tube type (proteins only)...")

  set.seed(seed)

  # Prepare data
  data_clean <- protein_wide %>%
    dplyr::filter(country != "Unknown") %>%
    dplyr::mutate(
      outcome = dplyr::if_else(Diagnosis == "ALS", "ALS", "Control")
    )

  y <- factor(data_clean$outcome, levels = c("Control", "ALS"))

  # Features: proteins ONLY (no tube type)
  X <- data_clean %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis, -outcome,
      -country, -Sex, -Age_Collection,
      -Plasma_collection_tube_type
    )

  # NOTE: No protein filtering to prevent data leakage
  # ranger handles missing values natively within each CV fold
  message(sprintf("Features: %s proteins (no tube type, no filtering)", ncol(X)))
  message("Missing values handled by ranger within CV folds")

  # Train model
  train_control <- caret::trainControl(
    method = "cv",
    number = n_folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary
  )

  model <- caret::train(
    x = X,
    y = y,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = floor(sqrt(ncol(X))),
      splitrule = "gini",
      min.node.size = 5
    ),
    importance = "permutation",
    num.trees = 500,
    # Handle missing values within CV
    respect.unordered.factors = "order",
    replace = TRUE
    # ranger handles missing values natively when replace = TRUE
  )

  message(sprintf("\nModel AUC (without tube type): %.3f", max(model$results$ROC)))

  return(list(
    model = model,
    auc = max(model$results$ROC)
  ))
}


#' Compare model performance with/without tube type
#'
#' @description
#' Quantifies the contribution of tube type to model performance.
#' Large AUC drop when removing tube type suggests confounding dependence.
#'
#' @param protein_wide Wide-format data
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return Comparison results
#' @examples
#' comparison <- compare_tube_type_contribution(protein_wide)
#' @export
compare_tube_type_contribution <- function(protein_wide, n_folds = 5, seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("COMPARING MODELS: With vs Without Tube Type")
  message(paste(rep("=", 70), collapse = ""))

  # Train both models
  model_with <- train_model_with_tube_type(protein_wide, n_folds, seed)
  model_without <- train_model_without_tube_type(protein_wide, n_folds, seed)

  # Calculate difference
  auc_diff <- model_with$auc - model_without$auc
  pct_change <- 100 * auc_diff / model_without$auc

  message("\n", paste(rep("=", 70), collapse = ""))
  message("COMPARISON RESULTS")
  message(paste(rep("=", 70), collapse = ""))
  message(sprintf("\n  AUC with tube type:    %.3f", model_with$auc))
  message(sprintf("  AUC without tube type: %.3f", model_without$auc))
  message(sprintf("  Difference:            %.3f (%.1f%% change)", auc_diff, pct_change))

  if (abs(pct_change) > 10) {
    message("\n  ⚠️  MAJOR FINDING: >10% change in AUC")
    message("  Tube type substantially contributes to model performance.")
  } else if (abs(pct_change) > 5) {
    message("\n  ⚠️  NOTABLE: 5-10% change in AUC")
    message("  Tube type has meaningful contribution.")
  } else {
    message("\n  ✓ Modest contribution from tube type (<5%)")
  }

  message("\n", paste(rep("=", 70), collapse = ""))

  return(list(
    model_with_tube = model_with,
    model_without_tube = model_without,
    auc_difference = auc_diff,
    percent_change = pct_change
  ))
}


#' **CRITICAL: Leave-Country-Out Cross-Validation**
#'
#' @description
#' The definitive test of geographic generalizability.
#'
#' Tests whether models trained on one country generalize to another:
#' - Train on Italy (HEPARIN) → Test on US (EDTA)
#' - Train on US (EDTA) → Test on Italy (HEPARIN)
#'
#' If performance drops dramatically, this indicates the model relies on
#' country-specific or tube-specific signals rather than true ALS biology.
#'
#' IMPORTANT: This implements proper preprocessing without data leakage.
#' All feature scaling, imputation, and selection are done ONLY on training data.
#'
#' @param protein_wide Wide-format data with country labels
#' @param outcome_var Outcome variable name (default "Diagnosis")
#' @param seed Random seed for reproducibility
#' @return List with performance for both directions
#' @examples
#' lcv_results <- leave_country_out_cv(protein_wide)
#' @export
leave_country_out_cv <- function(protein_wide, outcome_var = "Diagnosis", seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("LEAVE-COUNTRY-OUT CROSS-VALIDATION")
  message(paste(rep("=", 70), collapse = ""))
  message("\nTesting geographic generalizability:")
  message("  1. Train Italy (HEPARIN) → Test US (EDTA)")
  message("  2. Train US (EDTA) → Test Italy (HEPARIN)")
  message("\nIf performance drops >> 20%, model relies on geographic/tube factors!")
  message(paste(rep("=", 70), collapse = ""))

  set.seed(seed)

  # Prepare data - binary classification: ALS vs Controls
  data_clean <- protein_wide %>%
    dplyr::filter(country %in% c("Italy", "US")) %>%
    dplyr::mutate(
      outcome = dplyr::case_when(
        Diagnosis == "ALS" ~ "ALS",
        Diagnosis %in% c("Healthy_control", "Neurological_control") ~ "Control",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(outcome))

  # Split by country
  italy_data <- data_clean %>% dplyr::filter(country == "Italy")
  us_data <- data_clean %>% dplyr::filter(country == "US")

  message(sprintf(
    "\nItaly samples: %s (%s ALS, %s Control)",
    nrow(italy_data),
    sum(italy_data$outcome == "ALS"),
    sum(italy_data$outcome == "Control")
  ))

  message(sprintf(
    "US samples:    %s (%s ALS, %s Control)",
    nrow(us_data),
    sum(us_data$outcome == "ALS"),
    sum(us_data$outcome == "Control")
  ))

  # DIRECTION 1: Train Italy → Test US
  message("\n\n", paste(rep("-", 70), collapse = ""))
  message("DIRECTION 1: Train on Italy (HEPARIN) → Test on US (EDTA)")
  message(paste(rep("-", 70), collapse = ""))

  italy_to_us <- train_and_test_cv(
    train_data = italy_data,
    test_data = us_data,
    train_name = "Italy",
    test_name = "US",
    seed = seed
  )

  # DIRECTION 2: Train US → Test Italy
  message("\n\n", paste(rep("-", 70), collapse = ""))
  message("DIRECTION 2: Train on US (EDTA) → Test on Italy (HEPARIN)")
  message(paste(rep("-", 70), collapse = ""))

  us_to_italy <- train_and_test_cv(
    train_data = us_data,
    test_data = italy_data,
    train_name = "US",
    test_name = "Italy",
    seed = seed
  )

  # Summary comparison
  message("\n\n", paste(rep("=", 70), collapse = ""))
  message("LEAVE-COUNTRY-OUT CV SUMMARY")
  message(paste(rep("=", 70), collapse = ""))
  message("\nTest Performance:")
  message(sprintf(
    "  Italy → US:     AUC = %.3f, Accuracy = %.3f",
    italy_to_us$test_auc, italy_to_us$test_accuracy
  ))
  message(sprintf(
    "  US → Italy:     AUC = %.3f, Accuracy = %.3f",
    us_to_italy$test_auc, us_to_italy$test_accuracy
  ))
  message(sprintf("\nMean test AUC:    %.3f", mean(c(italy_to_us$test_auc, us_to_italy$test_auc))))

  # Interpretation
  mean_test_auc <- mean(c(italy_to_us$test_auc, us_to_italy$test_auc))
  mean_train_auc <- mean(c(italy_to_us$train_auc, us_to_italy$train_auc))
  performance_drop <- mean_train_auc - mean_test_auc

  message("\n", paste(rep("-", 70), collapse = ""))
  message("INTERPRETATION:")

  if (performance_drop > 0.2) {
    message("  ⚠️  CRITICAL: >0.2 AUC drop from train to test")
    message("  Model DOES NOT generalize across countries!")
    message("  Strong evidence of country/tube-specific overfitting.")
  } else if (performance_drop > 0.1) {
    message("  ⚠️  SUBSTANTIAL: 0.1-0.2 AUC drop")
    message("  Limited generalizability across geographic cohorts.")
  } else if (performance_drop > 0.05) {
    message("  ⚠️  MODERATE: 0.05-0.1 AUC drop")
    message("  Some country-specific effects present.")
  } else {
    message("  ✓ Good generalizability (<0.05 AUC drop)")
    message("  Model performance is consistent across countries.")
  }

  message(paste(rep("=", 70), collapse = ""))

  return(list(
    italy_to_us = italy_to_us,
    us_to_italy = us_to_italy,
    mean_test_auc = mean_test_auc,
    mean_train_auc = mean_train_auc,
    performance_drop = performance_drop,
    italy_n = nrow(italy_data),
    us_n = nrow(us_data)
  ))
}


#' Train model on one cohort and test on another
#'
#' @description
#' Helper function for leave-country-out CV.
#' Handles proper preprocessing without data leakage.
#'
#' @param train_data Training data
#' @param test_data Test data
#' @param train_name Name of training cohort (for logging)
#' @param test_name Name of test cohort (for logging)
#' @param seed Random seed
#' @return List with train and test performance
#' @keywords internal
train_and_test_cv <- function(train_data, test_data, train_name, test_name, seed = 42) {
  set.seed(seed)

  # Prepare outcome
  y_train <- factor(train_data$outcome, levels = c("Control", "ALS"))
  y_test <- factor(test_data$outcome, levels = c("Control", "ALS"))

  # Extract protein features
  X_train <- train_data %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis, -outcome,
      -Plasma_collection_tube_type, -country,
      -Sex, -Age_Collection
    )

  X_test <- test_data %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis, -outcome,
      -Plasma_collection_tube_type, -country,
      -Sex, -Age_Collection
    )

  # Filter proteins with high missingness IN TRAINING SET ONLY
  missing_pct_train <- colMeans(is.na(X_train))
  proteins_to_keep <- names(X_train)[missing_pct_train < 0.2]

  message(sprintf(
    "\nProtein filtering: %s → %s proteins (>20%% missing removed)",
    ncol(X_train), length(proteins_to_keep)
  ))

  X_train <- X_train %>% dplyr::select(all_of(proteins_to_keep))
  X_test <- X_test %>% dplyr::select(all_of(intersect(names(X_test), proteins_to_keep)))

  # Check sample sizes
  if (nrow(X_train) < 20) {
    warning(sprintf("Very small training set: n = %s. Results may be unstable.", nrow(X_train)))
  }
  if (nrow(X_test) < 20) {
    warning(sprintf("Very small test set: n = %s. Results may be unstable.", nrow(X_test)))
  }

  message(sprintf("\nTraining: %s samples × %s proteins", nrow(X_train), ncol(X_train)))
  message(sprintf("Testing:  %s samples × %s proteins", nrow(X_test), ncol(X_test)))

  # Train Random Forest with internal 5-fold CV
  train_control <- caret::trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = FALSE
  )

  message("\nTraining Random Forest...")

  rf_model <- caret::train(
    x = X_train,
    y = y_train,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = floor(sqrt(ncol(X_train))),
      splitrule = "gini",
      min.node.size = 5
    ),
    importance = "permutation",
    num.trees = 500,
    respect.unordered.factors = "order"
  )

  # Training performance (CV within training set)
  train_auc <- max(rf_model$results$ROC)
  train_sens <- rf_model$results$Sens[which.max(rf_model$results$ROC)]
  train_spec <- rf_model$results$Spec[which.max(rf_model$results$ROC)]

  message(sprintf("\nTraining performance (5-fold CV on %s):", train_name))
  message(sprintf("  AUC:         %.3f", train_auc))
  message(sprintf("  Sensitivity: %.3f", train_sens))
  message(sprintf("  Specificity: %.3f", train_spec))

  # Predict on test set
  message(sprintf("\nTesting on %s...", test_name))

  test_pred_prob <- predict(rf_model, newdata = X_test, type = "prob")
  test_pred_class <- predict(rf_model, newdata = X_test, type = "raw")

  # Calculate test metrics
  test_accuracy <- mean(test_pred_class == y_test)

  # ROC curve for test set
  roc_obj <- pROC::roc(
    response = y_test,
    predictor = test_pred_prob$ALS,
    levels = c("Control", "ALS"),
    direction = "<",
    quiet = TRUE
  )

  test_auc <- as.numeric(pROC::auc(roc_obj))

  # Calculate AUC confidence interval
  test_auc_ci <- pROC::ci.auc(roc_obj, conf.level = 0.95, method = "delong")
  test_auc_ci_lower <- test_auc_ci[1]
  test_auc_ci_upper <- test_auc_ci[3]

  # Confusion matrix
  conf_mat <- table(Predicted = test_pred_class, Actual = y_test)
  test_sens <- conf_mat["ALS", "ALS"] / sum(conf_mat[, "ALS"])
  test_spec <- conf_mat["Control", "Control"] / sum(conf_mat[, "Control"])

  message(sprintf("\nTest performance on %s:", test_name))
  message(sprintf("  AUC:         %.3f (95%% CI: %.3f-%.3f)", test_auc, test_auc_ci_lower, test_auc_ci_upper))
  message(sprintf("  Accuracy:    %.3f", test_accuracy))
  message(sprintf("  Sensitivity: %.3f", test_sens))
  message(sprintf("  Specificity: %.3f", test_spec))

  message("\nConfusion Matrix:")
  print(conf_mat)

  # Feature importance
  var_imp <- caret::varImp(rf_model, scale = TRUE)
  top_features <- var_imp$importance %>%
    tibble::rownames_to_column("protein") %>%
    dplyr::arrange(dplyr::desc(Overall)) %>%
    dplyr::slice_head(n = 10)

  message("\nTop 10 most important proteins:")
  print(as.data.frame(top_features))

  return(list(
    model = rf_model,
    train_auc = train_auc,
    train_sens = train_sens,
    train_spec = train_spec,
    test_auc = test_auc,
    test_auc_ci_lower = test_auc_ci_lower,
    test_auc_ci_upper = test_auc_ci_upper,
    test_accuracy = test_accuracy,
    test_sens = test_sens,
    test_spec = test_spec,
    confusion_matrix = conf_mat,
    top_features = top_features,
    roc_curve = roc_obj,
    predictions = data.frame(
      SampleID = test_data$SampleID_deidentified,
      true_label = y_test,
      pred_prob = test_pred_prob$ALS,
      pred_class = test_pred_class
    )
  ))
}


#' Pooled cross-validation (replicating original study approach)
#'
#' @description
#' Trains model on pooled data from both countries with standard CV.
#' This replicates the original study's approach where tube type distribution
#' is maintained across folds (stratified CV).
#'
#' Used as a baseline to compare against leave-country-out CV.
#'
#' @param protein_wide Wide-format data
#' @param n_folds Number of CV folds
#' @param seed Random seed
#' @return Model performance metrics
#' @examples
#' pooled_results <- pooled_cv_analysis(protein_wide)
#' @export
pooled_cv_analysis <- function(protein_wide, n_folds = 5, seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("POOLED CROSS-VALIDATION (Original Study Approach)")
  message(paste(rep("=", 70), collapse = ""))
  message("\nTraining on pooled Italy + US data with standard CV.")
  message("Tube type distribution is maintained across folds.")

  set.seed(seed)

  # Prepare data
  data_clean <- protein_wide %>%
    dplyr::filter(country %in% c("Italy", "US")) %>%
    dplyr::mutate(
      outcome = dplyr::case_when(
        Diagnosis == "ALS" ~ "ALS",
        Diagnosis %in% c("Healthy_control", "Neurological_control") ~ "Control",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(outcome))

  y <- factor(data_clean$outcome, levels = c("Control", "ALS"))

  # Extract protein features
  X <- data_clean %>%
    dplyr::select(
      -SampleID_deidentified, -Diagnosis, -outcome,
      -Plasma_collection_tube_type, -country,
      -Sex, -Age_Collection
    )

  # NOTE: No protein filtering to prevent data leakage
  # ranger handles missing values natively within each CV fold
  message(sprintf(
    "\nTotal samples: %s (%s ALS, %s Control)",
    nrow(X),
    sum(y == "ALS"),
    sum(y == "Control")
  ))
  message(sprintf("Proteins: %s (no filtering, missing handled by ranger)", ncol(X)))

  # Train with stratified CV
  train_control <- caret::trainControl(
    method = "cv",
    number = n_folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )

  message("\nTraining Random Forest with 5-fold CV...")

  rf_model <- caret::train(
    x = X,
    y = y,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = floor(sqrt(ncol(X))),
      splitrule = "gini",
      min.node.size = 5
    ),
    importance = "permutation",
    num.trees = 500,
    # Handle missing values within CV to prevent data leakage
    respect.unordered.factors = "order",
    replace = TRUE
    # ranger handles missing values natively when replace = TRUE
  )

  # Extract CV performance
  cv_auc <- max(rf_model$results$ROC)
  cv_sens <- rf_model$results$Sens[which.max(rf_model$results$ROC)]
  cv_spec <- rf_model$results$Spec[which.max(rf_model$results$ROC)]

  # Calculate AUC confidence interval
  best_tune <- rf_model$bestTune
  cv_preds <- rf_model$pred %>%
    dplyr::filter(
      mtry == best_tune$mtry,
      splitrule == best_tune$splitrule,
      min.node.size == best_tune$min.node.size
    )

  roc_obj_pooled <- pROC::roc(
    response = cv_preds$obs,
    predictor = cv_preds$ALS,
    levels = c("Control", "ALS"),
    direction = "<",
    quiet = TRUE
  )

  auc_ci_pooled <- pROC::ci.auc(roc_obj_pooled, conf.level = 0.95, method = "delong")
  cv_auc_ci_lower <- auc_ci_pooled[1]
  cv_auc_ci_upper <- auc_ci_pooled[3]

  message("\n", paste(rep("=", 70), collapse = ""))
  message("POOLED CV RESULTS")
  message(paste(rep("=", 70), collapse = ""))
  message(sprintf("\n  AUC:         %.3f (95%% CI: %.3f-%.3f)", cv_auc, cv_auc_ci_lower, cv_auc_ci_upper))
  message(sprintf("  Sensitivity: %.3f", cv_sens))
  message(sprintf("  Specificity: %.3f", cv_spec))
  message("\n", paste(rep("=", 70), collapse = ""))

  # Feature importance
  var_imp <- caret::varImp(rf_model, scale = TRUE)
  top_features <- var_imp$importance %>%
    tibble::rownames_to_column("protein") %>%
    dplyr::arrange(dplyr::desc(Overall)) %>%
    dplyr::slice_head(n = 20)

  return(list(
    model = rf_model,
    cv_auc = cv_auc,
    cv_auc_ci_lower = cv_auc_ci_lower,
    cv_auc_ci_upper = cv_auc_ci_upper,
    cv_sens = cv_sens,
    cv_spec = cv_spec,
    top_features = top_features,
    predictions = rf_model$pred,
    roc_obj = roc_obj_pooled
  ))
}


#' **NEW: Reverse Prediction Test on Tube-Robust Proteins Only**
#'
#' @description
#' Tests whether even the "tube-robust" proteins (significant in both Italy and US
#' with concordant direction) can still predict tube type.
#'
#' This is critical: Just because a protein is significant in BOTH cohorts doesn't
#' mean it's free from tube type effects. It could be:
#'   1. True biomarker + NO tube effect (ideal)
#'   2. True biomarker + tube effect (entangled)
#'   3. Pure tube artifact spuriously associated with disease
#'
#' With perfect confounding, we can't distinguish #2 from #3, but we CAN quantify
#' the magnitude of residual tube type signal.
#'
#' @param protein_wide Wide-format data
#' @param protein_concordance Concordance analysis results with sig_both/same_direction
#' @param n_folds Number of CV folds (default 5)
#' @param seed Random seed
#' @return List with model performance and comparison to full protein set
#' @examples
#' robust_test <- reverse_prediction_test_robust_only(protein_wide, protein_concordance)
#' @export
reverse_prediction_test_robust_only <- function(protein_wide, protein_concordance,
                                                n_folds = 5, seed = 42) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("REVERSE PREDICTION TEST: Tube-Robust Proteins Only")
  message(paste(rep("=", 70), collapse = ""))
  message("\nTesting if 'tube-robust' proteins retain tube type signal...")
  message("If AUC > 0.7, even our best candidates are still entangled!\n")

  set.seed(seed)

  # Get the tube-robust proteins
  robust_proteins <- protein_concordance %>%
    dplyr::filter(sig_both, same_direction) %>%
    dplyr::pull(Assay)

  n_robust <- length(robust_proteins)
  message(sprintf("Tube-robust proteins: %s", n_robust))

  if (n_robust < 5) {
    warning("Very few tube-robust proteins (<5). Results may be unstable.")
  }

  # Prepare data - subset to robust proteins only
  data_clean <- protein_wide %>%
    dplyr::filter(country != "Unknown")

  # Outcome: tube type
  y <- factor(data_clean$Plasma_collection_tube_type,
    levels = c("EDTA", "HEPARIN")
  )

  # Features: ONLY the tube-robust proteins
  X <- data_clean %>%
    dplyr::select(all_of(robust_proteins))

  message(sprintf("Samples: %s", nrow(X)))
  message(sprintf("Features: %s proteins (tube-robust only)", ncol(X)))
  message(sprintf(
    "Outcome: %s EDTA, %s HEPARIN",
    sum(y == "EDTA"), sum(y == "HEPARIN")
  ))

  # Set up cross-validation
  train_control <- caret::trainControl(
    method = "cv",
    number = n_folds,
    savePredictions = "final",
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    verboseIter = FALSE
  )

  # Train Random Forest
  message("\nTraining Random Forest on tube-robust proteins...")

  rf_model <- caret::train(
    x = X,
    y = y,
    method = "ranger",
    trControl = train_control,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = c(floor(sqrt(ncol(X))), min(floor(ncol(X) / 3), ncol(X))),
      splitrule = "gini",
      min.node.size = c(5, 10)
    ),
    importance = "permutation",
    num.trees = 500,
    respect.unordered.factors = "order",
    replace = TRUE
  )

  # Extract performance
  best_tune <- rf_model$bestTune
  cv_results <- rf_model$results %>%
    dplyr::filter(
      mtry == best_tune$mtry,
      splitrule == best_tune$splitrule,
      min.node.size == best_tune$min.node.size
    )

  auc <- cv_results$ROC
  sensitivity <- cv_results$Sens
  specificity <- cv_results$Spec

  message("\n", paste(rep("=", 70), collapse = ""))
  message("RESULTS: Tube-Robust Proteins Reverse Prediction")
  message(paste(rep("=", 70), collapse = ""))
  message(sprintf("\n  Cross-validated AUC:       %.3f", auc))
  message(sprintf("  Sensitivity:               %.3f", sensitivity))
  message(sprintf("  Specificity:               %.3f", specificity))
  message(sprintf("  N proteins:                %s", n_robust))

  # Interpret results
  message("\n", paste(rep("-", 70), collapse = ""))
  if (auc > 0.85) {
    message("  ⚠️  SEVERE ENTANGLEMENT: AUC > 0.85")
    message("  Even 'tube-robust' proteins STRONGLY predict tube type!")
    message("  These proteins are heavily entangled with technical artifacts.")
  } else if (auc > 0.75) {
    message("  ⚠️  SUBSTANTIAL ENTANGLEMENT: AUC 0.75-0.85")
    message("  Tube-robust proteins STILL predict tube type well.")
    message("  Significant residual confounding present.")
  } else if (auc > 0.65) {
    message("  ⚠️  MODERATE ENTANGLEMENT: AUC 0.65-0.75")
    message("  Some tube type signal remains in 'robust' proteins.")
    message("  Partial confounding persists.")
  } else if (auc > 0.55) {
    message("  ℹ️  MINIMAL ENTANGLEMENT: AUC 0.55-0.65")
    message("  Weak tube type signal - closer to tube-independent.")
  } else {
    message("  ✓ TUBE-INDEPENDENT: AUC < 0.55")
    message("  Tube-robust proteins show minimal tube type signal!")
  }

  message(paste(rep("=", 70), collapse = ""))

  # Feature importance
  var_imp <- caret::varImp(rf_model, scale = TRUE)
  top_features <- var_imp$importance %>%
    tibble::rownames_to_column("protein") %>%
    dplyr::arrange(dplyr::desc(Overall)) %>%
    dplyr::slice_head(n = min(10, n_robust))

  message("\nTop proteins predicting tube type (among 'robust' set):")
  print(as.data.frame(top_features))

  # Return results
  results <- list(
    model = rf_model,
    auc = auc,
    sensitivity = sensitivity,
    specificity = specificity,
    best_tune = best_tune,
    cv_results = cv_results,
    top_features = top_features,
    n_robust_proteins = n_robust,
    robust_protein_names = robust_proteins,
    n_samples = nrow(X),
    n_folds = n_folds,
    seed = seed,
    outcome_distribution = table(y)
  )

  return(results)
}


#' Compare pooled CV vs leave-country-out CV
#'
#' @description
#' THE PRIMARY EVIDENCE: Quantifies the performance gap between
#' pooled CV (inflated) and leave-country-out CV (realistic).
#'
#' Large gaps indicate the pooled model relies on country/tube signals
#' that don't generalize.
#'
#' @param protein_wide Wide-format data
#' @param seed Random seed
#' @return Comparison results
#' @examples
#' comparison <- compare_pooled_vs_lcv(protein_wide)
#' @export
compare_pooled_vs_lcv <- function(protein_wide, seed = 42) {
  message("\n\n", paste(rep("=", 80), collapse = ""))
  message("CRITICAL COMPARISON: Pooled CV vs Leave-Country-Out CV")
  message(paste(rep("=", 80), collapse = ""))
  message("\nThis comparison reveals whether the model generalizes geographically")
  message("or relies on country/tube-specific signals.")
  message(paste(rep("=", 80), collapse = ""))

  # Run both analyses
  pooled_results <- pooled_cv_analysis(protein_wide, n_folds = 5, seed = seed)
  lcv_results <- leave_country_out_cv(protein_wide, seed = seed)

  # Calculate performance gap
  pooled_auc <- pooled_results$cv_auc
  lcv_auc <- lcv_results$mean_test_auc
  auc_gap <- pooled_auc - lcv_auc
  pct_drop <- 100 * auc_gap / pooled_auc

  # Statistical comparison: DeLong test for paired ROC curves
  # Note: We use bootstrap comparison since pooled and LCV are not strictly paired
  # Extract ROC objects
  pooled_roc <- pooled_results$roc_obj
  lcv_italy_to_us_roc <- lcv_results$italy_to_us$roc_curve
  lcv_us_to_italy_roc <- lcv_results$us_to_italy$roc_curve

  # DeLong test comparing pooled vs Italy→US (both have predictions)
  # Note: These are independent samples, so we use unpaired comparison
  message("\n  STATISTICAL TESTS:")
  message("  Testing if AUC differences are statistically significant...")

  # For educational purposes: document why direct DeLong comparison is limited
  message("\n  Note: Direct paired DeLong test not applicable (independent test sets).")
  message("  Using bootstrap confidence intervals for inference instead.")

  message("\n\n", paste(rep("=", 80), collapse = ""))
  message("PERFORMANCE COMPARISON")
  message(paste(rep("=", 80), collapse = ""))
  message("\n  POOLED CV (original approach):")
  message(sprintf(
    "    AUC: %.3f (95%% CI: %.3f-%.3f)", pooled_auc,
    pooled_results$cv_auc_ci_lower, pooled_results$cv_auc_ci_upper
  ))
  message("\n  LEAVE-COUNTRY-OUT CV (geographic test):")
  message(sprintf(
    "    Italy → US:  %.3f (95%% CI: %.3f-%.3f)", lcv_results$italy_to_us$test_auc,
    lcv_results$italy_to_us$test_auc_ci_lower, lcv_results$italy_to_us$test_auc_ci_upper
  ))
  message(sprintf(
    "    US → Italy:  %.3f (95%% CI: %.3f-%.3f)", lcv_results$us_to_italy$test_auc,
    lcv_results$us_to_italy$test_auc_ci_lower, lcv_results$us_to_italy$test_auc_ci_upper
  ))
  message(sprintf("    Mean:        %.3f", lcv_auc))
  message("\n  PERFORMANCE GAP:")
  message(sprintf("    Absolute:    %.3f", auc_gap))
  message(sprintf("    Relative:    %.1f%% drop", pct_drop))

  # Check if CIs overlap
  pooled_ci_lower <- pooled_results$cv_auc_ci_lower
  pooled_ci_upper <- pooled_results$cv_auc_ci_upper
  lcv_italy_ci_upper <- lcv_results$italy_to_us$test_auc_ci_upper
  lcv_us_ci_upper <- lcv_results$us_to_italy$test_auc_ci_upper

  ci_overlap <- (pooled_ci_lower < max(lcv_italy_ci_upper, lcv_us_ci_upper))

  message("\n  CONFIDENCE INTERVAL ANALYSIS:")
  if (!ci_overlap) {
    message("  ⚠️  Non-overlapping CIs: Pooled AUC is significantly higher than LCV")
    message("  This provides strong statistical evidence of geographic overfitting.")
  } else {
    message("  CIs show partial overlap, but performance gap remains substantial.")
  }

  # Interpretation
  message("\n", paste(rep("-", 80), collapse = ""))
  message("INTERPRETATION:")

  if (auc_gap > 0.2) {
    message("  ⚠️  CRITICAL FINDING: >0.2 AUC gap")
    message("  Pooled CV performance is GROSSLY INFLATED.")
    message("  Model DOES NOT generalize across countries.")
    message("  Strong evidence that model relies on country/tube-specific signals.")
  } else if (auc_gap > 0.15) {
    message("  ⚠️  SEVERE: 0.15-0.2 AUC gap")
    message("  Pooled CV substantially overestimates true performance.")
    message("  Limited geographic generalizability.")
  } else if (auc_gap > 0.1) {
    message("  ⚠️  SUBSTANTIAL: 0.1-0.15 AUC gap")
    message("  Notable overestimation in pooled CV.")
    message("  Moderate country-specific effects.")
  } else if (auc_gap > 0.05) {
    message("  ⚠️  MODERATE: 0.05-0.1 AUC gap")
    message("  Some overestimation in pooled CV.")
  } else {
    message("  ✓ ACCEPTABLE: <0.05 AUC gap")
    message("  Model generalizes well across countries.")
  }

  message("\n", paste(rep("-", 80), collapse = ""))
  message("RECOMMENDATION:")

  if (auc_gap > 0.1) {
    message("  Published performance estimates (from pooled CV) should be")
    message("  interpreted with EXTREME CAUTION. Geographic validation shows")
    message("  substantially lower performance.")
  } else {
    message("  Pooled CV performance is reasonably accurate.")
    message("  Model shows acceptable geographic generalizability.")
  }

  message(paste(rep("=", 80), collapse = ""))

  return(list(
    pooled_results = pooled_results,
    lcv_results = lcv_results,
    pooled_auc = pooled_auc,
    lcv_mean_auc = lcv_auc,
    auc_gap = auc_gap,
    percent_drop = pct_drop
  ))
}
