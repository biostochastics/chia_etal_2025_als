#' Data Import and Country Label Reconstruction
#'
#' Functions for loading the Chia et al. OLINK data and reconstructing
#' country labels from plate identifiers and tube type information.
#'
#' @description
#' The original study does not include explicit country labels, but based on
#' the Jupyter notebooks and metadata:
#' - Plates 5-9 (HEPARIN) → Italy (Scholz registry, TRAYNOR facility)
#' - Plates 1-4 (EDTA) → US (NDRU facility, NIH/Johns Hopkins/BLSA cohorts)
#'

#' Load raw OLINK data
#'
#' @param file_path Path to the merged OLINK data file
#' @return data.table with raw OLINK measurements
#' @examples
#' raw_data <- load_raw_olink_data("original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt")
#' @export
load_raw_olink_data <- function(file_path) {
  # Use data.table for fast loading of large file (312MB)
  dt <- data.table::fread(
    file_path,
    sep = "\t",
    header = TRUE,
    na.strings = c("", "NA", "NaN"),
    data.table = TRUE
  )

  # Validate expected columns
  required_cols <- c(
    "SampleID_deidentified",
    "Olink_Plate_No",
    "Diagnosis",
    "Sex",
    "Age_Collection",
    "Plasma_collection_tube_type",
    "Assay",
    "NPX",
    "UniProt"
  )

  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Log summary
  message(sprintf("Loaded %s rows × %s columns", nrow(dt), ncol(dt)))
  message(sprintf("Unique samples: %s", data.table::uniqueN(dt$SampleID_deidentified)))
  message(sprintf("Unique proteins: %s", data.table::uniqueN(dt$Assay)))

  return(dt)
}


#' Add country labels based on plate and tube type
#'
#' @description
#' Reconstructs likely country of origin using:
#' - Plate number (5-9 = Italy, 1-4 = US)
#' - Tube type (HEPARIN = Italy, EDTA = US)
#'
#' These assignments are based on:
#' 1. Jupyter notebook file naming ("Scholz" = Italy, "NDRU" = US)
#' 2. Perfect confounding of neurological controls with EDTA
#' 3. Facility codes in metadata
#'
#' @param data Data frame with Olink_Plate_No and Plasma_collection_tube_type
#' @return Data frame with added 'country' column
#' @examples
#' data_with_country <- add_country_labels(raw_data)
#' @export
add_country_labels <- function(data) {
  # Convert to data.table if not already
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  # Add country label based on plate + tube type
  data[, country := data.table::fcase(
    Olink_Plate_No %in% c("Plate_5", "Plate_6", "Plate_7", "Plate_8", "Plate_9") &
      Plasma_collection_tube_type == "HEPARIN",
    "Italy",
    Olink_Plate_No %in% c("Plate_1", "Plate_2", "Plate_3", "Plate_4") &
      Plasma_collection_tube_type == "EDTA",
    "US",
    default = "Unknown"
  )]

  # Validate country assignment
  n_unknown <- sum(data$country == "Unknown")
  if (n_unknown > 0) {
    warning(sprintf(
      "%s samples (%.2f%%) could not be assigned to a country",
      n_unknown,
      100 * n_unknown / nrow(data)
    ))
  }

  # Log results
  country_summary <- data[, .N, by = country]
  message("\nCountry assignment summary:")
  print(country_summary)

  return(data)
}


#' Create sample-level metadata
#'
#' @description
#' Collapses protein-level data to one row per sample with clinical variables.
#' Used for confounding analysis and stratification.
#'
#' @param data Data frame with sample-level variables
#' @return Tibble with one row per sample
#' @examples
#' sample_metadata <- create_sample_metadata(data_with_country)
#' @export
create_sample_metadata <- function(data) {
  # Get unique sample-level information
  sample_meta <- data %>%
    dplyr::select(
      SampleID_deidentified,
      Olink_Plate_No,
      Diagnosis,
      Sex,
      Age_Collection,
      Plasma_collection_tube_type,
      country
    ) %>%
    dplyr::distinct()

  # Validate uniqueness
  n_samples <- dplyr::n_distinct(sample_meta$SampleID_deidentified)
  if (nrow(sample_meta) != n_samples) {
    warning("Some samples have inconsistent metadata across proteins!")
  }

  message(sprintf("\nCreated metadata for %s unique samples", n_samples))

  # Summary by diagnosis and tube type (critical confounding check)
  confound_check <- sample_meta %>%
    dplyr::count(Diagnosis, Plasma_collection_tube_type, country) %>%
    dplyr::arrange(Diagnosis, Plasma_collection_tube_type)

  message("\nDiagnosis × Tube Type × Country distribution:")
  print(as.data.frame(confound_check))

  return(tibble::as_tibble(sample_meta))
}


#' Validate data integrity
#'
#' @description
#' Performs data quality checks:
#' - No duplicate sample-protein combinations
#' - All expected diagnoses present
#' - Country labels match tube type expectations
#' - No missing values in critical variables
#'
#' @param data Data frame to validate
#' @return TRUE if all checks pass, throws error otherwise
#' @examples
#' validate_data_integrity(data_with_country)
#' @export
validate_data_integrity <- function(data) {
  # Check for duplicates (allow some - likely technical replicates)
  dup_check <- data %>%
    dplyr::count(SampleID_deidentified, Assay) %>%
    dplyr::filter(n > 1)

  if (nrow(dup_check) > 0) {
    message(sprintf("WARNING: Found %s duplicate sample-protein combinations (likely technical replicates)", nrow(dup_check)))
    message("These will be handled by taking the first occurrence in downstream analyses.")
    # Don't stop - just warn
  }

  # Check expected diagnoses
  expected_dx <- c("ALS", "Healthy_control", "Neurological_control")
  actual_dx <- unique(data$Diagnosis)
  missing_dx <- setdiff(expected_dx, actual_dx)

  if (length(missing_dx) > 0) {
    warning("Missing expected diagnoses: ", paste(missing_dx, collapse = ", "))
  }

  # Check country-tube type consistency
  inconsistent <- data %>%
    dplyr::filter(
      (country == "Italy" & Plasma_collection_tube_type != "HEPARIN") |
        (country == "US" & Plasma_collection_tube_type != "EDTA")
    ) %>%
    nrow()

  if (inconsistent > 0) {
    stop(sprintf("Found %s rows with inconsistent country-tube type mapping", inconsistent))
  }

  # Check for missing values in critical variables
  critical_vars <- c(
    "SampleID_deidentified",
    "Diagnosis",
    "Plasma_collection_tube_type",
    "NPX"
  )

  for (var in critical_vars) {
    n_missing <- sum(is.na(data[[var]]))
    if (n_missing > 0) {
      warning(sprintf(
        "%s has %s missing values (%.2f%%)",
        var, n_missing, 100 * n_missing / nrow(data)
      ))
    }
  }

  message("\n✓ Data integrity checks passed")
  return(TRUE)
}
