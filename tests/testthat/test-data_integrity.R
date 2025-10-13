# Test Data Integrity and Country Reconstruction
# Critical tests for perfect confounding verification

test_that("perfect confounding is verified - 100% neurological controls are EDTA", {
  skip_on_cran()
  skip_if_not(file.exists("../../original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt"))

  # Load minimal data for testing
  data <- data.table::fread(
    "../../original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt",
    nrows = 100000  # Sample for speed
  )

  data_with_country <- add_country_labels(data)
  sample_meta <- create_sample_metadata(data_with_country)

  # CRITICAL TEST: 100% of neurological controls must be EDTA
  neuro <- sample_meta %>%
    dplyr::filter(Diagnosis == "Neurological_control")

  if (nrow(neuro) > 0) {
    expect_true(all(neuro$Plasma_collection_tube_type == "EDTA"),
                info = "ALL neurological controls must use EDTA tubes")
    expect_equal(nrow(neuro), sum(neuro$Plasma_collection_tube_type == "EDTA"),
                 info = "Count mismatch in neurological EDTA tubes")
  }
})

test_that("country labels are consistent with tube type", {
  skip_on_cran()
  skip_if_not(file.exists("../../original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt"))

  data <- data.table::fread(
    "../../original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt",
    nrows = 100000
  )

  data_with_country <- add_country_labels(data)

  # Italy samples must all be HEPARIN
  italy <- data_with_country %>%
    dplyr::filter(country == "Italy")

  if (nrow(italy) > 0) {
    expect_true(all(italy$Plasma_collection_tube_type == "HEPARIN"),
                info = "All Italy samples must use HEPARIN tubes")
  }

  # US samples must all be EDTA
  us <- data_with_country %>%
    dplyr::filter(country == "US")

  if (nrow(us) > 0) {
    expect_true(all(us$Plasma_collection_tube_type == "EDTA"),
                info = "All US samples must use EDTA tubes")
  }
})

test_that("ALS samples are predominantly HEPARIN (>80%)", {
  skip_on_cran()
  skip_if_not(file.exists("../../original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt"))

  data <- data.table::fread(
    "../../original/Chia_et_al_OLINK_bridge_normalized_data_merged_deidentified_updated.txt",
    nrows = 1000000
  )

  data_with_country <- add_country_labels(data)
  sample_meta <- create_sample_metadata(data_with_country)

  als <- sample_meta %>% dplyr::filter(Diagnosis == "ALS")

  if (nrow(als) > 0) {
    pct_heparin <- 100 * sum(als$Plasma_collection_tube_type == "HEPARIN") / nrow(als)

    expect_gt(pct_heparin, 80, info = "At least 80% of ALS should be HEPARIN")
    expect_lt(pct_heparin, 90, info = "At most 90% of ALS should be HEPARIN")
  }
})

test_that("plate numbers map correctly to countries", {
  # Test the country reconstruction logic
  test_data <- tibble::tibble(
    SampleID_deidentified = paste0("S", 1:9),
    Olink_Plate_No = paste0("Plate_", 1:9),
    Plasma_collection_tube_type = c(rep("EDTA", 4), rep("HEPARIN", 5)),
    Diagnosis = "ALS"
  ) %>%
    data.table::as.data.table()

  result <- add_country_labels(test_data)

  # Plates 1-4 should be US
  expect_equal(result$country[1:4], rep("US", 4))

  # Plates 5-9 should be Italy
  expect_equal(result$country[5:9], rep("Italy", 5))
})

test_that("no Unknown country labels are created with valid data", {
  test_data <- tibble::tibble(
    SampleID_deidentified = c("S1", "S2"),
    Olink_Plate_No = c("Plate_1", "Plate_5"),
    Plasma_collection_tube_type = c("EDTA", "HEPARIN"),
    Diagnosis = "ALS"
  ) %>%
    data.table::as.data.table()

  result <- add_country_labels(test_data)

  expect_false(any(result$country == "Unknown"),
               info = "Valid plate-tube combinations should not produce Unknown country")
})
