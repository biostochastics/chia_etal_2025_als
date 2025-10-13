# Test Suite - ALS Confounding Investigation

**Status:** 14 tests across 3 test files
**Coverage Target:** >80% line coverage
**Date:** 2025-10-13

---

## Quick Start

### Run All Tests

```r
# From R console in project root
testthat::test_dir("tests/testthat")
```

### Run Individual Test Files

```r
# Data integrity tests (critical confounding verification)
testthat::test_file("tests/testthat/test-data_integrity.R")

# ML evaluation tests (data leakage prevention, missing value handling)
testthat::test_file("tests/testthat/test-ml_evaluation.R")

# Differential analysis tests (limma correctness, meta-analysis)
testthat::test_file("tests/testthat/test-differential_analysis.R")
```

---

## Test Coverage Summary

### test-data_integrity.R (6 tests)

**Purpose:** Validates critical assumptions about confounding structure

✅ **Test 1: Perfect confounding verification**
- Ensures 100% of neurological controls are EDTA
- CRITICAL for investigation validity

✅ **Test 2: Country-tube consistency**
- Italy = HEPARIN, US = EDTA
- Validates country reconstruction logic

✅ **Test 3: ALS tube distribution**
- Checks 80-90% HEPARIN in ALS samples
- Confirms confounding structure

✅ **Test 4: Plate-country mapping**
- Plates 1-4 → US, Plates 5-9 → Italy
- Tests reconstruction algorithm

✅ **Test 5: No unknown countries**
- All samples mapped to Italy or US
- No ambiguous classifications

✅ **Test 6: Sample ID uniqueness**
- No duplicate sample IDs within proteins
- Data integrity check

**Note:** Some tests marked `skip()` - require full dataset. Run manually after loading data.

### test-ml_evaluation.R (4 tests)

**Purpose:** Validates ML methods and prevents data leakage

✅ **Test 1: Reproducibility with fixed seed**
- Template for manual testing
- Verifies deterministic results

✅ **Test 2: Leave-country-out CV correctness**
- Confirms independent feature selection
- Validates no data leakage in LCV

✅ **Test 3: AUC calculations**
- Tests with synthetic data
- Perfect prediction → AUC = 1.0
- Random prediction → AUC ≈ 0.5

✅ **Test 4: ranger missing value handling**
- Validates `na.action = na.pass` works
- Tests with simulated missing data
- **CRITICAL**: Confirms fix for data leakage

✅ **Test 5: Model comparison edge cases**
- Small AUC differences
- Percentage change calculations

**Note:** Tests 1-2 require full dataset. Tests 3-5 run immediately with synthetic data.

### test-differential_analysis.R (4 tests)

**Purpose:** Validates limma implementation and meta-analysis

✅ **Test 1: limma results structure**
- Required columns present
- FDR control working
- P-values in valid range [0, 1]

✅ **Test 2: Meta-analysis Fisher's method**
- Directional concordance detection
- Combined p-values calculated correctly
- Chi-square statistic validation

✅ **Test 3: Set operation correctness**
- Overlap calculations (Italy ∩ US)
- Pooled-only proteins (setdiff)
- Replication percentages

✅ **Test 4: Significance thresholds**
- FDR < 0.05 applied correctly
- |logFC| > 0.5 filtering works
- Combined threshold logic

**Note:** Test 1 requires full dataset. Tests 2-4 run immediately with mock data.

---

## Running Tests with Full Dataset

Some tests are marked `skip()` and require loading the full dataset. To run these manually:

```r
# Load targets pipeline
targets::tar_load(c(
  protein_wide,
  protein_long,
  sample_metadata,
  differential_italy,
  differential_us
))

# Uncomment skip() lines in test files
# Then run:
testthat::test_file("tests/testthat/test-data_integrity.R")
testthat::test_file("tests/testthat/test-ml_evaluation.R")
testthat::test_file("tests/testthat/test-differential_analysis.R")
```

---

## Test Philosophy

### What We Test

1. **Critical Assumptions:**
   - Perfect confounding (100% neuro controls = EDTA)
   - Country reconstruction accuracy
   - Sample distribution patterns

2. **Methodological Correctness:**
   - No data leakage (protein filtering in CV)
   - Missing value handling
   - Cross-validation stratification

3. **Statistical Validity:**
   - AUC calculations
   - Fisher's combined p-values
   - Set operations for overlaps

4. **Edge Cases:**
   - Small AUC differences
   - Missing data patterns
   - Threshold boundary conditions

### What We Don't Test

- **Specific AUC values:** These will change after fixes (0.999 → 0.95-0.995)
- **Exact protein counts:** Biology-dependent, not deterministic
- **Plot aesthetics:** Visual validation done manually

---

## Expected Test Results

### Before Pipeline Rerun (Current State)

Many tests will **skip** because they require full dataset that hasn't been loaded yet.

**Synthetic data tests should PASS:**
- ✅ AUC calculation correctness
- ✅ ranger missing value handling
- ✅ Meta-analysis Fisher's method
- ✅ Set operation logic
- ✅ Significance threshold application

### After Pipeline Rerun (`targets::tar_make()`)

All tests should be runnable after uncommenting `skip()` directives.

**Expected results:**
- ✅ All data integrity tests PASS (confounding structure confirmed)
- ✅ All ML evaluation tests PASS (no data leakage, correct methods)
- ✅ All differential analysis tests PASS (valid statistics)

**If any test FAILS:**
1. Check that pipeline completed successfully
2. Verify data was loaded correctly with `targets::tar_load()`
3. Inspect test failure message for specific issue
4. Re-run failed test file individually for detailed output

---

## Continuous Integration

### Pre-commit Checks

```r
# Before committing code changes
devtools::test()           # Run all tests
covr::package_coverage()   # Check coverage
devtools::check()          # R CMD check
```

### Test-Driven Development Workflow

```r
# 1. Write test first
testthat::test_file("tests/testthat/test-new_function.R")

# 2. Write function to pass test
devtools::load_all()

# 3. Iterate until passing
testthat::test_file("tests/testthat/test-new_function.R")
```

---

## Adding New Tests

### Template Structure

```r
# tests/testthat/test-new_module.R

test_that("descriptive test name", {
  # Arrange: Set up test data
  test_data <- tibble::tibble(x = 1:10, y = 11:20)

  # Act: Call function being tested
  result <- my_function(test_data)

  # Assert: Check expectations
  expect_equal(result$x, 1:10)
  expect_true(all(result$y > result$x))
})
```

### Best Practices

1. **One concept per test:** Each test should validate one specific behavior
2. **Descriptive names:** `test_that("reverses vector correctly")` not `test_that("test1")`
3. **Use synthetic data:** Create minimal examples rather than loading full dataset
4. **Test edge cases:** Empty data, NA values, extreme values
5. **Document assumptions:** Use `info` parameter to explain what should happen

---

## Test Dependencies

Tests require the following packages:

```r
# Testing framework
testthat

# Data manipulation
tidyverse, data.table

# ML methods
caret, ranger, pROC

# Statistics
limma

# Utilities
DescTools (Cramér's V), broom
```

All dependencies are managed by `renv` and will be automatically installed with `renv::restore()`.

---

## Coverage Targets

**Current Status:** Test suite created (14 tests)
**Line Coverage Target:** >80% (aim for 90%)

### Priority for Additional Tests

1. **High Priority:**
   - Integration tests (full pipeline end-to-end)
   - Edge cases in differential_analysis.R functions
   - Visualization function structure checks

2. **Medium Priority:**
   - Additional ML evaluation scenarios
   - More complex confounding patterns
   - Performance benchmarks

3. **Low Priority:**
   - Helper function utilities
   - Data export functions
   - Report rendering checks

---

## Troubleshooting

### Error: "cannot find function 'my_function'"

**Solution:** Source R files before testing

```r
devtools::load_all()
# Or manually:
source("R/04_ml_evaluation.R")
```

### Error: "object 'protein_wide' not found"

**Solution:** Load targets pipeline data

```r
targets::tar_load(protein_wide)
```

### Tests skip unexpectedly

**Cause:** Tests marked with `skip()` for full dataset

**Solution:** Uncomment `skip()` line after loading data

---

## Contact

**Questions about tests?**
Contact: Sergey A. Kornilov, PhD
Email: sergey.kornilov@biostochastics.com

---

**Version:** 1.0
**Last Updated:** 2025-10-13
**Status:** Test suite complete, ready for pipeline rerun validation
