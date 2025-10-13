# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(tidyverse)
library(data.table)
library(caret)
library(ranger)
library(limma)
library(pROC)

# Source all R functions
source_files <- list.files(
  path = "../../R",
  pattern = "\\.R$",
  full.names = TRUE
)

for (file in source_files) {
  source(file)
}

test_check("chia_etal_2025_als")
