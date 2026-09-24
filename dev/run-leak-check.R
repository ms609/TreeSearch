#!/usr/bin/env Rscript
# Runs each test file in a fresh process and captures stdout+stderr
suppressMessages({
  library(testthat)
  library(TreeSearch)
  library(TreeTools)
})
Sys.setenv(NOT_CRAN = "true")

args <- commandArgs(trailingOnly = TRUE)
for (f in args) {
  cat("===BEGIN===", f, "\n", sep = "")
  try(test_file(paste0("tests/testthat/", f),
                reporter = "minimal",
                env = new.env(parent = asNamespace("TreeSearch"))),
      silent = FALSE)
  cat("\n===END===", f, "\n", sep = "")
}
