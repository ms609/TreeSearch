# Correctness gate for the sectorial micro-levers + TBR getenv hoists.
# Run against the hot-swapped .agent-sect lib (absolute lib.loc per testthat CWD gotcha).
lib <- normalizePath(".agent-sect")
suppressMessages(library(TreeSearch, lib.loc = lib))
suppressMessages(library(TreeTools))
suppressMessages(library(testthat))
Sys.setenv(NOT_CRAN = "true")

# testthat auto-sources helper-*.R only via test_dir; source them manually.
helpers <- list.files("tests/testthat", pattern = "^helper.*[.]R$", full.names = TRUE)
for (h in helpers) sys.source(h, envir = globalenv())

files <- c("test-ts-sector.R", "test-ts-sector-resolve.R", "test-ts-conflict-sector.R",
           "test-ts-tbr-search.R", "test-ts-tbr-dirty-rescore.R", "test-ts-tbr-symmetry.R",
           "test-ts-ratchet-search.R", "test-ts-drift-search.R")
for (f in files) {
  p <- file.path("tests/testthat", f)
  if (!file.exists(p)) { cat("(skip missing", f, ")\n"); next }
  cat("\n===========", f, "===========\n")
  tryCatch(
    test_file(p, reporter = "summary"),
    error = function(e) cat("FILE ERROR:", conditionMessage(e), "\n")
  )
}
