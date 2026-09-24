# Run the IW-specific testthat files against a chosen lib (TS_LIB, default
# .agent-tbr), loading helper-*.R the way testthat does so the tests actually
# execute (test_file alone does not auto-source helpers).
lib <- normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"), winslash = "/")
library(TreeSearch, lib.loc = lib)
library(testthat)
helpers <- list.files("tests/testthat", pattern = "^helper", full.names = TRUE)
for (h in helpers) sys.source(h, envir = globalenv())

files <- c("test-ts-iw.R", "test-iw-scoring.R", "test-ts-nni-iw-rescore.R",
           "test-ts-xpiwe.R", "test-ts-iw-profile-red10.R")
tp <- 0; tf <- 0; te <- 0; ts <- 0
for (f in files) {
  p <- file.path("tests/testthat", f)
  if (!file.exists(p)) { cat(sprintf("-- %-32s MISSING\n", f)); next }
  r <- tryCatch(as.data.frame(test_file(p, reporter = "silent")),
                error = function(e) { cat("FILE ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  pa <- sum(r$passed); fl <- sum(r$failed)
  er <- sum(as.integer(r$error), na.rm = TRUE)
  sk <- if ("skipped" %in% names(r)) sum(as.integer(r$skipped), na.rm = TRUE) else 0L
  tp <- tp + pa; tf <- tf + fl; te <- te + er; ts <- ts + sk
  cat(sprintf("-- %-32s passed=%d failed=%d error=%d skipped=%d\n", f, pa, fl, er, sk))
}
cat(sprintf("\nIW TESTS TOTAL: passed=%d failed=%d error=%d skipped=%d\n", tp, tf, te, ts))
