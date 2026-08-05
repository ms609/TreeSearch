# Run the features test file under LocationReporter so the LAST "Start" line
# printed before an abnormal exit pinpoints the crashing expectation.
# Lib chosen via TS_LIB env (.agent-l1 = lever-1, .agent-base = unmodified).
lib <- Sys.getenv("TS_LIB", ".agent-l1")
suppressMessages(library(TreeSearch, lib.loc = lib))
library(testthat)
cat(sprintf(">>> testing lib=%s\n", lib)); flush.console()
test_file("tests/testthat/test-MaximizeParsimony-features.R",
          reporter = LocationReporter$new())
cat(">>> FILE COMPLETED WITHOUT CRASH\n")
