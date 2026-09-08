# Regression check: default path must be intact after the directional-vroot
# merge + the opt-in unrooted reroot mechanism (which defaults off).
suppressMessages({
  library(testthat)
  library(TreeSearch, lib.loc = ".agent-tbr")
})
files <- list.files("tests/testthat", pattern = "^test-", full.names = TRUE)
keep <- grepl("ts-tbr|ts-sector|ts-driven|ts-ratchet|ts-drift|ts-spr|SPR|wagner|MaximizeParsimony|SearchControl",
              files, ignore.case = TRUE)
files <- files[keep]
fail <- 0L
for (f in files) {
  cat("---", basename(f), "---\n")
  r <- tryCatch(as.data.frame(test_file(f, reporter = "silent")),
                error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) {
    nf <- sum(r$failed); ne <- sum(r$error %in% TRUE)
    cat(sprintf("   passed=%d failed=%d error=%d\n", sum(r$passed), nf, ne))
    fail <- fail + nf + ne
  } else fail <- fail + 1L
}
cat(sprintf("\nTOTAL failures/errors: %d\n", fail))
