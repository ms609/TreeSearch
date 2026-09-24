# Run the search-related testthat files against the lever-1 lib (.agent-l1).
library(TreeSearch, lib.loc = ".agent-l1")
library(testthat)
res <- test_dir(
  "tests/testthat",
  filter = "fitch|tbr|wagner|sector|ratchet|drift|MaximizeParsimony|CustomSearch|SearchControl",
  reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n=== TOTALS: pass=%d fail=%d warn=%d skip=%d ===\n",
            sum(df$passed), sum(df$failed), sum(df$warning), sum(df$skipped)))
if (sum(df$failed) > 0) {
  print(df[df$failed > 0, c("file", "test", "failed")])
  quit(status = 1)
}
