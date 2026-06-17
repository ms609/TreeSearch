setwd("C:/Users/pjjg18/GitHub/TreeSearch")
source("dev/benchmarks/bench_tnt_settings.R")
results <- tnt_settings_full()
message("Survey complete. ", nrow(results), " rows collected.")
