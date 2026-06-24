# Aggregate the two-island sweep shards into the two decision tables.
# Usage (on Hamilton, after the array finishes):
#   Rscript dev/benchmarks/hamilton_two_island_aggregate.R [results_dir]
# Default results_dir: /nobackup/$USER/TreeSearch/two_island_results
args <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(args)) args[[1]] else
        file.path("/nobackup", Sys.getenv("USER"), "TreeSearch/two_island_results")
files <- list.files(dir, "^two_island_sweep_shard.*\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
df <- do.call(rbind, lapply(files, read.csv))
write.csv(df, file.path(dir, "two_island_sweep_ALL.csv"), row.names = FALSE)
cat(sprintf("combined %d shards -> %d rows\n", length(files), nrow(df)))

A <- df[df$part == "A", ]
B <- df[df$part == "B", ]

# ---- PART A: completeness. Goal = fraction capturing BOTH islands >= 0.90 ----
cat("\n================ PART A: two-island recovery ================\n")
cat("frac(both islands) and mean island2 / replicates, by variant x budget:\n")
agg <- aggregate(cbind(both, island2, replicates, elapsed) ~ variant + budget, A,
                 function(x) round(mean(x), 2))
agg <- agg[order(agg$budget, -agg$both[,1]), ]
print(agg, row.names = FALSE)
cat("\nPASS threshold = frac(both) >= 0.90. Winning (variant,budget) = cheapest that clears it.\n")

# Throughput hypothesis: does recovery track completed replicates?
cat("\nthroughput check -- frac(both) vs mean replicates (per variant x budget):\n")
A$bothn <- as.integer(A$both)
thr <- aggregate(cbind(bothn, replicates) ~ variant + budget, A, mean)
print(thr[order(thr$replicates), ], row.names = FALSE)

# ---- PART B: score/throughput regression + wagnerStarts merge ----------------
cat("\n================ PART B: score regression / merge ================\n")
cat("median(score - target) by dataset x budget x variant (negative = better):\n")
mb <- aggregate(over ~ dataset + budget + variant, B, median)
w <- reshape(mb, idvar = c("dataset", "budget"), timevar = "variant",
             direction = "wide")
names(w) <- sub("over\\.", "", names(w))
print(w, row.names = FALSE)

cat("\nmedian replicates + elapsed by variant x budget (throughput cost):\n")
cost <- aggregate(cbind(replicates, elapsed) ~ variant + budget, B,
                  function(x) round(median(x), 1))
print(cost[order(cost$budget, cost$variant), ], row.names = FALSE)

cat("\nMERGE TEST (thorough vs intensive): compare variant 'base' (wagnerStarts=3)\n")
cat("vs 'ws5' (=intensive). Merge is safe if ws5 does not regress median over-target\n")
cat("at matched budget on any dataset (esp. Zanol2014/Giles2015, the documented +1).\n")
