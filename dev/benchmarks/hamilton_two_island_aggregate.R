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
A$both <- as.logical(A$both)

# ---- PART A: completeness. Goal = fraction capturing BOTH islands >= 0.90 ----
# NB use split()+data.frame, NOT aggregate(cbind(...) ~ ., scalar-FUN): the latter
# returns vector (not matrix) response columns, so the earlier `agg$both[,1]`
# indexing crashed with "incorrect number of dimensions".
cat("\n================ PART A: two-island recovery ================\n")
cat("frac(both islands), mean island2/16, replicates, by variant x budget:\n")
spA <- split(A, list(A$variant, A$budget), drop = TRUE)
tabA <- do.call(rbind, lapply(spA, function(d) data.frame(
  variant = d$variant[1], budget = d$budget[1], n = nrow(d),
  frac_both = round(mean(d$both), 3),
  mean_isl2 = round(mean(d$island2), 1),
  mean_main = round(mean(d$main), 1),
  mean_reps = round(mean(d$replicates), 1),
  mean_s    = round(mean(d$elapsed), 1))))
tabA <- tabA[order(tabA$budget, -tabA$frac_both), ]
print(tabA, row.names = FALSE)
cat("\nPASS = frac_both >= 0.90. Winner = cheapest (variant,budget) that clears it.\n")

# Throughput hypothesis: does recovery track completed replicates?
cat("\nthroughput check -- frac_both vs mean_reps (sorted by reps):\n")
print(tabA[order(tabA$mean_reps), c("variant", "budget", "mean_reps", "frac_both")],
      row.names = FALSE)

# ---- PART B: score/throughput regression + wagnerStarts merge ----------------
cat("\n================ PART B: score regression / merge ================\n")
cat("median(score - target) by dataset x budget x variant (0 = at target, + = worse):\n")
mb <- aggregate(over ~ dataset + budget + variant, B, median)
w <- reshape(mb, idvar = c("dataset", "budget"), timevar = "variant",
             direction = "wide")
names(w) <- sub("over\\.", "", names(w))
print(w[order(w$dataset, w$budget), ], row.names = FALSE)

cat("\nmedian replicates + elapsed by variant x budget (throughput cost):\n")
cost <- aggregate(cbind(replicates, elapsed) ~ variant + budget, B,
                  function(x) round(median(x), 1))
print(cost[order(cost$budget, cost$variant), ], row.names = FALSE)

# MERGE TEST: base (wagnerStarts=3) vs ws5 (=intensive). Median AND worst-case,
# so an occasional +1 trade-off (the documented Zanol2014/Giles2015 concern)
# cannot hide behind the median.
cat("\nMERGE TEST (thorough vs intensive): base (ws=3) vs ws5 (=intensive)\n")
mtB <- B[B$variant %in% c("base", "ws5"), ]
mtab <- do.call(rbind, lapply(
  split(mtB, list(mtB$dataset, mtB$budget, mtB$variant), drop = TRUE),
  function(d) data.frame(dataset = d$dataset[1], budget = d$budget[1],
                         variant = d$variant[1], med_over = median(d$over),
                         worst_over = max(d$over))))
mw <- reshape(mtab[, c("dataset", "budget", "variant", "med_over")],
              idvar = c("dataset", "budget"), timevar = "variant", direction = "wide")
names(mw) <- sub("med_over\\.", "", names(mw))
mw$ws5_minus_base <- mw$ws5 - mw$base
print(mw[order(mw$dataset, mw$budget), ], row.names = FALSE)
cat("Merge safe if ws5_minus_base <= 0 everywhere (esp. Zanol2014/Giles2015).\n")
