# T-300 NA wall-clock A/B
# Usage:  Rscript dev/profiling/t300_na_bench.R <lib_dir> <reps> <label>
# Runs the same Zhu2013 ratchet workload that produced the 18.2% VTune share
# on c504ea87 (full_rescore accept path), now under HEAD's dirty-set NA path.

args <- commandArgs(trailingOnly = TRUE)
lib_dir <- args[1]
reps    <- as.integer(if (length(args) >= 2) args[2] else 5L)
label   <- if (length(args) >= 3) args[3] else basename(lib_dir)

library(TreeSearch, lib.loc = lib_dir)

dataset  <- inapplicable.phyData[["Zhu2013"]]
at       <- attributes(dataset)
contrast <- at$contrast
tip_data <- matrix(unlist(dataset, use.names = FALSE),
                   nrow = length(dataset), byrow = TRUE)
weight   <- TreeSearch:::.ScaleWeight(at$weight)
levels   <- at$levels

set.seed(5813)
starting_edge <- ape::rtree(length(dataset), tip.label = names(dataset),
                             rooted = FALSE)
starting_edge <- ape::root(starting_edge, 1L, resolve.root = TRUE)[["edge"]]
stopifnot(starting_edge[1L, 1L] > length(dataset))

run_once <- function() {
  t0 <- proc.time()
  last <- NULL
  for (rep in seq_len(12L)) {
    set.seed(rep)
    last <- TreeSearch:::ts_ratchet_search(
      edge        = starting_edge,
      contrast    = contrast,
      tip_data    = tip_data,
      weight      = weight,
      levels      = levels,
      nCycles     = 12L,
      perturbProb = 0.04,
      maxHits     = 1L
    )
  }
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  list(elapsed = elapsed, score = last$score)
}

# Warmup (touch caches, JIT, etc.)
invisible(run_once())

times <- numeric(reps)
scores <- integer(reps)
for (r in seq_len(reps)) {
  res <- run_once()
  times[r]  <- res$elapsed
  scores[r] <- res$score
  cat(sprintf("[%s] rep %d/%d: %.2f s  (score %d)\n",
              label, r, reps, res$elapsed, res$score))
}

cat(sprintf("\n[%s] reps=%d median=%.3f s  mean=%.3f s  min=%.3f s  max=%.3f s\n",
            label, reps, median(times), mean(times), min(times), max(times)))
cat(sprintf("[%s] scores: min=%d max=%d\n",
            label, min(scores), max(scores)))

# Emit machine-readable summary for downstream A/B aggregation
saveRDS(list(label = label, lib = lib_dir, reps = reps,
             times = times, scores = scores,
             median = median(times), mean = mean(times)),
        file = file.path("dev/profiling",
                         paste0("t300_na_bench_", label, ".rds")))
