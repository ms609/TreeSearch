#!/usr/bin/env Rscript
# A/B: original (converging reduced-TBR) vs optimised (5-move limit + code fixes)
# Both libraries must already be installed.
#   orig_lib = .vtune-lib  (tbr_max_moves=0, no build_postorder fix)
#   opt_lib  = .agent-Eopt (tbr_max_moves=5, build_postorder deferred)
#
# Usage: Rscript dev/benchmarks/ab_pr_bench.R

orig_lib <- ".vtune-lib"
opt_lib  <- ".agent-Eopt"
datasets <- c("Zhu2013", "Dikow2009")
seeds    <- 1:5
budget   <- 20L

run_one <- function(lib, label, ds_name, seed, budget) {
  # Write temp script to avoid shell quoting issues
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', lib),
    'library(TreeSearch)',
    sprintf('ds <- inapplicable.phyData[["%s"]]', ds_name),
    sprintf('set.seed(%d)', seed),
    't0 <- proc.time()',
    sprintf(
      'res <- MaximizeParsimony(ds, maxSeconds=%dL, strategy="auto",',
      budget),
    '  pruneReinsertCycles=5L, pruneReinsertDrop=0.10,',
    '  driftCycles=0L, nniPerturbCycles=0L, verbosity=0L, nThreads=1L)',
    sprintf(
      'cat(sprintf("%s|%s|%d|%%g|%%d|%%.2f\\n",',
      label, ds_name, seed),
    '  attr(res,"score"), attr(res,"replicates"), (proc.time()-t0)[3])'
  ), tmp)
  out <- system2("Rscript", c("--no-save", tmp),
                 stdout = TRUE, stderr = FALSE)
  unlink(tmp)
  trimws(tail(out, 1))
}

results <- data.frame(
  label = character(), dataset = character(), seed = integer(),
  score = numeric(), reps = integer(), wall = numeric(),
  stringsAsFactors = FALSE
)

for (ds in datasets) {
  cat(sprintf("\n=== %s ===\n", ds))
  for (s in seeds) {
    for (cfg in list(list("orig", orig_lib), list("opt5", opt_lib))) {
      line <- run_one(cfg[[2]], cfg[[1]], ds, s, budget)
      cat(line, "\n")
      parts <- strsplit(line, "\\|")[[1]]
      if (length(parts) == 6) {
        results <- rbind(results, data.frame(
          label   = parts[1], dataset = parts[2], seed = as.integer(parts[3]),
          score   = as.numeric(parts[4]), reps = as.integer(parts[5]),
          wall    = as.numeric(parts[6]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

cat("\n=== Median summary (5 seeds) ===\n")
agg <- aggregate(cbind(score, reps, wall) ~ label + dataset, results, median)
print(agg[order(agg$dataset, agg$label), ], row.names = FALSE)

cat("\n=== Delta: opt5 vs orig ===\n")
for (ds in datasets) {
  orig <- agg[agg$label == "orig" & agg$dataset == ds, ]
  opt  <- agg[agg$label == "opt5" & agg$dataset == ds, ]
  cat(sprintf("%s: score %+.0f, reps %+.0f (%.0f%%), wall %+.2fs\n",
    ds,
    opt$score - orig$score,
    opt$reps  - orig$reps,
    (opt$reps - orig$reps) / orig$reps * 100,
    opt$wall  - orig$wall))
}
