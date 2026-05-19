#!/usr/bin/env Rscript
# Sequential A/B: original vs optimised PR, run in subprocesses one at a time

orig_lib <- "C:/Users/pjjg18/GitHub/TreeSearch-a/.vtune-lib"
opt_lib  <- "C:/Users/pjjg18/GitHub/TreeSearch-a/.agent-Eopt"

run_bench <- function(lib, label, ds_name, seed, budget = 20L) {
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', gsub("\\\\", "/", lib)),
    'suppressPackageStartupMessages(library(TreeSearch))',
    sprintf('ds <- inapplicable.phyData[["%s"]]', ds_name),
    sprintf('set.seed(%dL)', seed),
    't0 <- proc.time()',
    sprintf('res <- MaximizeParsimony(ds, maxSeconds = %dL,', budget),
    '  strategy = "auto", pruneReinsertCycles = 5L,',
    '  pruneReinsertDrop = 0.10, driftCycles = 0L,',
    '  nniPerturbCycles = 0L, verbosity = 0L, nThreads = 1L)',
    'elapsed <- (proc.time() - t0)[[3]]',
    'cat(attr(res, "score"), attr(res, "replicates"), round(elapsed, 2))'
  ), tmp)
  out <- system2("Rscript", c("--no-save", tmp), stdout = TRUE, stderr = FALSE)
  unlink(tmp)
  vals <- as.numeric(strsplit(trimws(paste(out, collapse = " ")), "\\s+")[[1]])
  if (length(vals) >= 3) {
    cat(sprintf("  %s | %s | seed=%d | score=%g  reps=%d  wall=%.2fs\n",
                label, ds_name, seed, vals[1], vals[2], vals[3]))
    return(data.frame(label=label, dataset=ds_name, seed=seed,
                      score=vals[1], reps=as.integer(vals[2]), wall=vals[3]))
  } else {
    cat(sprintf("  %s | %s | seed=%d | FAILED\n", label, ds_name, seed))
    return(NULL)
  }
}

results <- list()
for (ds in c("Zhu2013", "Dikow2009")) {
  cat(sprintf("\n=== Dataset: %s ===\n", ds))
  for (s in 1:5) {
    results[[length(results)+1]] <- run_bench(orig_lib, "orig(cvg)", ds, s)
    Sys.sleep(0.5)  # let DLL unload
    results[[length(results)+1]] <- run_bench(opt_lib,  "opt(5mv)",  ds, s)
    Sys.sleep(0.5)
  }
}

df <- do.call(rbind, Filter(Negate(is.null), results))

cat("\n=== Median over 5 seeds ===\n")
agg <- aggregate(cbind(score, reps, wall) ~ label + dataset, df, median)
for (ds in unique(df$dataset)) {
  sub <- agg[agg$dataset == ds, ]
  orig <- sub[sub$label == "orig(cvg)", ]
  opt  <- sub[sub$label == "opt(5mv)",  ]
  cat(sprintf("\n%s:\n", ds))
  cat(sprintf("  orig(cvg): score=%.0f  reps=%d  wall=%.2fs\n",
              orig$score, orig$reps, orig$wall))
  cat(sprintf("  opt(5mv):  score=%.0f  reps=%d  wall=%.2fs\n",
              opt$score, opt$reps, opt$wall))
  cat(sprintf("  Delta:     score=%+.0f  reps=%+.0f (%+.0f%%)  wall=%+.2fs\n",
              opt$score - orig$score,
              opt$reps  - orig$reps,
              (opt$reps  - orig$reps) / orig$reps * 100,
              opt$wall  - orig$wall))
}
