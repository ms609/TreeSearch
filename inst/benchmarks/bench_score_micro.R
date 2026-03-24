# Micro-benchmark: just Fitch scoring, no search
# Usage: Rscript inst/benchmarks/bench_score_micro.R <lib_path>
args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[1] else ".agent-pgo"

library(TreeSearch, lib.loc = lib_path)
library(TreeTools)

data("inapplicable.phyData")

prep_ds <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  storage.mode(contrast) <- "double"
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  storage.mode(tip_data) <- "integer"
  weight <- at$weight
  levels <- at$levels
  min_steps <- apply(contrast, 2, function(x) sum(x > 0)) - 1L
  min_steps <- pmax(min_steps, 0L)
  list(contrast = contrast, tip_data = tip_data, weight = weight,
       levels = levels, min_steps = min_steps)
}

for (nm in c("Agnarsson2004", "Dikow2009")) {
  ds <- inapplicable.phyData[[nm]]
  ds_args <- prep_ds(ds)
  
  set.seed(7294)
  tree <- RandomTree(names(ds), root = TRUE)
  edge <- tree$edge
  
  # Time many scoring calls
  n_iter <- 500L
  t0 <- system.time({
    for (i in seq_len(n_iter)) {
      TreeSearch:::ts_fitch_score(
        edge, ds_args$contrast, ds_args$tip_data,
        ds_args$weight, ds_args$levels, ds_args$min_steps
      )
    }
  })
  cat(nm, ": ", n_iter, " scores in ", t0["elapsed"], "s (",
      round(t0["elapsed"] / n_iter * 1000, 2), " ms/score)\n", sep = "")
}

# TBR phase breakdown
for (nm in c("Agnarsson2004", "Dikow2009")) {
  ds <- inapplicable.phyData[[nm]]
  ds_args <- prep_ds(ds)
  
  set.seed(7294)
  edge <- RandomTree(names(ds), root = TRUE)$edge
  
  r <- TreeSearch:::ts_bench_tbr_phases(
    edge, ds_args$contrast, ds_args$tip_data,
    ds_args$weight, ds_args$levels, ds_args$min_steps
  )
  cat(nm, " TBR: indirect=", r$time_indirect_us, "us, clip_incr=",
      r$time_clip_incr_us, "us, total_candidates=", r$n_candidates, "\n", sep = "")
}
