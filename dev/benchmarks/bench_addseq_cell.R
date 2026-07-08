# Reachability A/B: does the TNT-style addition-sequence variety (CLOSEST/
# FURTHEST/INFORMATIVE wagnerBias modes) reach a given score in fewer
# replicates than RANDOM, at matched seed?
#
# Fixed-replicate-budget cells across dataset x wagnerBias mode x seed,
# mirroring bench_cell.R's replicate-bounded convention (never maxSeconds --
# under shared-node contention, replicates-completed depends on load, so a
# time-based stop is not reproducible; see dev/expertise/fast-iteration.md).
# Records the FULL per-replicate (replicate, elapsed, best_score) trajectory
# via progressCallback: replicate index is the deterministic/load-invariant
# axis; elapsed is recorded too but is a secondary, load-sensitive reading.
#
# The wall-clock COST of the extra order-generation work (the actual
# mission-relevant question) is measured separately in
# bench_addseq_ordergen.R, which times Wagner-tree construction alone in a
# single-tenant run. Net verdict = reachability gain (this script) minus
# order-gen cost (that script).
#
# adaptiveStart is forced FALSE and wagnerStarts forced to 1L for every arm
# (including bias=0 RANDOM) so the *only* difference between arms is the
# addition-sequence criterion -- otherwise the bandit or extra Wagner starts
# would dilute/select against the very mechanism under test (T-addseq).
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID (0-based) into
# expand.grid(dataset, bias, seed).
# Env: TS_LIB, TS_DATASETS, TS_BIASES, TS_SEEDS, TS_REPS, PARTIAL_DIR.
# Local test: Rscript dev/benchmarks/bench_addseq_cell.R 0

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
idx  <- as.integer(if (length(args) >= 1L) args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
reps <- as.integer(Sys.getenv("TS_REPS", "30"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5 6 7 8")), "\\s+")[[1]])
biases <- as.integer(strsplit(trimws(Sys.getenv("TS_BIASES", "0 3 4 5")), "\\s+")[[1]])
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
        "Longrich2010 Vinther2008 Sansom2010 Wortley2006 Eklund2004 Zanol2014 Zhu2013 Dikow2009")),
        "\\s+")[[1]]
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_addseq")

grid <- expand.grid(dataset = dsN, bias = biases, seed = seeds,
                    stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid)) {
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
}
row <- grid[idx + 1L, ]

data("inapplicable.phyData", package = "TreeSearch")
# Datasets resolve from inapplicable.phyData by name; anything not packaged
# (large off-corpus matrices: lobo, project175) is loaded from a directory of
# pre-converted .rds objects (TS_EXTRADIR). Pre-converting once -- rather than
# ReadAsPhyDat-ing raw .nex inside every one of the array's cells -- keeps the
# fragile native-inapplicable contrast build off the hot, 96x-replicated path
# (na-validation-alignment-gotcha) and lets it be verified a single time.
extradir <- Sys.getenv("TS_EXTRADIR", "")
d <- inapplicable.phyData[[row$dataset]]
if (is.null(d)) {
  f <- file.path(extradir, paste0(row$dataset, ".rds"))
  if (nzchar(extradir) && file.exists(f)) {
    d <- readRDS(f)
  } else {
    stop("Dataset not found: ", row$dataset,
         " (absent from inapplicable.phyData; no .rds at ", f, ")")
  }
}

traj <- list()
cb <- function(pi) {
  if (identical(pi$phase, "replicate")) {
    traj[[length(traj) + 1L]] <<- data.frame(
      replicate = pi$replicate, elapsed = pi$elapsed, best_score = pi$best_score
    )
  }
}

set.seed(row$seed)
r <- suppressWarnings(MaximizeParsimony(d,
  control = SearchControl(wagnerBias = row$bias, wagnerBiasTemp = 0.3,
                          adaptiveStart = FALSE, wagnerStarts = 1L),
  maxReplicates = reps, targetHits = 999L, maxSeconds = 0,
  nThreads = 1L, verbosity = 0L, progressCallback = cb))

out <- do.call(rbind, traj)
if (is.null(out)) out <- data.frame(replicate = integer(0), elapsed = numeric(0),
                                     best_score = numeric(0))
out$dataset <- row$dataset
out$bias <- row$bias
out$seed <- row$seed
out$final_score <- attr(r, "score")

dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("cell_%04d.csv", idx)), row.names = FALSE)
cat(sprintf("cell %d: %s bias=%d seed=%d -> final=%g (%d reps captured)\n",
            idx, row$dataset, row$bias, row$seed, out$final_score[1], nrow(out)))
