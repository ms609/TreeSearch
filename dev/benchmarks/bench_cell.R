# Single-cell runner for SLURM job arrays (and local testing).
#
# Runs ONE (dataset x seed) cell of a panel and writes one partial CSV, so a
# panel can fan out across a SLURM --array (one task per cell). Replicate-bounded
# (deterministic candidates), nThreads=1. Merge partials with hamilton_merge.sh.
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID (0-based) into expand.grid(dataset, seed).
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_REPS, PARTIAL_DIR.
# Local test: Rscript dev/benchmarks/bench_cell.R 0

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
idx  <- as.integer(if (length(args) >= 1L) args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
reps <- as.integer(Sys.getenv("TS_REPS", "20"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")), "\\s+")[[1]])
dsN  <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009")), "\\s+")[[1]]
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials")

grid <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid))
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
row <- grid[idx + 1L, ]

data("inapplicable.phyData", package = "TreeSearch")
m <- PhyDatToMatrix(inapplicable.phyData[[row$dataset]], ambigNA = FALSE)
m[m == "-"] <- "?"
d <- MatrixToPhyDat(m)
set.seed(row$seed)
r <- suppressWarnings(MaximizeParsimony(d, maxReplicates = reps, targetHits = 999L,
                                        maxSeconds = 0, nThreads = 1L, verbosity = 0L))
out <- data.frame(dataset = row$dataset, seed = row$seed, score = attr(r, "score"),
                  candidates = attr(r, "candidates_evaluated"), stringsAsFactors = FALSE)
dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("cell_%04d.csv", idx)), row.names = FALSE)
cat(sprintf("cell %d: %s seed %d -> score %g, candidates %s\n",
            idx, row$dataset, row$seed, out$score, format(out$candidates, big.mark = ",")))
