# NA incremental-rescore mission-wall A/B cell (Hamilton job-array).
#
# Measures whether the incremental exact_verify rescore (commit 50046bfa,
# TS_NA_INCR) translates to MISSION wall on native NA. exact_verify is ~95% of
# native-NA mission wall and the incremental speeds it DIRECTLY (unlike the
# extract-fusion, which relocated work and washed) -- so the ~11-21% element win
# should carry. Paired arms, same node/seed:
#   off = legacy full_rescore per candidate (TS_NA_INCR unset)
#   on  = TS_NA_INCR=1 (incremental dirty rescore)
# Native NA (NOT recoded); XPIWE (MaximizeParsimony default); concavity=10.
# Identical score per cell also re-validates byte-identity at mission scale.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
                                              winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
idx  <- as.integer(if (length(args) >= 1L) args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
reps <- as.integer(Sys.getenv("TS_REPS", "5"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")), "\\s+")[[1]])
dsN  <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Dikow2009 Giles2015 Zanol2014 Zhu2013")), "\\s+")[[1]]
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_na_incr")

grid <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid))
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
row <- grid[idx + 1L, ]

data("inapplicable.phyData", package = "TreeSearch")
d <- inapplicable.phyData[[row$dataset]]   # NATIVE NA -- do NOT recode

run_arm <- function(arm) {
  if (arm == "on") Sys.setenv(TS_NA_INCR = "1") else Sys.unsetenv("TS_NA_INCR")
  set.seed(row$seed)
  t <- system.time(r <- suppressWarnings(MaximizeParsimony(
        d, concavity = 10, maxReplicates = reps, targetHits = 999L,
        maxSeconds = 0, nThreads = 1L, verbosity = 0L)))
  Sys.unsetenv("TS_NA_INCR")
  list(score = min(attr(r, "score")), wall = as.numeric(t[["elapsed"]]))
}
off <- run_arm("off")
on  <- run_arm("on")

out <- data.frame(dataset = row$dataset, seed = row$seed, reps = reps,
                  score_off = off$score, score_on = on$score,
                  wall_off = off$wall, wall_on = on$wall,
                  speedup = off$wall / on$wall,
                  identical = isTRUE(all.equal(off$score, on$score)),
                  stringsAsFactors = FALSE)
dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("naincr_%04d.csv", idx)), row.names = FALSE)
cat(sprintf("cell %d: %s seed %d | score %.4f/%.4f (%s) | wall %.1f/%.1f s => %.3fx\n",
            idx, row$dataset, row$seed, off$score, on$score,
            if (out$identical) "ok" else "DIFF", off$wall, on$wall, out$speedup))
