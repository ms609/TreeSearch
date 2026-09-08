# NA-IW x4 mission-wall A/B cell runner (Hamilton job-array; also local-testable).
#
# Measures whether the NA-IW 4-wide reroot batch (indirect_na_iw_cached_flat_x4,
# commit b34df50a) translates to WALL on the production MaximizeParsimony pipeline
# over the NATIVE inapplicable-bearing corpus. Runs BOTH arms per cell (same node,
# same seed) to control machine variance:
#   off = TS_IW_X4 unset (x4 disabled -> scalar indirect_na_iw_length_cached)
#   on  = TS_IW_X4=1     (x4 batch active; default-off opt-in)
#
# DELIBERATELY NATIVE NA (no "-"->"?" recode): the x4 NA branch is the whole point
# -- recoding would measure the no-NA path instead. The dirty-region opt is
# !has_na-gated so it is INACTIVE here; toggling TS_IW_X4 alone isolates the x4
# kernel cleanly (unlike the entangled no-NA x4+dirty A/B). concavity=10; the
# MaximizeParsimony default is extended IW => ScoringMode::XPIWE = the production
# path. Bounded by maxReplicates (fixed work, deterministic) so identical score
# validates byte-identity at scale and the wall ratio is the realized speedup.
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID (0-based) into expand.grid(dataset, seed).
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_REPS, PARTIAL_DIR.
# Local test: TS_REPS=1 TS_DATASETS=Zanol2014 TS_SEEDS=1 Rscript dev/benchmarks/bench_iw_na_wall_cell.R 0
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
                                              winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
idx  <- as.integer(if (length(args) >= 1L) args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
reps <- as.integer(Sys.getenv("TS_REPS", "8"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3 4 5")), "\\s+")[[1]])
dsN  <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Dikow2009 Giles2015 Zanol2014 Zhu2013")), "\\s+")[[1]]
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_iw_na")

grid <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid))
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
row <- grid[idx + 1L, ]

data("inapplicable.phyData", package = "TreeSearch")
d <- inapplicable.phyData[[row$dataset]]   # NATIVE NA -- do NOT recode

run_arm <- function(arm) {
  if (arm == "off") Sys.unsetenv("TS_IW_X4") else Sys.setenv(TS_IW_X4 = "1")
  set.seed(row$seed)
  t <- system.time(r <- suppressWarnings(MaximizeParsimony(
        d, concavity = 10, maxReplicates = reps, targetHits = 999L,
        maxSeconds = 0, nThreads = 1L, verbosity = 0L)))
  Sys.unsetenv("TS_IW_X4")
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
write.csv(out, file.path(partdir, sprintf("nacell_%04d.csv", idx)), row.names = FALSE)
cat(sprintf("cell %d: %s seed %d | score %.4f/%.4f (%s) | wall %.1f/%.1f s => %.3fx\n",
            idx, row$dataset, row$seed, off$score, on$score,
            if (out$identical) "ok" else "DIFF", off$wall, on$wall, out$speedup))
