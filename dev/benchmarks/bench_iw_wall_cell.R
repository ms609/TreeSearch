# IW mission-wall A/B cell runner (Hamilton job-array; also local-testable).
#
# Confirms the per-element IW opts (x4 reroot batch + extract_char_steps
# dirty-region) translate to WALL on the full MaximizeParsimony pipeline. Runs
# BOTH arms per cell (same node, same seed) to control machine variance:
#   off = TS_IW_NOX4=1 TS_IW_NODIRTY=1 (opts disabled, original kernel)
#   on  = opts enabled (default)
# PURE-IW only: the opts are gated !has_na, so the data MUST be recoded
# ("-"->"?") or the A/B measures ~0. concavity=10 (implied weights). Bounded by
# maxReplicates (fixed work, deterministic) so identical score validates
# byte-identity at scale and the wall ratio is the realized speedup.
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID (0-based) into expand.grid(dataset, seed).
# Env: TS_LIB, TS_DATASETS, TS_SEEDS, TS_REPS, PARTIAL_DIR.
# Local test: TS_REPS=1 TS_DATASETS=Vinther2008 TS_SEEDS=1 Rscript dev/benchmarks/bench_iw_wall_cell.R 0
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
                                              winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
idx  <- as.integer(if (length(args) >= 1L) args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
reps <- as.integer(Sys.getenv("TS_REPS", "10"))
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
dsN  <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Wortley2006 Zanol2014 Giles2015 Dikow2009")), "\\s+")[[1]]
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_iw")

grid <- expand.grid(dataset = dsN, seed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid))
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
row <- grid[idx + 1L, ]

data("inapplicable.phyData", package = "TreeSearch")
m <- PhyDatToMatrix(inapplicable.phyData[[row$dataset]], ambigNA = FALSE)
m[m == "-"] <- "?"                       # pure-IW (has_na = FALSE)
d <- MatrixToPhyDat(m)

run_arm <- function(arm) {
  if (arm == "off") { Sys.setenv(TS_IW_NOX4 = "1"); Sys.setenv(TS_IW_NODIRTY = "1") }
  else              { Sys.unsetenv("TS_IW_NOX4");   Sys.unsetenv("TS_IW_NODIRTY") }
  set.seed(row$seed)
  t <- system.time(r <- suppressWarnings(MaximizeParsimony(
        d, concavity = 10, maxReplicates = reps, targetHits = 999L,
        maxSeconds = 0, nThreads = 1L, verbosity = 0L)))
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
write.csv(out, file.path(partdir, sprintf("iwcell_%04d.csv", idx)), row.names = FALSE)
cat(sprintf("cell %d: %s seed %d | score %.4f/%.4f (%s) | wall %.1f/%.1f s => %.3fx\n",
            idx, row$dataset, row$seed, off$score, on$score,
            if (out$identical) "ok" else "DIFF", off$wall, on$wall, out$speedup))
