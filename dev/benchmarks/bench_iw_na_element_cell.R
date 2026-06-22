# NA-IW x4 ELEMENT-level A/B cell (dilution-free): times a direct ts_tbr_search
# climb from a RANDOM (poor) start to convergence -- the regime that exercises the
# full reroot loop where the 4-wide batch fires, with NO Wagner/ratchet/sectorial
# overhead to dilute it. This is the NA analog of bench_iw_realized.R and the
# decisive test of whether the NA-IW x4 is a real ELEMENT lever (like the no-NA
# x4's 1.1-1.2x direct-climb win) or genuinely at-limit even at the kernel.
#
# Paired arms on identical work (same start, set.seed before EACH call => byte-
# identical trajectory): off = TS_IW_NOX4=1 (scalar), on = unset (x4). Native NA,
# IW (concavity=10); the x4 kernel is shared by IW/XPIWE so IW-NA is representative.
# Each arm = median of TS_REPS repeated identical climbs (timing stability).
#
# Cell index: arg[1] or $SLURM_ARRAY_TASK_ID into expand.grid(dataset, startseed).
# Local test: TS_REPS=2 TS_DATASETS=Zhu2013 TS_SEEDS=1 Rscript dev/benchmarks/bench_iw_na_element_cell.R 0
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
partdir <- Sys.getenv("PARTIAL_DIR", "dev/benchmarks/partials_iw_na_elem")

grid <- expand.grid(dataset = dsN, startseed = seeds, stringsAsFactors = FALSE)
if (idx < 0L || idx >= nrow(grid))
  stop(sprintf("cell index %d out of range [0, %d)", idx, nrow(grid)))
row <- grid[idx + 1L, ]

data("inapplicable.phyData", package = "TreeSearch")
phy <- inapplicable.phyData[[row$dataset]]    # NATIVE NA
ct  <- attr(phy, "contrast"); lv <- attr(phy, "levels")
tip <- matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE)
wt  <- attr(phy, "weight")
mins <- as.integer(MinimumLength(phy, compress = TRUE))
set.seed(row$startseed)
redge <- RandomTree(names(phy), root = TRUE)$edge

time_arm <- function(x4on) {
  if (x4on) Sys.unsetenv("TS_IW_NOX4") else Sys.setenv(TS_IW_NOX4 = "1")
  walls <- numeric(reps); sc <- NA_real_
  for (r in seq_len(reps)) {
    set.seed(row$startseed)                   # identical trajectory each rep & arm
    t <- system.time(res <- TreeSearch:::ts_tbr_search(
           redge, ct, tip, wt, lv, maxHits = 1L, min_steps = mins, concavity = 10))
    walls[r] <- as.numeric(t[["elapsed"]]); sc <- res$score
  }
  Sys.unsetenv("TS_IW_NOX4")
  list(score = sc, wall = median(walls))
}
off <- time_arm(FALSE)
on  <- time_arm(TRUE)

out <- data.frame(dataset = row$dataset, startseed = row$startseed, reps = reps,
                  score_off = off$score, score_on = on$score,
                  wall_off = off$wall, wall_on = on$wall,
                  speedup = off$wall / on$wall,
                  identical = isTRUE(all.equal(off$score, on$score)),
                  stringsAsFactors = FALSE)
dir.create(partdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(partdir, sprintf("naelem_%04d.csv", idx)), row.names = FALSE)
cat(sprintf("cell %d: %s start %d | score %.4f/%.4f (%s) | wall %.2f/%.2f s => %.3fx\n",
            idx, row$dataset, row$startseed, off$score, on$score,
            if (out$identical) "ok" else "DIFF", off$wall, on$wall, out$speedup))
