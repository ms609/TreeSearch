#!/usr/bin/env Rscript
# Matched-wall A/B for TNT-style in-sector drifting (godrift, HTU-PINNED variant)
# on the hard sector-resolve / floor corpus. Isolates the DRIFT contribution: both
# arms allow large sectors (sectorMaxSize raised) and use rasStarts>=2 so the
# sector-resolve retention channel is live (see memory sect-slack-ablation /
# tnt-sectorial-recipe); they differ ONLY in whether large sectors are solved by
# drift. The drift is HTU-pinned (content sector_mask) so ~88% of solves are kept
# and TBR-polished (see memory sector-drift-htu-float); this A/B answers the open
# question that the local revert-rate metric CANNOT: does pinned drift move the
# anytime curve (time-to-best), not just fire.
#
# Mission = wall-clock time-to-optimum, so this reports the ANYTIME curve (score
# vs elapsed) per cell, not just the final score. Drift is expensive and runs per
# sector pick, so the likely failure mode is "better per-sector quality, WORSE
# time-to-optimum" — hence the threshold × cycles SWEEP, not a single on/off point.
#
# Hamilton SLURM array: one cell per (dataset × seed × arm). No %throttle.
# Env: TERM=dumb OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1, dedicated isolated lib.
#
# Usage (one array cell): Rscript sector_drift_ab.R <cell_index>
suppressMessages({ library("TreeTools"); library("TreeSearch", lib.loc = Sys.getenv("TS_LIB")) })

# ---- Grid ----------------------------------------------------------------
datasets <- c("Zanol2014", "Zhu2013", "Wortley2006", "Giles2015")  # sector-resolve corpus
seeds    <- 1:5
# arm 0 = baseline (RAS+TBR only, drift force-OFF); arms 1.. sweep threshold×cycles.
arms <- rbind(
  data.frame(name = "base",         goDrift =  0L, cycles = 0L, forceOff = TRUE),
  data.frame(name = "drift40_c3",   goDrift = 40L, cycles = 3L, forceOff = FALSE),
  data.frame(name = "drift40_c6",   goDrift = 40L, cycles = 6L, forceOff = FALSE),
  data.frame(name = "drift25_c3",   goDrift = 25L, cycles = 3L, forceOff = FALSE),
  data.frame(name = "drift55_c3",   goDrift = 55L, cycles = 3L, forceOff = FALSE)
)
grid <- expand.grid(dataset = datasets, seed = seeds, arm = seq_len(nrow(arms)),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

cell <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
stopifnot(cell >= 1L, cell <= nrow(grid))
g   <- grid[cell, ]
arm <- arms[g$arm, ]

datPath <- file.path("data-raw", paste0(g$dataset, ".nex"))
dat <- ReadAsPhyDat(datPath)

# Matched wall: fixed maxSeconds budget, both arms identical except drift.
# rasStarts=3 (TNT) in BOTH arms so the mechanism is exercised, not inert.
ctrl <- SearchControl(
  xssRounds = 3L, rssRounds = 2L, cssRounds = 0L,
  sectorMinSize = 6L, sectorMaxSize = 80L,     # large sectors allowed in BOTH arms
  rasStarts = 3L,
  sectorGoDrift = arm$goDrift,
  sectorDriftCycles = if (arm$cycles > 0L) arm$cycles else 5L
)

if (isTRUE(arm$forceOff)) Sys.setenv(TS_SECT_DRIFT = "0") else Sys.unsetenv("TS_SECT_DRIFT")
set.seed(g$seed)

# Anytime trajectory via progressCallback (best_score, elapsed) — the mission
# metric. The callback fires per replicate/phase (coarse but matched across arms).
traj <- list()
cb <- function(info) {
  traj[[length(traj) + 1L]] <<- data.frame(t = as.double(info$elapsed),
                                            score = as.double(info$best_score))
  invisible(TRUE)
}
t0 <- Sys.time()
tr <- MaximizeParsimony(dat, maxReplicates = 200L, maxSeconds = 90,
                        nThreads = 1L, verbosity = 0L, strategy = "none",
                        control = ctrl, progressCallback = cb)
wall <- as.double(difftime(Sys.time(), t0, units = "secs"))

out <- data.frame(cell = cell, dataset = g$dataset, seed = g$seed,
                  arm = arm$name, goDrift = arm$goDrift, cycles = arm$cycles,
                  final_score = attr(tr, "score")[1], wall = wall)
outDir <- Sys.getenv("AB_OUT", "dev/benchmarks/sector_drift/out")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(outDir, sprintf("cell_%03d_%s_s%d_%s.csv",
          cell, g$dataset, g$seed, arm$name)), row.names = FALSE)
if (length(traj)) {
  tj <- do.call(rbind, traj); tj$cell <- cell; tj$arm <- arm$name
  write.csv(tj, file.path(outDir, sprintf("traj_%03d_%s_s%d_%s.csv",
            cell, g$dataset, g$seed, arm$name)), row.names = FALSE)
}
cat(sprintf("cell %d %s s%d [%s] final=%.1f wall=%.1fs\n",
            cell, g$dataset, g$seed, arm$name, out$final_score, wall))
