#!/usr/bin/env Rscript
# Matched-wall A/B isolating the MARGINAL value of the gocomb FUSE sub-step over
# drift-only combined analysis. Pre-committed decision rule (set BEFORE running, to
# avoid open-ended tuning -- cf. the auto-tuner NO-GO):
#   * fuse arm moves the anytime curve below its matched drift-only arm (lower score
#     at equal wall, no dataset regressing) -> KEEP gocomb-fuse opt-in.
#   * fuse arm ties/does-not-improve its drift-only arm -> mechanism stays IN,
#     default-OFF, documented "implemented, no measured benefit over drift-only".
# Given pool>=2 is ~1/20 comb solves locally (RAS+drift starts converge), the tie
# outcome is the expectation, not a surprise.
#
# The comparison is clean regardless of engagement: within a combStarts level the
# fuse/drift arms are IDENTICAL except TS_SECT_FUSE, so any delta is fuse alone.
# Two combStarts levels (3, 6): more starts -> more start-diversity -> more chances
# for fuse to find >=2 distinct donors, at higher per-sector cost (matched wall
# absorbs the cost -> the anytime curve is the honest judge).
#
# Hamilton SLURM array, one cell per (dataset x seed x arm). No %throttle.
# Env per arm sets TS_SECT_FUSE / TS_SECT_DRIFT; rasStarts=3 in ALL arms.
# Usage (one array cell): Rscript sector_fuse_ab.R <cell_index>
suppressMessages({ library("TreeTools"); library("TreeSearch", lib.loc = Sys.getenv("TS_LIB")) })

datasets <- c("Zanol2014", "Zhu2013", "Wortley2006", "Giles2015")  # sector-resolve corpus
seeds    <- 1:5
# arm: goComb (0 => base/off), combStarts, fuse (TS_SECT_FUSE), forceOff (TS_SECT_DRIFT=0)
arms <- rbind(
  data.frame(name = "base",           goComb =  0L, combStarts = 3L, fuse = FALSE, forceOff = TRUE),
  data.frame(name = "comb_drift_s3",  goComb = 30L, combStarts = 3L, fuse = FALSE, forceOff = FALSE),
  data.frame(name = "comb_fuse_s3",   goComb = 30L, combStarts = 3L, fuse = TRUE,  forceOff = FALSE),
  data.frame(name = "comb_drift_s6",  goComb = 30L, combStarts = 6L, fuse = FALSE, forceOff = FALSE),
  data.frame(name = "comb_fuse_s6",   goComb = 30L, combStarts = 6L, fuse = TRUE,  forceOff = FALSE)
)
grid <- expand.grid(dataset = datasets, seed = seeds, arm = seq_len(nrow(arms)),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

cell <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
stopifnot(cell >= 1L, cell <= nrow(grid))
g   <- grid[cell, ]
arm <- arms[g$arm, ]

dat <- ReadAsPhyDat(file.path("data-raw", paste0(g$dataset, ".nex")))

ctrl <- SearchControl(
  xssRounds = 3L, rssRounds = 2L, cssRounds = 0L,
  sectorMinSize = 6L, sectorMaxSize = 80L,   # large sectors allowed in all arms
  rasStarts = 3L,
  sectorGoComb = arm$goComb, sectorCombStarts = arm$combStarts,
  sectorDriftCycles = 3L, sectorFuseRounds = 3L
)

# Env gates: force-off drift in base; disable only the fuse step in the drift arms.
if (isTRUE(arm$forceOff)) Sys.setenv(TS_SECT_DRIFT = "0") else Sys.unsetenv("TS_SECT_DRIFT")
if (isTRUE(arm$fuse)) Sys.unsetenv("TS_SECT_FUSE") else Sys.setenv(TS_SECT_FUSE = "0")
set.seed(g$seed)

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
                  arm = arm$name, goComb = arm$goComb, combStarts = arm$combStarts,
                  fuse = arm$fuse, final_score = attr(tr, "score")[1], wall = wall)
outDir <- Sys.getenv("AB_OUT", "dev/benchmarks/sector_fuse/out")
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
