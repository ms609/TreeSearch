# Stage 2 (HAMILTON): from the shared TNT T0, run OUR sectorial (ratchet/drift OFF)
# across two config arms x a BUDGET LADDER, recording score AND wall-clock per cell.
# The anytime curve (score vs wall) answers the MISSION GATE: does the selectem-
# config sectorial escape the frozen T0 at a budget CHEAP enough to beat the ratchet
# it would displace, or only at high budget (=> mission-dead)?
#
# Arms:
#   baseline  = current default sector config (band [6,~0.65n], rasStarts=1, no
#               collapse, auto picks). Expect ~0 escape (the documented null).
#   selectem  = TNT-faithful coll30-style: large clade band + sub-clade collapse +
#               rasStarts=3 + rssPicks=20 (sequential picks between global TBRs).
# Budget ladder: rssRounds in {1,2,4,8,16,32}. Wall is recorded so the two arms'
# different picks/round are still comparable on the score-vs-wall anytime curve.
# ratchet_ref = ONE reference run (default strategy, ratchet ON) from the same T0,
#   to stamp the wall/score yardstick the mission gate compares against (NOT a
#   recipe A/B — a single calibration point).
#
# One SLURM array task = one (dataset, seed) cell (reads manifest row TASK_ID).
# Env: TS_LIB, T0_DIR, OUT_DIR, TASK_ID (1-based; falls back to SLURM_ARRAY_TASK_ID).

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".missiongate-lib"),
                                              winslash = "/"))
  library(TreeTools)
})
t0dir  <- Sys.getenv("T0_DIR",  "dev/benchmarks/missiongate/t0")
outdir <- Sys.getenv("OUT_DIR", "dev/benchmarks/missiongate/out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tid <- as.integer(Sys.getenv("TASK_ID", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))

man <- read.csv(file.path(t0dir, "manifest.csv"), stringsAsFactors = FALSE)
if (tid < 1 || tid > nrow(man)) stop(sprintf("TASK_ID %d out of range 1..%d", tid, nrow(man)))
row <- man[tid, ]
nm <- row$dataset; sd <- row$seed; nTip <- row$nTip
start <- row$start; tnt_sect <- row$tnt_sect

data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fitch(inapplicable.phyData[[nm]])
t0  <- ape::read.tree(file.path(t0dir, row$t0_file))

# Confirm the loaded T0 scores == the manifest start (label-safe round-trip guard).
chk <- TreeLength(t0, phy)
if (abs(chk - start) > 1e-6)
  cat(sprintf("WARN %s s%d: T0 reloaded scores %.0f != manifest start %.0f\n", nm, sd, chk, start))

# Budget grid. The mission gate needs selectem sampled DOWN to ratchet-competitive
# wall (~0.3s), not just at high rounds -- so sweep picks x ras x rounds, spanning
# cheap (picks=5,ras=1,rounds=1) to expensive (picks=20,ras=3,rounds=16). Points
# are treated as wall-vs-escape Pareto SAMPLES (each an independent run from the
# fixed T0), robust to the per-point strict-descent noise. baseline + one ratchet
# yardstick complete the frontier.
SEL_PICKS  <- c(5L, 20L)
SEL_RAS    <- c(1L, 3L)
SEL_ROUNDS <- c(1L, 4L, 16L)
BASE_ROUNDS <- c(1L, 4L, 16L)

# One sectorial run from T0 at a given arm-config + budget. Ratchet/drift/xss/fuse
# OFF so the ONLY escape engine is RSS sectorial. Wall = elapsed of the call.
# ts_seed varies the sectorial RNG (the T0 is fixed): replicates separate a REAL
# cheap escape from a lucky random-sector draw (the per-point strict-descent noise).
run_sect <- function(minS, maxS, collapse, ras, picks, rounds, ts_seed) {
  set.seed(ts_seed)
  el <- system.time(
    r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
      nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
      ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
      wagnerStarts = 1L, fuseInterval = 9999L,
      sectorMinSize = minS, sectorMaxSize = maxS, sectorCollapseTarget = collapse,
      rasStarts = ras, rssPicks = picks, rssRounds = rounds))
  )["elapsed"]
  list(score = as.double(attr(r, "score")),
       cand = as.double(attr(r, "candidates_evaluated")), wall = as.double(el))
}

REPS <- as.integer(Sys.getenv("TS_REPS", "5"))   # sectorial RNG replicates per config

rows <- list()
add <- function(arm, rounds, res, picks, ras, ts_seed) {
  rows[[length(rows) + 1]] <<- data.frame(
    dataset = nm, seed = sd, nTip = nTip, arm = arm, rounds = rounds,
    rssPicks = picks, rasStarts = ras, ts_seed = ts_seed,
    start = start, tnt_sect = tnt_sect,
    score = res$score, escape = res$score - start, gap_to_tnt = res$score - tnt_sect,
    wall_s = round(res$wall, 3), Mcand = round(res$cand / 1e6, 3),
    stringsAsFactors = FALSE)
}

# Config bands (n-scaled). selectem = large backbone clades collapsed to ~0.4n units.
base_min <- 6L;                         base_max <- as.integer(round(nTip * 0.65))
sel_min  <- as.integer(round(nTip * 0.40)); sel_max <- as.integer(nTip)
sel_coll <- as.integer(round(nTip * 0.40))

summ <- function(tag, i0) {   # median/best escape + median wall over the REPS just added
  sl <- do.call(rbind, rows[(i0 + 1):length(rows)])
  cat(sprintf("%-12s s%d  %-28s esc med=%+d best=%+d  wall med=%.2fs\n",
              nm, sd, tag, as.integer(median(sl$escape)), as.integer(min(sl$escape)),
              median(sl$wall_s)))
}
for (b in BASE_ROUNDS) {
  i0 <- length(rows)
  for (ts in 1:REPS) add("baseline", b, run_sect(base_min, base_max, 0L, 1L, 0L, b, ts), 0L, 1L, ts)
  summ(sprintf("baseline rounds=%d", b), i0)
}
for (pk in SEL_PICKS) for (rs in SEL_RAS) for (b in SEL_ROUNDS) {
  i0 <- length(rows)
  for (ts in 1:REPS) add("selectem", b, run_sect(sel_min, sel_max, sel_coll, rs, pk, b, ts), pk, rs, ts)
  summ(sprintf("selectem pk=%d ras=%d rounds=%d", pk, rs, b), i0)
}

# ratchet_ref: default strategy (ratchet ON) from the same T0 — the wall/score
# yardstick the mission gate compares against. maxReplicates=1 so it's a single-
# start ratchet escape; replicate over ts_seed like the sectorial arms.
i0 <- length(rows)
for (ts in 1:REPS) {
  set.seed(ts)
  el <- system.time(
    rr <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
      nThreads = 1L, strategy = "default", maxSeconds = 0, verbosity = 0L))
  )["elapsed"]
  rows[[length(rows) + 1]] <- data.frame(
    dataset = nm, seed = sd, nTip = nTip, arm = "ratchet_ref", rounds = NA,
    rssPicks = NA, rasStarts = NA, ts_seed = ts, start = start, tnt_sect = tnt_sect,
    score = as.double(attr(rr, "score")),
    escape = as.double(attr(rr, "score")) - start,
    gap_to_tnt = as.double(attr(rr, "score")) - tnt_sect,
    wall_s = round(as.double(el), 3),
    Mcand = round(as.double(attr(rr, "candidates_evaluated")) / 1e6, 3),
    stringsAsFactors = FALSE)
}
summ("ratchet_ref", i0)

D <- do.call(rbind, rows)
of <- file.path(outdir, sprintf("cell_%02d_%s_s%d.csv", tid, nm, sd))
write.csv(D, of, row.names = FALSE)
cat(sprintf("Wrote %s\n", of))
