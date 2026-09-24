# DSS SIZE-AXIS PROBE (zero-code precondition for the diameter-limited-sector build).
#
# The dss hypothesis: SMALLER, locally-bounded sectors cut per-sector rebuild cost,
# so sectorial escape becomes cheaper per candidate -- possibly cheap enough to FLIP
# the rev5 ranking (ratchet -1.34 vs sectorial -0.25 steps/Mcand on Zanol).
#
# The PRIMARY driver of "smaller sector" is TIP-COUNT, controllable via the EXISTING
# sectorMaxSize knob -- and rev5's sweep NEVER varied it small (baseline held maxS at
# 0.65n, selectem at 1.0n; only picks/ras/rounds were swept). Diameter-limiting is a
# strictly-MORE-restrictive geometry ON TOP of size, so if shrinking maxS does not
# climb escape/candidate toward ratchet's -1.34, the diameter build cannot rescue it
# and the thread is mission-dead. If small maxS DOES climb, the build is justified and
# this tells us the size regime to centre the diameter sweep on. (rev4/rev5 discipline:
# config-first, existing knobs, no new src; write the routine only if the proxy shows
# the effect.)
#
# Arms, ALL from the identical fixed TNT T0, RSS-only escape (ratchet/drift/xss/fuse
# OFF -- exactly sweep_missiongate.R's run_sect), reps over the sectorial RNG:
#   smallmax        maxS in {12,20,30}, minS=6, collapse OFF, ras=3, picks=auto.
#                   The NEW small-LOCAL-sector regime rev5 never measured. All three
#                   maxS values sit below rev5's baseline band (0.65n) on the larger
#                   datasets. ras=3 because rev5 found it is the ONLY ingredient that
#                   moves escape at all; auto picks give each sector size its natural
#                   per-round pick budget (2n/avg_size).
#   selectem_large  minS=0.40n, maxS=n, collapse=0.40n, ras=3, picks=auto. The rev5
#                   best-efficiency LARGE-sector recipe, re-run on THIS engine as the
#                   within-run "large sector" anchor (~ -0.25 steps/Mcand on Zanol).
#   ratchet_ref     default strategy, maxReplicates=1. The wall/candidate yardstick
#                   the mission gate compares against (~ -1.34 steps/Mcand on Zanol).
# Budget ladder rounds in {1,4,16} traces the escape-vs-candidate curve so efficiency
# is read across budgets, not at one point.
#
# GUARDRAIL (rev2 + advisor): candidates_evaluated goes NEGATIVE when differenced
# across configs -- NEVER subtract Mcand across arms. Each run reports its OWN escape
# and Mcand from the fixed T0; efficiency = escape/Mcand computed WITHIN a run (or, for
# comparability with rev5's table, esc_med/Mcand_med in the aggregate). We never diff.
#
# One SLURM array task = one (dataset, seed) manifest row (TASK_ID).
# Env: TS_LIB, T0_DIR, OUT_DIR, TASK_ID (falls back to SLURM_ARRAY_TASK_ID), TS_REPS.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".missiongate-lib"),
                                              winslash = "/"))
  library(TreeTools)
})
t0dir  <- Sys.getenv("T0_DIR",  "dev/benchmarks/missiongate/t0")
outdir <- Sys.getenv("OUT_DIR", "dev/benchmarks/missiongate/dss_out")
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

chk <- TreeLength(t0, phy)
if (abs(chk - start) > 1e-6)
  cat(sprintf("WARN %s s%d: T0 reloaded scores %.0f != manifest start %.0f\n", nm, sd, chk, start))

# One sectorial run from T0. RSS is the ONLY escape engine (all else OFF). Identical
# to sweep_missiongate.R::run_sect so numbers are directly comparable to rev5.
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

REPS <- as.integer(Sys.getenv("TS_REPS", "5"))

SMALL_MAXS <- c(12L, 20L, 30L)       # small LOCAL sectors (below rev5's 0.65n / 1.0n bands)
SMALL_ROUNDS <- c(1L, 4L, 16L, 64L)  # smallmax budget ladder (64 = per-cell reach ceiling)
ROUNDS       <- c(1L, 4L, 16L)       # anchor/ladder for selectem_large

sel_min  <- as.integer(round(nTip * 0.40))
sel_max  <- as.integer(nTip)
sel_coll <- as.integer(round(nTip * 0.40))

rows <- list()
add <- function(arm, maxS, minS, collapse, ras, picks, rounds, res, ts_seed) {
  rows[[length(rows) + 1]] <<- data.frame(
    dataset = nm, seed = sd, nTip = nTip, arm = arm,
    minS = minS, maxS = maxS, collapse = collapse, rasStarts = ras,
    rssPicks = picks, rounds = rounds, ts_seed = ts_seed,
    start = start, tnt_sect = tnt_sect,
    score = res$score, escape = res$score - start, gap_to_tnt = res$score - tnt_sect,
    wall_s = round(res$wall, 3), Mcand = round(res$cand / 1e6, 4),
    stringsAsFactors = FALSE)
}
summ <- function(tag, i0) {
  sl <- do.call(rbind, rows[(i0 + 1):length(rows)])
  eff <- if (median(sl$Mcand) > 0) median(sl$escape) / median(sl$Mcand) else NA
  cat(sprintf("%-12s s%d  %-30s esc med=%+d best=%+d  Mcand med=%.2f  eff=%.2f st/Mc  wall med=%.2fs\n",
              nm, sd, tag, as.integer(median(sl$escape)), as.integer(min(sl$escape)),
              median(sl$Mcand), eff, median(sl$wall_s)))
}

# smallmax arm (the new small-LOCAL-sector regime)
for (mx in SMALL_MAXS) for (b in SMALL_ROUNDS) {
  i0 <- length(rows)
  for (ts in 1:REPS) add("smallmax", mx, 6L, 0L, 3L, 0L, b,
                         run_sect(6L, mx, 0L, 3L, 0L, b, ts), ts)
  summ(sprintf("smallmax maxS=%d rounds=%d", mx, b), i0)
}

# selectem_large anchor (rev5 large-collapsed recipe, this engine)
for (b in ROUNDS) {
  i0 <- length(rows)
  for (ts in 1:REPS) add("selectem_large", sel_max, sel_min, sel_coll, 3L, 0L, b,
                         run_sect(sel_min, sel_max, sel_coll, 3L, 0L, b, ts), ts)
  summ(sprintf("selectem_large rounds=%d", b), i0)
}

# ratchet_ref yardstick (default strategy, single-start)
i0 <- length(rows)
for (ts in 1:REPS) {
  set.seed(ts)
  el <- system.time(
    rr <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
      nThreads = 1L, strategy = "default", maxSeconds = 0, verbosity = 0L))
  )["elapsed"]
  rows[[length(rows) + 1]] <- data.frame(
    dataset = nm, seed = sd, nTip = nTip, arm = "ratchet_ref",
    minS = NA, maxS = NA, collapse = NA, rasStarts = NA, rssPicks = NA, rounds = NA,
    ts_seed = ts, start = start, tnt_sect = tnt_sect,
    score = as.double(attr(rr, "score")),
    escape = as.double(attr(rr, "score")) - start,
    gap_to_tnt = as.double(attr(rr, "score")) - tnt_sect,
    wall_s = round(as.double(el), 3),
    Mcand = round(as.double(attr(rr, "candidates_evaluated")) / 1e6, 4),
    stringsAsFactors = FALSE)
}
summ("ratchet_ref", i0)

D <- do.call(rbind, rows)
of <- file.path(outdir, sprintf("cell_%02d_%s_s%d.csv", tid, nm, sd))
write.csv(D, of, row.names = FALSE)
cat(sprintf("Wrote %s\n", of))
