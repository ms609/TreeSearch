#!/usr/bin/env Rscript
# sprFirst gate -- JOB 1: the warmup -> first-TBR HANDOFF (fail-fast discriminator).
#
# Now that spr_search is an EXACT hill-climb (PR #273), is an SPR warmup a better
# handoff to TBR than the incumbent NNI warmup?  NNI c SPR c TBR, so TBR subsumes
# the SPR warmup; the ONLY place a better SPR basin can survive is where TBR is
# too expensive to fully wash it -- i.e. LARGE trees.  So the go/no-go is: on the
# large end, does SPR-first reach a better end-of-first-TBR score than NNI-first?
# (Advisor 2026-07-16: if the warmup delta is washed by end-of-first-TBR on the
# large end, the outer cycles wash it too -> NO-GO without the full pipeline.)
#
# Measurement: MaximizeParsimony with the ENTIRE perturbation pipeline stripped
# (ratchet/drift/sector/nni-perturb/prune/anneal/fuse all 0, outerCycles=1,
# maxReplicates=1, wagnerStarts=1).  The pipeline then reduces to exactly
#   Wagner start -> (warmup) -> TBR to convergence -> stop
# so the returned score IS the warmup->first-TBR handoff, via the real production
# code path (validates the actual sprFirst plumbing, not a reassembled proxy).
#
# wagnerStarts=1 is deliberate: it is SPR's BEST case.  At ws>1 the incumbent NNI
# path optimizes EVERY start then selects the best, whereas sprFirst selects by
# raw Wagner score then SPR-optimizes only the winner -- an asymmetry that only
# widens NNI's advantage.  If SPR-first can't win unhandicapped at ws=1, it won't
# at ws=3.
#
# Arms (all reachable on lib-driftexact; exact-SPR is its DEFAULT scorer):
#   nni    : nniFirst=TRUE , sprFirst=FALSE                 (incumbent / default)
#   spr_ex : nniFirst=FALSE, sprFirst=TRUE                  (exact SPR warmup)
#   spr_un : nniFirst=FALSE, sprFirst=TRUE , TS_SPR_UNION=1 (historical CONTROL --
#            reproduces the old "washed by TBR" test; if spr_ex ~ spr_un the
#            exactness fix reopened nothing -> fastest NO-GO with a clean why)
#   none   : nniFirst=FALSE, sprFirst=FALSE                 (no warmup baseline)
#
# Same seed per arm => IDENTICAL Wagner start (built before warmup diverges), so
# the four arms are paired on the starting topology.  EW-recoded (- -> ?) so the
# pure-EW spr_search path is exercised (NA/IW keep their own scorers, untouched).
# TRAINING split only (validation-set-sequestered); harness refuses non-training.
#
# Env knobs:
#   GATE_KEY      catalogue key (one per invocation / array task)
#   GATE_CAT      catalogue csv  (default /nobackup/$USER/floor/mbank_catalogue.csv)
#   GATE_MATRICES matrices dir   (resolved from ../neotrans/... if empty)
#   GATE_LIB      package library (default .agent-hj)
#   GATE_SEEDS    number of seeds / starting trees (default 5)
#   GATE_OUT      output csv      (default dev/profiling/sprfirst_<key>.csv)
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})

key     <- Sys.getenv("GATE_KEY")
catpath <- Sys.getenv("GATE_CAT", "/nobackup/pjjg18/floor/mbank_catalogue.csv")
mdir    <- Sys.getenv("GATE_MATRICES", "")
seeds   <- seq_len(as.integer(Sys.getenv("GATE_SEEDS", "5")))
outfile <- Sys.getenv("GATE_OUT", sprintf("dev/profiling/sprfirst_%s.csv", key))

cat0 <- read.csv(catpath, stringsAsFactors = FALSE)
row <- cat0[cat0$key == key, ]
if (nrow(row) != 1) stop("key not found / ambiguous: ", key)
if (row$split != "training")
  stop("REFUSING non-training split for tuning: ", key, " is ", row$split)
if (mdir == "") for (d in c("../neotrans/inst/matrices", "../neotrans/matrices"))
  if (dir.exists(d)) { mdir <- normalizePath(d); break }
nex <- file.path(mdir, row$filename)
if (!file.exists(nex)) stop("matrix not found: ", nex)

fitchPhy <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fitchPhy(ReadAsPhyDat(nex)); nt <- length(phy)
cat(sprintf("KEY %s : %d tips, %d chars (recoded EW), split=%s\n",
            key, nt, ncol(PhyDatToMatrix(phy)), row$split))

# Pipeline fully stripped: returned score == warmup -> first-TBR handoff.
# Pinned strategy "sprint" (explicit, not "auto") so there is NO size-based
# preset swap across the size-stratified keys; the control below overrides every
# field that matters (control merges over the preset, explicit fields win --
# control-clobbers-strategy-preset is FIXED on cpp-search).
handoffCtrl <- function(nniFirst, sprFirst) SearchControl(
  wagnerStarts   = 1L,
  nniFirst       = nniFirst,
  sprFirst       = sprFirst,
  tbrMaxHits     = 1L,          # production TBR behaviour
  ratchetCycles  = 0L, driftCycles = 0L, nniPerturbCycles = 0L,
  pruneReinsertCycles = 0L, annealCycles = 0L,
  xssRounds = 0L, rssRounds = 0L, cssRounds = 0L, sectorGoDrift = 0L,
  outerCycles = 1L, maxOuterResets = 0L, fuseInterval = 0L,
  adaptiveLevel = FALSE, adaptiveStart = FALSE, consensusStableReps = 0L
)

arms <- list(
  nni    = list(nniFirst = TRUE,  sprFirst = FALSE, union = FALSE),
  spr_ex = list(nniFirst = FALSE, sprFirst = TRUE,  union = FALSE),
  spr_un = list(nniFirst = FALSE, sprFirst = TRUE,  union = TRUE),
  none   = list(nniFirst = FALSE, sprFirst = FALSE, union = FALSE)
)

rows <- list()
for (s in seeds) for (anm in names(arms)) {
  a <- arms[[anm]]
  if (a$union) Sys.setenv(TS_SPR_UNION = "1") else Sys.unsetenv("TS_SPR_UNION")
  set.seed(s)                                    # identical Wagner start per arm
  t0 <- Sys.time()
  res <- MaximizeParsimony(phy, strategy = "sprint",
                           control = handoffCtrl(a$nniFirst, a$sprFirst),
                           maxReplicates = 1L, maxSeconds = 0, nThreads = 1L,
                           verbosity = 0L)
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  Sys.unsetenv("TS_SPR_UNION")
  sc <- as.double(attr(res, "score"))
  rows[[length(rows) + 1L]] <- data.frame(
    key = key, nTip = nt, nChar = ncol(PhyDatToMatrix(phy)),
    seed = s, arm = anm, score = sc, wall_s = round(wall, 3))
  cat(sprintf("%-14s seed=%d %-6s handoff=%.0f wall=%.2fs\n", key, s, anm, sc, wall))
  write.csv(do.call(rbind, rows), outfile, row.names = FALSE)  # incremental
}
cat(sprintf("Wrote %s (%d rows)\n", outfile, length(rows)))

# --- paired summary: does SPR-first beat NNI-first at end-of-first-TBR? ---
df <- do.call(rbind, rows)
ms <- function(a) mean(df$score[df$arm == a]); mw <- function(a) mean(df$wall_s[df$arm == a])
cat("\n== HANDOFF summary (mean over seeds) ==\n")
for (a in names(arms)) cat(sprintf("  %-6s score=%.1f  wall=%.2fs\n", a, ms(a), mw(a)))
cat(sprintf("\n  spr_ex - nni : dScore=%+.1f  dWall=%+.2fs\n", ms("spr_ex")-ms("nni"), mw("spr_ex")-mw("nni")))
cat(sprintf("  spr_ex - spr_un (exactness effect): dScore=%+.1f\n", ms("spr_ex")-ms("spr_un")))
# paired per-seed win count spr_ex vs nni
w <- sapply(seeds, function(s) sign(df$score[df$arm=="nni"&df$seed==s] - df$score[df$arm=="spr_ex"&df$seed==s]))
cat(sprintf("  spr_ex beats nni on %d/%d seeds (ties %d)\n", sum(w>0), length(seeds), sum(w==0)))
