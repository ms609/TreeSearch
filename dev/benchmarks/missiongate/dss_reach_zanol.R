# DECISIVE reach-to-target check (advisor's cheapest disconfirmer).
# Question, binary: does SMALL-sector sectorial reach the Zanol s1 -13 target at ANY
# budget, or does it plateau above it? If it plateaus, "escape cheaper than the ratchet
# it would displace" is moot regardless of steps/Mcand -- the thread reframes to
# "cheap SHALLOW first-stage that can't reach the hard target" (~ rev5's verdict, and
# the rev2 whole-tree-barrier prior: a sector <= maxS tips can't emit a split spanning
# all 74). ratchet reached -13 (best) in the pilot; can small sectors?
#
# Best-of-many restart trajectories at solid budget (restart diversity is how barriers
# get crossed), plus an extreme-budget probe. RSS-only escape from the fixed T0.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".dss-lib"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
phy <- fitch(inapplicable.phyData[["Zanol2014"]])
t0  <- ape::read.tree("dev/benchmarks/missiongate/t0/Zanol2014_s1.tre")
start <- 1275; target <- 1262   # tnt_sect target => escape -13

run_sect <- function(maxS, rounds, seed) {
  set.seed(seed)
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L,
    sectorMinSize = 6L, sectorMaxSize = maxS, sectorCollapseTarget = 0L,
    rasStarts = 3L, rssPicks = 0L, rssRounds = rounds))
  as.double(attr(r, "score"))
}

grid <- rbind(
  expand.grid(maxS = c(12L, 20L, 30L), rounds = 64L, reps = 12L),
  expand.grid(maxS = 20L,              rounds = 256L, reps = 4L))

all_best <- Inf
for (i in seq_len(nrow(grid))) {
  mx <- grid$maxS[i]; rd <- grid$rounds[i]; rp <- grid$reps[i]
  sc <- vapply(seq_len(rp), function(s) run_sect(mx, rd, s), numeric(1))
  best <- min(sc); all_best <- min(all_best, best)
  cat(sprintf("maxS=%-3d rounds=%-3d reps=%-2d  best_score=%.0f  best_escape=%+d  reached_target(-13)=%s\n",
              mx, rd, rp, best, as.integer(best - start), best <= target))
}
cat(sprintf("\nOVERALL best small-sector score across all budgets: %.0f (escape %+d); target=%d (-13)\n",
            all_best, as.integer(all_best - start), target))
cat(sprintf("=> small-sector sectorial %s the -13 target at any tested budget.\n",
            if (all_best <= target) "REACHES" else "PLATEAUS ABOVE"))
