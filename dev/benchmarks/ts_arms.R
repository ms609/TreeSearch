# TreeSearch shared-start arms vs the TNT ratchet-off target (define_target.R).
# Starts from the IDENTICAL canonical T0 (dev/benchmarks/t0/<nm>.tre), ratchet/drift/
# fuse OFF, rss-only.  Verifies TreeLength(T0)==expected before searching.
#
# Arms (TS_ARMS env, space-sep; default "base coll30"):
#   base         defaults [6,50] ras1 coll0                  -- current behaviour
#   coll30       [31,99] ras3 break-big collapse->30 units   -- the established null
#   freezeDet    [31,99] ras3 freeze-big (cap15 thr8) DET     -- H2: large movable units
#   freezeRand   [31,99] ras3 freeze-big (cap15 thr8) RANDOM  -- H1: per-pass diversity
#   freezeRand20 freezeRand + rss_picks=20                    -- + pick count
# Freeze arms route through build_reduced_dataset_freeze (TS_FREEZE_COLLAPSE);
# byte-identical to current code when unset.  Point TS_LIB at the freeze build.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband2"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Wortley2006 Zhu2013 Giles2015")), "\\s+")[[1]]
arms   <- strsplit(trimws(Sys.getenv("TS_ARMS", "base coll30")), "\\s+")[[1]]
ROUNDS <- as.integer(Sys.getenv("TS_RSSROUNDS", "15"))
SEEDS  <- as.integer(strsplit(Sys.getenv("TS_SEEDS", "1 2 3"), "\\s+")[[1]])
target <- c(Zanol2014 = 1261, Wortley2006 = 480, Zhu2013 = 624, Giles2015 = 670)
t0dir  <- "dev/benchmarks/t0"

# arm = list(min, max, ras, coll, eq, freeze, rand, cap, thresh, picks)
cfg <- list(
  base         = list(6L,  50L, 1L,  0L, FALSE, 0L, 0L,  0L, 0L,  0L),
  coll30       = list(31L, 99L, 3L, 30L, FALSE, 0L, 0L,  0L, 0L,  0L),
  base20       = list(6L,  50L, 1L,  0L, FALSE, 0L, 0L,  0L, 0L, 20L),  # budget-matched null
  coll30_20    = list(31L, 99L, 3L, 30L, FALSE, 0L, 0L,  0L, 0L, 20L),  # budget-matched null
  largeOnly    = list(31L, 99L, 1L,  0L, FALSE, 0L, 0L,  0L, 0L, 20L),  # large band, ras1, NO collapse
  largeRas3    = list(31L, 99L, 3L,  0L, FALSE, 0L, 0L,  0L, 0L, 20L),  # large band, ras3, NO collapse
  freezeDet    = list(31L, 99L, 3L,  0L, FALSE, 1L, 0L, 15L, 8L,  0L),
  freezeRand   = list(31L, 99L, 3L,  0L, FALSE, 1L, 1L, 15L, 8L,  0L),
  freezeRand20 = list(31L, 99L, 3L,  0L, FALSE, 1L, 1L, 15L, 8L, 20L),
  freezeHT     = list(31L, 99L, 3L,  0L, FALSE, 1L, 1L, 33L, 28L, 20L),  # high thr, overshoot
  freezeHT2    = list(31L, 99L, 3L,  0L, FALSE, 1L, 1L, 40L, 30L, 20L),
  freezeHTdet  = list(31L, 99L, 3L,  0L, FALSE, 1L, 0L, 33L, 28L, 20L),  # H1/H2 ablation: DET
  # n-scaled (negative field = percent of NTip): min .42n max .99n cap .45n thr .38n
  freezeScaled = list(-42L, -99L, 3L, 0L, FALSE, 1L, 1L, -45L, -38L, 20L),
  freezeScalDet= list(-42L, -99L, 3L, 0L, FALSE, 1L, 0L, -45L, -38L, 20L)
)

run_arm <- function(phy, t0, a, seed) {
  set.seed(seed)
  n <- NTip(phy)
  res <- function(v) if (v < 0L) as.integer(round(n * (-v) / 100)) else v   # neg = pct of n
  a[[1]] <- res(a[[1]]); a[[2]] <- res(a[[2]]); a[[8]] <- res(a[[8]]); a[[9]] <- res(a[[9]])
  if (a[[6]] > 0) {
    Sys.setenv(TS_FREEZE_COLLAPSE = "1",
               TS_FREEZE_CAP = as.character(a[[8]]),
               TS_FREEZE_THRESH = as.character(a[[9]]))
    if (a[[7]] > 0) Sys.setenv(TS_FREEZE_RANDOM = "1") else Sys.unsetenv("TS_FREEZE_RANDOM")
  } else {
    Sys.unsetenv("TS_FREEZE_COLLAPSE"); Sys.unsetenv("TS_FREEZE_RANDOM")
    Sys.unsetenv("TS_FREEZE_CAP"); Sys.unsetenv("TS_FREEZE_THRESH")
  }
  if (a[[10]] > 0) Sys.setenv(TS_RSS_PICKS = as.character(a[[10]])) else Sys.unsetenv("TS_RSS_PICKS")
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
        maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
        xssRounds = 0L, cssRounds = 0L, rssRounds = ROUNDS, wagnerStarts = 1L,
        fuseInterval = 9999L, sectorMinSize = a[[1]], sectorMaxSize = a[[2]],
        rasStarts = a[[3]], sectorCollapseTarget = a[[4]], sectorAcceptEqual = a[[5]]))
  Sys.unsetenv("TS_FREEZE_COLLAPSE"); Sys.unsetenv("TS_FREEZE_RANDOM")
  Sys.unsetenv("TS_FREEZE_CAP"); Sys.unsetenv("TS_FREEZE_THRESH"); Sys.unsetenv("TS_RSS_PICKS")
  min(as.double(attr(r, "score")))
}

for (nm in dsN) {
  phy <- readRDS(file.path(t0dir, paste0(nm, ".phy.rds")))
  t0  <- ape::read.tree(file.path(t0dir, paste0(nm, ".tre")))
  t0len <- TreeLength(t0, phy); tgt <- target[[nm]]
  cat(sprintf("\n==== %s | T0=%.0f  target=%d (gap %+.0f) ====\n", nm, t0len, tgt, tgt - t0len))
  for (an in arms) {
    a <- cfg[[an]]
    sc <- vapply(SEEDS, function(s) run_arm(phy, t0, a, s), double(1))
    best <- min(sc)
    cat(sprintf("  %-12s seeds[%s] -> %s | best %.0f (%+.0f vs T0, %+.0f vs target)%s\n",
                an, paste(SEEDS, collapse = ","), paste(format(sc), collapse = " "),
                best, best - t0len, best - tgt, if (best <= tgt) "  <== REACHED" else ""))
  }
}
