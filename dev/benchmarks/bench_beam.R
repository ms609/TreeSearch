# Beam sectorial vs the single-tree baseline, from the canonical T0.
# Tests whether running RSS over a RETAINED diverse buffer (beam) escapes the
# frozen T0 where single-tree sectorial plateaus (~1267 on Zanol).  See
# dev/plans/2026-06-18-beam-sectorial.md.
#
# Budget is MATCHED: both arms run rssRounds x TS_RSS_PICKS = 30 x 20 = 600
# sector searches.  The ONLY differences in the beam arm: (a) each round starts
# from a beam-picked tree, not the cumulative single tree; (b) accept_equal is
# forced ON inside beam_sectorial (the diversity engine); (c) results written to
# a shared buffer.  coll30_20 single-tree reached only ~-4 (1267) at this budget.
#
# Arms (TS_BMARMS env, space-sep; default "single beam"):
#   single   coll30_20 single-tree (TS_BEAM unset, accept_equal FALSE) -- baseline
#   beam     same params, TS_BEAM=1 (best-equal buffer, accept_equal forced ON)
#   beamWide beam + TS_BEAM_SUBOPT band + TS_BEAM_PICKALL (Claim B, second test)
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem"),
            winslash = "/"))
  library(TreeTools)
})
dsN    <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zanol2014 Zhu2013 Wortley2006 Giles2015")), "\\s+")[[1]]
arms   <- strsplit(trimws(Sys.getenv("TS_BMARMS", "single beam")), "\\s+")[[1]]
ROUNDS <- as.integer(Sys.getenv("TS_RSSROUNDS", "30"))
PICKS  <- Sys.getenv("TS_RSS_PICKS", "20")
SEEDS  <- as.integer(strsplit(Sys.getenv("TS_SEEDS", "1 2 3"), "\\s+")[[1]])
SUBOPT <- Sys.getenv("TS_BEAM_SUBOPT", "5")   # band width for beamWide arm
target <- c(Zanol2014 = 1261, Wortley2006 = 480, Zhu2013 = 624, Giles2015 = 670)
t0dir  <- "dev/benchmarks/t0"

run_arm <- function(phy, t0, arm, seed) {
  set.seed(seed)
  Sys.setenv(TS_RSS_PICKS = PICKS)
  if (arm == "single") {
    Sys.unsetenv("TS_BEAM"); Sys.unsetenv("TS_BEAM_SUBOPT"); Sys.unsetenv("TS_BEAM_PICKALL")
  } else if (arm == "beam") {
    Sys.setenv(TS_BEAM = "1"); Sys.unsetenv("TS_BEAM_SUBOPT"); Sys.unsetenv("TS_BEAM_PICKALL")
  } else if (arm == "beamWide") {
    Sys.setenv(TS_BEAM = "1", TS_BEAM_SUBOPT = SUBOPT, TS_BEAM_PICKALL = "1")
  } else if (arm == "beamMulti") {
    # Full TNT analog: K diverse RAS+TBR seeds + wide buffer (retains them as
    # best drops) + pick-all (re-solves each seed's distinct descent).
    Sys.setenv(TS_BEAM = "1", TS_BEAM_SEEDS = Sys.getenv("TS_BEAM_SEEDS", "10"),
               TS_BEAM_SUBOPT = SUBOPT, TS_BEAM_PICKALL = "1")
  }
  r <- suppressWarnings(MaximizeParsimony(phy, tree = t0, maxReplicates = 1L, nThreads = 1L,
        maxSeconds = 0, verbosity = 0L, ratchetCycles = 0L, driftCycles = 0L,
        xssRounds = 0L, cssRounds = 0L, rssRounds = ROUNDS, wagnerStarts = 1L,
        fuseInterval = 9999L, sectorMinSize = 31L, sectorMaxSize = 99L,
        rasStarts = 3L, sectorCollapseTarget = 30L, sectorAcceptEqual = FALSE))
  Sys.unsetenv("TS_BEAM"); Sys.unsetenv("TS_BEAM_SUBOPT"); Sys.unsetenv("TS_BEAM_PICKALL")
  Sys.unsetenv("TS_RSS_PICKS")
  min(as.double(attr(r, "score")))
}

for (nm in dsN) {
  phy <- readRDS(file.path(t0dir, paste0(nm, ".phy.rds")))
  t0  <- ape::read.tree(file.path(t0dir, paste0(nm, ".tre")))
  t0len <- TreeLength(t0, phy); tgt <- target[[nm]]
  cat(sprintf("\n==== %s | T0=%.0f  target=%d (gap %+.0f) | budget %dx%s ====\n",
              nm, t0len, tgt, tgt - t0len, ROUNDS, PICKS))
  for (an in arms) {
    sc <- vapply(SEEDS, function(s) run_arm(phy, t0, an, s), double(1))
    best <- min(sc)
    cat(sprintf("  %-9s seeds[%s] -> %s | best %.0f (%+.0f vs T0, %+.0f vs target)%s\n",
                an, paste(SEEDS, collapse = ","), paste(format(sc), collapse = " "),
                best, best - t0len, best - tgt, if (best <= tgt) "  <== REACHED" else ""))
  }
}
