# Step 3: which diversity lever lets TS's NORMAL search (random starts, no
# external seed) reach the disconnected 16-tree island?  Measure recovery of
# the 16 TNT-only trees under different SearchControl configs.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(TreeTools))
base <- "C:/Users/pjjg18/GitHub/wide-sample/dev/zhu2013"
dat0 <- ReadAsPhyDat(file.path(base, "zhu2013_orig.nex"))
cmp <- readRDS("dev/profiling/zhu_tnt_ts_compare.rds")
og <- c("Galeaspida", "Osteostraci")
key1 <- function(t) ape::write.tree(SortTree(RootTree(t, og)))
k_tnt_only <- vapply(cmp$tnt_only, key1, character(1))   # the 16
k_tnt_all  <- vapply(cmp$tnt, key1, character(1))         # all 224

run <- function(label, ctl, maxRep, secs) {
  set.seed(2013L)
  res <- MaximizeParsimony(dat0, concavity = Inf, inapplicable = "missing",
    strategy = "intensive", maxReplicates = maxRep, targetHits = 9999L,
    maxSeconds = secs, control = ctl, verbosity = 0L)
  kr <- vapply(res, key1, character(1))
  cat(sprintf("[%-16s] pool=%4d | recovers %2d/16 island | total-vs-TNT shared=%d\n",
      label, length(kr), sum(k_tnt_only %in% kr), sum(k_tnt_all %in% kr)))
}

base_ctl <- function(...) SearchControl(poolMaxSize = 20000L, tbrMaxHits = 50L,
                                        enumTimeFraction = 0.4, ...)

# A: baseline intensive, more replicates (does raw start count help?)
run("more-reps",      base_ctl(),                              40L, 150)
# B: enable tree drifting (uphill tunnelling between islands; default OFF)
run("drift",          base_ctl(driftCycles = 15L),            40L, 150)
# C: TNT-style diverse retained set (accept-equal in sector + fuse)
run("diverse-set",    base_ctl(sectorAcceptEqual = TRUE, fuseAcceptEqual = TRUE,
                                cssRounds = 2L),                40L, 150)
