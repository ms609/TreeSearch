# Step 2: Is the 16-tree TNT-only set a separate TBR-plateau island?
# Seed TS from a TNT-only tree with ratchet/sector/drift/fuse DISABLED, so the
# only thing that runs is TBR (which cannot move an optimal tree) + the
# accept-equal MPT-enumeration plateau walk.  The returned pool is then the
# pure TBR-plateau closure of the island containing the seed.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(TreeTools))
base <- "C:/Users/pjjg18/GitHub/wide-sample/dev/zhu2013"
dat0 <- ReadAsPhyDat(file.path(base, "zhu2013_orig.nex"))
cmp <- readRDS("dev/profiling/zhu_tnt_ts_compare.rds")
og <- c("Galeaspida", "Osteostraci")
key1 <- function(t) ape::write.tree(SortTree(RootTree(t, og)))
k_tnt_only <- vapply(cmp$tnt_only, key1, character(1))
k_ts       <- vapply(cmp$ts, key1, character(1))

# minimal pipeline: TBR + enumeration only
minctl <- SearchControl(
  poolMaxSize = 20000L, tbrMaxHits = 200L, enumTimeFraction = 0.6,
  ratchetCycles = 0L, xssRounds = 0L, rssRounds = 0L, cssRounds = 0L,
  driftCycles = 0L, nniPerturbCycles = 0L, fuseInterval = 0L)

closure_from <- function(seed_tree, label) {
  set.seed(7L)
  res <- MaximizeParsimony(dat0, tree = seed_tree, concavity = Inf,
    inapplicable = "missing", maxReplicates = 1L, maxSeconds = 90,
    control = minctl, verbosity = 0L)
  kr <- vapply(res, key1, character(1))
  cat(sprintf("[%s] closure=%d | seed retained=%s | in TNT-only=%d/16 | in TS-256=%d/256 | new=%d\n",
      label, length(kr), key1(seed_tree) %in% kr,
      sum(kr %in% k_tnt_only), sum(kr %in% k_ts),
      sum(!(kr %in% k_tnt_only) & !(kr %in% k_ts))))
  kr
}

# Probe 3 different TNT-only seeds
for (i in c(1L, 6L, 12L)) {
  closure_from(cmp$tnt_only[[i]], paste0("tnt_only#", i))
}
