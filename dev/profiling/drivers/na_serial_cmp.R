# NA single-threaded bit-identical check: lever-1's change touches the shared
# Wagner/sector compute_insertion_edge_sets, which the NA path also uses.
# nThreads=1 sidesteps the pre-existing parallel race so we can compare the
# NA scoring path itself. Prints score + candidates_evaluated for one lib;
# run for both libs and diff.
lib <- Sys.getenv("TS_LIB", ".agent-l1")
suppressMessages({ library(TreeSearch, lib.loc = lib); library(TreeTools) })
data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]
for (sd in 1:4) {
  set.seed(sd)
  r <- suppressWarnings(MaximizeParsimony(ds, maxReplicates = 3L, targetHits = 2L,
                                          nThreads = 1L, strategy = "auto", verbosity = 0L))
  cat(sprintf("%s seed%d  score=%s  cand=%.0f\n", lib, sd,
              attr(r, "score"), as.numeric(attr(r, "candidates_evaluated"))))
}
