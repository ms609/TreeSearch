# Real-matrix probe: is Zhu 256 the plateau closure, and do those trees
# collapse to far fewer once zero-length (unsupported) branches are removed?
suppressMessages(devtools::load_all(".", quiet = TRUE))
nex <- "C:/Users/pjjg18/GitHub/wide-sample/dev/zhu2013/zhu2013_orig.nex"
dat <- TreeTools::ReadAsPhyDat(nex)

set.seed(2013L)
trees <- MaximizeParsimony(
  dat, concavity = Inf, inapplicable = "missing",
  strategy = "intensive", maxReplicates = 15L, maxSeconds = 45,
  control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 50L),
  verbosity = 2L)

cat("\n==== RESULT ====\n")
cat("returned resolved MPTs:", length(trees), "\n")
cat("score:", attr(trees, "score"), "\n")
resolved <- vapply(trees, function(t)
  ape::write.tree(TreeTools::SortTree(ape::unroot(t))), character(1))
cat("distinct resolved (sorted Newick):", length(unique(resolved)), "\n")
saveRDS(trees, "dev/profiling/zhu_probe_trees.rds")
cat("saved trees to dev/profiling/zhu_probe_trees.rds\n")
