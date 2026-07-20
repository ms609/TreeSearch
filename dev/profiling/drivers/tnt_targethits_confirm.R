# Confirm: the default stopping rule (targetHits = nTip/5) halts the search
# once the best SCORE is confirmed, long before the full MPT SET is enumerated.
# Completeness needs many independent replicates -> raise targetHits.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(TreeTools))
base <- "C:/Users/pjjg18/GitHub/wide-sample/dev/zhu2013"
dat0 <- ReadAsPhyDat(file.path(base, "zhu2013_orig.nex"))
cmp <- readRDS("dev/profiling/zhu_tnt_ts_compare.rds")
og <- c("Galeaspida", "Osteostraci")
key1 <- function(t) ape::write.tree(SortTree(RootTree(t, og)))
k16 <- vapply(cmp$tnt_only, key1, character(1))

run <- function(label, th) {
  set.seed(2013L)
  res <- MaximizeParsimony(dat0, concavity = Inf, inapplicable = "missing",
    strategy = "intensive", maxReplicates = 60L, targetHits = th,
    maxSeconds = 200,
    control = SearchControl(poolMaxSize = 20000L, tbrMaxHits = 50L,
                            enumTimeFraction = 0.4),
    verbosity = 0L)
  kr <- vapply(res, key1, character(1))
  cat(sprintf("[%-14s targetHits=%-5s] reps_run=%s | pool=%d | recovers %d/16 | perturb_stop=%s timed_out=%s\n",
      label, th, attr(res, "replicates"), length(kr), sum(k16 %in% kr),
      attr(res, "perturb_stop"), attr(res, "timed_out")))
}

run("default-stop", max(10L, as.integer(NTip(dat0) / 5)))  # = 15, the default
run("no-early-stop", 9999L)
