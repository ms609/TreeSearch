# Isolate "union formula wrong" vs "incremental final_ stale": run the
# directional kernel with TS_WAGNER_UNION (production union-of-finals formula
# + FULL fitch_score each step) vs my directional combine, vs brute-force
# oracle, all on the SAME addition order.
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
})
nm   <- Sys.getenv("DS", "Zanol2014")
phy  <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
taxa <- names(phy); n <- length(taxa)
at <- attributes(phy); contrast <- at$contrast
tipData <- matrix(unlist(phy, use.names = FALSE), nrow = n, byrow = TRUE)
weight  <- TreeSearch:::.ScaleWeight(at$weight); levels <- at$levels

Dir <- function(ord) {
  TreeSearch:::ts_wagner_tree_dir(contrast, tipData, weight, levels,
                                  addition_order = as.integer(match(ord, taxa)))$score
}
Brute <- function(ord) {
  tr <- PectinateTree(ord[1:3])
  for (k in 4:n) {
    tip <- ord[k]; phySub <- phy[ord[1:k]]
    nNodeNow <- 2L * (k - 1L) - 1L; best <- NULL; bestLen <- Inf
    for (w in 0:nNodeNow) {
      cand <- tryCatch(AddTip(tr, where = w, label = tip), error = function(e) NULL)
      if (is.null(cand)) next
      L <- TreeLength(cand, phySub); if (L < bestLen) { bestLen <- L; best <- cand }
    }
    tr <- best
  }
  TreeLength(tr, phy)
}

cat(sprintf("== %s | union(full final_) vs combine vs brute, same order ==\n", nm))
for (s in 1:6) {
  set.seed(2000 + s); ord <- sample(taxa)
  Sys.setenv(TS_WAGNER_UNION = "1"); u <- Dir(ord)
  Sys.unsetenv("TS_WAGNER_UNION"); d <- Dir(ord)
  b <- Brute(ord)
  cat(sprintf("  seed %d: union=%.0f  combine=%.0f  brute=%.0f\n", s, u, d, b))
}
