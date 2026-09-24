# tbr_oracle_na.R -- GROUND-TRUTH completeness oracle on a REAL (inapplicable)
# dataset, scored by the kernel's own ts_fitch_score (NA 3-pass, with concavity).
# For each of several random starts, run the kernel TBR to convergence (direct
# OR physical), then assert no TBR/SPR neighbour (TreeTools enumerators, two
# rootings) scores strictly lower.  This does NOT assume physical is complete --
# it is the independent oracle that decides whether EITHER path reaches true
# unrooted-TBR optima for NA.
#
# Usage: Rscript dev/benchmarks/tbr_oracle_na.R [dataset] [concavity] [phys 0/1] [nStart]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"), winslash = "/"))
  library(TreeTools)
})
a <- commandArgs(trailingOnly = TRUE)
dsn  <- if (length(a) >= 1) a[[1]] else "Aria2015"
conc <- if (length(a) >= 2) as.numeric(a[[2]]) else -1
phys <- if (length(a) >= 3) as.logical(as.integer(a[[3]])) else FALSE
nStart <- if (length(a) >= 4) as.integer(a[[4]]) else 12L
eps <- 1e-6
data("inapplicable.phyData", package = "TreeSearch")
dataset <- inapplicable.phyData[[dsn]]
at <- attributes(dataset)
d <- list(contrast = at$contrast,
          tip_data = matrix(unlist(dataset, use.names = FALSE), nrow = length(dataset), byrow = TRUE),
          weight = at$weight, levels = at$levels, labels = names(dataset))
scoreTree <- function(tree) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  TreeSearch:::ts_fitch_score(edge, d$contrast, d$tip_data, d$weight, d$levels, concavity = conc)
}
kernelTbr <- function(tree) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  if (phys) Sys.setenv(TS_PHYS_REROOT = "1") else Sys.unsetenv("TS_PHYS_REROOT")
  res <- TreeSearch:::ts_tbr_diagnostics(edge, d$contrast, d$tip_data, d$weight, d$levels,
           maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, concavity = conc, unrooted = TRUE)
  Sys.unsetenv("TS_PHYS_REROOT")
  structure(list(edge = res$edge, Nnode = length(d$labels) - 1L, tip.label = d$labels), class = "phylo")
}
bestImproving <- function(tree) {
  base <- scoreTree(tree); cand <- list()
  for (rt in d$labels[1:2]) { rr <- Preorder(RootTree(tree, rt)); cand <- c(cand, TBRMoves(rr), SPRMoves(rr)) }
  ls <- vapply(cand, scoreTree, double(1))
  if (min(ls) < base - eps) list(len = min(ls), base = base) else NULL
}
cat(sprintf("=== NA oracle: %s conc=%s path=%s, %d starts ===\n", dsn, conc, if (phys) "PHYSICAL" else "DIRECT", nStart))
fails <- 0L; ff <- NULL
for (i in seq_len(nStart)) {
  set.seed(7000L + i); start <- RandomTree(dataset, root = TRUE)
  set.seed(i); conv <- kernelTbr(start)
  imp <- bestImproving(conv)
  if (!is.null(imp)) { fails <- fails + 1L; if (is.null(ff)) ff <- list(i = i, conv = scoreTree(conv), imp = imp) }
}
cat(sprintf("RESULT: %d / %d converged trees had an improving neighbour (FAILURES)\n", fails, nStart))
if (fails == 0L) cat("=> COMPLETE on this sample (kernel-scored).\n") else
  cat(sprintf("=> INCOMPLETE. First: start #%d, converged %.4f, neighbour %.4f (miss %.4f)\n",
              ff$i, ff$conv, ff$imp$len, ff$conv - ff$imp$len))
