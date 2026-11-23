# Validate the directional edge-set Wagner kernel (ts_wagner_tree_dir) against
# the exact-insertion oracle (brute-force full-rescore argmin).  Strict gate:
#   (1) score DISTRIBUTION matches the oracle on >=2 datasets;
#   (2) on the SAME addition order, directional == brute-force (modulo early
#       tie-break divergence). A consistently-higher directional score would
#       indicate a close-but-wrong formula.
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
})
datasets <- strsplit(Sys.getenv("DS", "Zanol2014 Zhu2013"), "\\s+")[[1]]

Mats <- function(phy) {
  n <- length(phy); at <- attributes(phy)
  list(n = n, contrast = at$contrast,
       tipData = matrix(unlist(phy, use.names = FALSE), nrow = n, byrow = TRUE),
       weight = TreeSearch:::.ScaleWeight(at$weight), levels = at$levels)
}

# Brute-force exact-insertion Wagner for an EXPLICIT addition order.
BruteWagner <- function(phy, ord) {
  taxa <- names(phy); n <- length(taxa)
  tr <- PectinateTree(ord[1:3])
  for (k in 4:n) {
    tip <- ord[k]; phySub <- phy[ord[1:k]]
    nNodeNow <- 2L * (k - 1L) - 1L
    best <- NULL; bestLen <- Inf
    for (w in 0:nNodeNow) {
      cand <- tryCatch(AddTip(tr, where = w, label = tip), error = function(e) NULL)
      if (is.null(cand)) next
      L <- TreeLength(cand, phySub)
      if (L < bestLen) { bestLen <- L; best <- cand }
    }
    tr <- best
  }
  TreeLength(tr, phy)
}

for (nm in datasets) {
  phy  <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
  taxa <- names(phy); M <- Mats(phy)
  cat(sprintf("\n================ %s (n=%d) ================\n", nm, M$n))

  Dir <- function(ord = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    ao <- if (is.null(ord)) integer(0) else as.integer(match(ord, taxa))
    TreeSearch:::ts_wagner_tree_dir(M$contrast, M$tipData, M$weight, M$levels,
                                    addition_order = ao)$score
  }
  Buggy <- function(seed) {
    set.seed(seed)
    TreeSearch:::ts_random_wagner_tree(M$contrast, M$tipData, M$weight, M$levels)$score
  }

  # (1) distribution: directional random vs buggy random
  dirRand <- vapply(1:12, function(s) Dir(seed = s), double(1))
  buggy   <- vapply(1:12, function(s) Buggy(s),       double(1))
  cat(sprintf("  directional(random) mean=%.1f sd=%.1f min=%.0f max=%.0f\n",
              mean(dirRand), sd(dirRand), min(dirRand), max(dirRand)))
  cat(sprintf("  buggy union-formula mean=%.1f sd=%.1f min=%.0f max=%.0f\n",
              mean(buggy), sd(buggy), min(buggy), max(buggy)))

  # oracle reference (3 brute-force trees)
  oracle <- vapply(1:3, function(s) { set.seed(s); BruteWagner(phy, sample(taxa)) }, double(1))
  cat(sprintf("  oracle brute-force  : %s\n", paste(round(oracle), collapse = " ")))

  # (2) same-order: directional vs brute-force on identical addition orders
  cat("  -- same-order (directional vs brute-force) --\n")
  for (s in 1:6) {
    set.seed(1000 + s); ord <- sample(taxa)
    d <- Dir(ord = ord); b <- BruteWagner(phy, ord)
    cat(sprintf("    seed %d: directional=%.0f  brute=%.0f  diff=%+.0f%s\n",
                s, d, b, d - b, if (d > b) "  <-- dir WORSE" else ""))
  }
}
