# Is exact-cost greedy stepwise addition highly sensitive to TIE-BREAKING?
# Brute-force exact greedy with first-min ('<') vs last-min ('<=') tie-break,
# same addition orders.  If they swing by ~100+, the kernel's 1482-vs-1305 gap
# is tie-break (kernel correct, just a poor deterministic tie-break), and the
# real lever is a good/random tie-break (cf. TNT rseed[).
suppressMessages({ library(TreeSearch); library(TreeTools) })
nm  <- Sys.getenv("DS", "Zanol2014")
phy <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
taxa <- names(phy); n <- length(taxa)

Brute <- function(ord, lastMin = FALSE) {
  tr <- PectinateTree(ord[1:3])
  for (k in 4:n) {
    tip <- ord[k]; phySub <- phy[ord[1:k]]
    nNodeNow <- 2L * (k - 1L) - 1L; best <- NULL; bestLen <- Inf
    for (w in 0:nNodeNow) {
      cand <- tryCatch(AddTip(tr, where = w, label = tip), error = function(e) NULL)
      if (is.null(cand)) next
      L <- TreeLength(cand, phySub)
      take <- if (lastMin) (L <= bestLen) else (L < bestLen)
      if (take) { bestLen <- L; best <- cand }
    }
    tr <- best
  }
  TreeLength(tr, phy)
}

cat(sprintf("== %s | brute greedy tie-break sensitivity ==\n", nm))
for (s in 1:6) {
  set.seed(3000 + s); ord <- sample(taxa)
  cat(sprintf("  seed %d: firstMin=%.0f  lastMin=%.0f  diff=%+.0f\n",
              s, Brute(ord, FALSE), Brute(ord, TRUE),
              Brute(ord, TRUE) - Brute(ord, FALSE)))
}
