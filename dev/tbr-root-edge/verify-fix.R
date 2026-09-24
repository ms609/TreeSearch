# Post-fix verification for agent-issues/TreeSearch#147.
# Usage: Rscript dev/tbr-root-edge/verify-fix.R [lib]

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[[1]] else ".agent-147"
library("TreeSearch", lib.loc = normalizePath(lib))
library("TreeTools", quietly = TRUE)
source("dev/tbr-root-edge/oracle.R")

Head <- function(x) cat("\n== ", x, " ==\n", sep = "")
Keys <- function(trees) sort(unique(vapply(trees, TreeKey, character(1))))

# ---- 1. TBRMoves == the complete TBR neighbourhood ------------------------
Head("TBRMoves vs oracle, and SPRMoves containment")
Check <- function(label, tr) {
  orc <- TbrOracle(tr)
  tbr <- Keys(TBRMoves(tr))
  spr <- Keys(SPRMoves(tr))
  ok <- identical(orc, tbr) && !length(setdiff(spr, tbr))
  cat(sprintf("%-16s oracle %4d | TBR %4d | missing %2d | spurious %2d | SPR\\TBR %2d  %s\n",
              label, length(orc), length(tbr),
              length(setdiff(orc, tbr)), length(setdiff(tbr, orc)),
              length(setdiff(spr, tbr)), if (ok) "OK" else "**FAIL**"))
  ok
}
ok <- c(
  Check("Balanced(7)", BalancedTree(7)),
  Check("Pectinate(7)", PectinateTree(7)),
  Check("Balanced(8)", BalancedTree(8)),
  Check("Pectinate(8)", PectinateTree(8)),
  Check("Balanced(9)", BalancedTree(9)),
  Check("Pectinate(9)", PectinateTree(9))
)
set.seed(2)
for (i in 1:4) ok <- c(ok, Check(sprintf("random8 #%d", i), RandomTree(8, root = TRUE)))
set.seed(3)
for (i in 1:3) ok <- c(ok, Check(sprintf("random9 #%d", i), RandomTree(9, root = TRUE)))
set.seed(4)
for (i in 1:2) ok <- c(ok, Check(sprintf("random11 #%d", i), RandomTree(11, root = TRUE)))
stopifnot(all(ok))

# ---- 2. exhaustive over every 7-leaf unrooted topology --------------------
Head("Every unrooted 7-leaf topology")
n7 <- NUnrooted(7)
bad <- 0L
for (i in seq_len(n7) - 1L) {
  tr <- as.phylo(i, 7)
  if (!identical(TbrOracle(tr), Keys(TBRMoves(tr)))) bad <- bad + 1L
}
cat(sprintf("%d / %d topologies match the oracle exactly\n", n7 - bad, n7))
stopifnot(bad == 0L)

# ---- 3. the root-edge identity: TBR on edge 2 == SPR on edge 2 ------------
Head("all_tbr(e, 2) == all_spr(e, 2)")
for (n in 5:12) {
  set.seed(n)
  tr <- Preorder(RootTree(RandomTree(n, root = TRUE), 1))
  e <- tr[["edge"]]
  a <- TreeSearch:::all_tbr(e, 2L)
  b <- TreeSearch:::.all_spr(e, 2L)
  # Theory: bisecting tip 1's pendant edge leaves 2(n-1)-3 re-rooting sites,
  # one of which recreates the starting tree.
  expected <- 2 * (n - 1) - 3 - 1
  cat(sprintf("n=%2d  tbr %2d  spr %2d  identical %s  expected %2d %s\n",
              n, length(a), length(b), identical(a, b), expected,
              if (length(a) == expected) "OK" else "**FAIL**"))
  stopifnot(identical(a, b), length(a) == expected)
}

cat("\nAll checks passed.\n")
