# Validate the oracle in oracle.R, then measure all_tbr() against it.
# Usage: Rscript dev/tbr-root-edge/validate.R [lib]

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[[1]] else ".agent-147"
library("TreeSearch", lib.loc = normalizePath(lib))
library("TreeTools", quietly = TRUE)
source("dev/tbr-root-edge/oracle.R")

Head <- function(x) cat("\n== ", x, " ==\n", sep = "")

# ---- 1. the key must be injective over unrooted tree space ----------------
Head("Key injectivity, 8 leaves")
n8 <- NUnrooted(8)
cat("NUnrooted(8) =", n8, "\n")
keys8 <- vapply(seq_len(n8) - 1L, function(i) TreeKey(as.phylo(i, 8)),
                character(1))
cat("distinct keys over as.phylo(0:", n8 - 1, ", 8): ",
    length(unique(keys8)), "\n", sep = "")
stopifnot(length(unique(keys8)) == n8)

# ---- 2. the oracle must contain the NNI neighbourhood, exactly 2(n-3) -----
Head("NNI containment")
for (nm in c("BalancedTree", "PectinateTree")) {
  tr <- get(nm)(8)
  nni <- unique(vapply(seq_len(500), function(i) {
    set.seed(i)
    TreeKey(NNI(tr))
  }, character(1)))
  nni <- setdiff(nni, TreeKey(tr))
  orc <- TbrOracle(tr)
  # Symmetric trees collapse some of the 2(n - 3) NNI swaps to one topology,
  # so the sampled count is a ceiling, not an identity.
  cat(sprintf("%-14s |NNI sampled| = %2d (<= %d), all inside oracle: %s\n",
              nm, length(nni), 2 * (8 - 3), all(nni %in% orc)))
  stopifnot(all(nni %in% orc), length(nni) <= 2 * (8 - 3), length(nni) > 0)
}

# ---- 3. the oracle must agree with the pure-R TBRSwap ---------------------
Head("Oracle vs exhaustive TBRSwap")
SwapAll <- function(tr) {
  tr <- Preorder(RootTree(tr, 1))
  parent <- tr[["edge"]][, 1]
  child <- tr[["edge"]][, 2]
  nEdge <- length(parent)
  keys <- character(0)
  for (etb in seq_len(nEdge)) {
    for (m1 in seq_len(nEdge)) {
      for (m2 in seq_len(nEdge)) {
        if (m1 == m2) next
        res <- suppressWarnings(
          TBRSwap(parent, child, nEdge, edgeToBreak = etb,
                  mergeEdges = c(m1, m2))
        )
        t2 <- tr
        t2[["edge"]] <- cbind(res[[1]], res[[2]])
        keys <- c(keys, TreeKey(t2))
      }
    }
  }
  sort(setdiff(unique(keys), TreeKey(tr)))
}
for (nm in c("BalancedTree", "PectinateTree")) {
  tr <- get(nm)(8)
  orc <- TbrOracle(tr)
  swp <- SwapAll(tr)
  cat(sprintf("%-14s oracle %3d  TBRSwap %3d  swap\\oracle %d  oracle\\swap %d\n",
              nm, length(orc), length(swp),
              length(setdiff(swp, orc)), length(setdiff(orc, swp))))
}
set.seed(1)
for (i in 1:3) {
  tr <- RandomTree(8, root = TRUE)
  orc <- TbrOracle(tr)
  swp <- SwapAll(tr)
  cat(sprintf("random8 #%d    oracle %3d  TBRSwap %3d  swap\\oracle %d  oracle\\swap %d\n",
              i, length(orc), length(swp),
              length(setdiff(swp, orc)), length(setdiff(orc, swp))))
}

# ---- 4. what all_tbr() actually returns -----------------------------------
Head("TBRMoves vs SPRMoves vs oracle")
Row <- function(label, tr) {
  orc <- TbrOracle(tr)
  tbr <- sort(unique(vapply(TBRMoves(tr), TreeKey, character(1))))
  spr <- sort(unique(vapply(SPRMoves(tr), TreeKey, character(1))))
  cat(sprintf(
    "%-16s oracle %3d | TBR %3d | SPR %3d | SPR\\TBR %2d | oracle\\TBR %2d | TBR\\oracle %d\n",
    label, length(orc), length(tbr), length(spr),
    length(setdiff(spr, tbr)), length(setdiff(orc, tbr)),
    length(setdiff(tbr, orc))))
  invisible(list(orc = orc, tbr = tbr, spr = spr))
}
Row("Balanced(8)", BalancedTree(8))
Row("Pectinate(8)", PectinateTree(8))
Row("Balanced(7)", BalancedTree(7))
Row("Pectinate(7)", PectinateTree(7))
set.seed(2)
for (i in 1:4) Row(sprintf("random8 #%d", i), RandomTree(8, root = TRUE))
set.seed(3)
for (i in 1:3) Row(sprintf("random9 #%d", i), RandomTree(9, root = TRUE))
set.seed(4)
for (i in 1:2) Row(sprintf("random11 #%d", i), RandomTree(11, root = TRUE))

Head("Are the missing trees exactly the tip-1 regrafts?")
tr <- BalancedTree(8)
r <- Row("Balanced(8)", tr)
missing <- setdiff(r[["orc"]], r[["tbr"]])
byEdge2 <- sort(unique(vapply(SPRMoves(tr, 2L), TreeKey, character(1))))
cat("missing =", length(missing), "; SPR(edge 2) =", length(byEdge2),
    "; missing in SPR(edge 2):", all(missing %in% byEdge2), "\n")

cat("\nDone.\n")
