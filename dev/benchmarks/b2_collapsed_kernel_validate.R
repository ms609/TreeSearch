# Validate the KERNEL's aggressive collapsed flags against a brute-force MPR
# oracle (2026-06-22).  The criterion (final[p] & final[c] != 0) was already
# validated in R (b2_minlength_oracle.R, 0/3192); this confirms the C++ port —
# final_ freshness for EW + the transposed-bitset reduction — end-to-end.
#
# Method: small random tree + random data -> kernel ts_collapsed_flags_debug
# returns 0-based parent[] + per-node flags -> run a GENERIC brute-force oracle
# on the kernel's OWN parent[] (no numbering translation) -> compare.  A length
# cross-check (oracle min summed over chars == Fitch TreeLength) guards tip
# alignment (na-validation-alignment-gotcha).
#
# Usage: TS_LIB=.agent-fuse Rscript dev/benchmarks/b2_collapsed_kernel_validate.R

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-fuse"),
                                              winslash = "/", mustWork = TRUE))
  library(TreeTools)
})

nState <- 3L

# Build a USER phyDat + kernel data bundle from a tip x char matrix.
bundleFromMatrix <- function(m) {
  tips <- paste0("t", seq_len(nrow(m)))
  rownames(m) <- tips
  phy <- phangorn::phyDat(m, type = "USER", levels = as.character(0:(nState - 1L)))
  at <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, labels = names(phy), m = m[names(phy), , drop = FALSE])
}

# Generic brute-force soft-min-length-0 on kernel 0-based parent[] (root = n_tip,
# self-parent).  tipState indexed by kernel tip id (0..n_tip-1) = row of mAligned.
oracleKernel <- function(parent0, nTip, nNode, mAligned) {
  root <- nTip                              # kernel root id (0-based)
  isTip <- function(v) v < nTip
  internal <- setdiff(0:(nNode - 1L), 0:(nTip - 1L))   # nTip .. nNode-1 (incl root)
  grid <- as.matrix(expand.grid(rep(list(0:(nState - 1L)), length(internal))))
  nonRootEdges <- setdiff(0:(nNode - 1L), root)        # every node except root has a parent edge
  zero <- rep(TRUE, nNode); names(zero) <- 0:(nNode - 1L)
  Msum <- 0
  for (cc in seq_len(ncol(mAligned))) {
    tipState <- mAligned[, cc]
    cost <- function(g) {
      st <- integer(nNode)
      st[1:nTip] <- tipState
      st[internal + 1L] <- g
      ch <- 0L
      for (v in 0:(nNode - 1L)) {
        p <- parent0[v + 1L]
        if (v != root && p != v && st[v + 1L] != st[p + 1L]) ch <- ch + 1L
      }
      ch
    }
    nc <- apply(grid, 1L, cost); M <- min(nc); Msum <- Msum + M
    optRows <- which(nc == M)
    for (c0 in nonRootEdges) {
      p0 <- parent0[c0 + 1L]
      feas <- FALSE
      for (r in optRows) {
        st <- integer(nNode); st[1:nTip] <- tipState; st[internal + 1L] <- grid[r, ]
        if (st[c0 + 1L] == st[p0 + 1L]) { feas <- TRUE; break }
      }
      if (!feas) zero[c0 + 1L] <- FALSE
    }
  }
  list(zero = zero, Msum = Msum)
}

# A character is parsimony-informative iff >=2 states each appear >=2 times.
# The kernel simplifies away uninformative characters (constants AND
# autapomorphies), which would make terminal-edge collapsibility differ from a
# raw-column oracle.  Generate informative-only columns so simplification is a
# no-op and the two views see identical characters.
informativeCol <- function(col) sum(table(col) >= 2L) >= 2L
infoCol <- function(nt) { repeat { c <- sample(0:(nState - 1L), nt, TRUE); if (informativeCol(c)) return(c) } }

set.seed(7)
tot <- 0L; mism <- 0L; fp <- 0L; fn <- 0L; lenBad <- 0L
totInt <- 0L; mismInt <- 0L
for (trial in seq_len(120L)) {
  nt <- sample(4:6, 1L); nChar <- sample(4:7, 1L)
  m <- vapply(seq_len(nChar), function(j) infoCol(nt), integer(nt))
  b <- bundleFromMatrix(m)
  set.seed(1000L + trial)
  tr <- Preorder(RenumberTips(RandomTree(b$phy, root = TRUE), b$labels))
  edge <- tr[["edge"]]
  res <- TreeSearch:::ts_collapsed_flags_debug(
    edge, b$contrast, b$tip_data, b$weight, b$levels, aggressive = TRUE)
  orc <- oracleKernel(res$parent, res$n_tip, res$n_node, b$m)
  # length cross-check: oracle min summed over chars (unweighted; all weight 1
  # for USER data here) must equal Fitch tree length
  tl <- TreeLength(tr, b$phy)
  if (abs(orc$Msum - tl) > 0.5) lenBad <- lenBad + 1L
  # compare flags on every non-root edge (kernel sets root-children to 0; oracle
  # may mark a root-child collapsible — exclude root-children to match the
  # kernel's deliberate skip, as the conservative path also skips them)
  for (c0 in 0:(res$n_node - 1L)) {
    if (c0 == res$n_tip) next                      # root
    if (res$parent[c0 + 1L] == res$n_tip) next      # root-child (kernel skips by design)
    k <- res$collapsed[c0 + 1L] == 1L
    o <- orc$zero[c0 + 1L]
    isInternal <- c0 >= res$n_tip
    tot <- tot + 1L
    if (isInternal) totInt <- totInt + 1L
    if (k != o) { mism <- mism + 1L; if (k && !o) fp <- fp + 1L else fn <- fn + 1L
      if (isInternal) mismInt <- mismInt + 1L }
  }
}
cat(sprintf("Validated %d non-root-child edges over 120 kernel trees (nt 4-6, informative chars)\n", tot))
cat(sprintf("  length cross-check failures (alignment): %d / 120 trees\n", lenBad))
cat(sprintf("  kernel-vs-oracle mismatches: %d (fp=%d over-collapse, fn=%d under-collapse)\n",
            mism, fp, fn))
cat(sprintf("  internal-edge-only mismatches: %d / %d\n", mismInt, totInt))
# Gate on INTERNAL edges: the aggressive criterion flags internal branches only
# (TNT collapse-3 semantics; terminal/pendant edges are never collapsed), so the
# internal-edge match is the verdict.  Terminal-edge "mismatches" are the
# deliberately-unflagged tips (under-collapse) + the autapomorphy-simplification
# confound, neither a correctness concern.
cat(if (mismInt == 0L && lenBad == 0L)
      "=> KERNEL aggressive flags MATCH the oracle on ALL internal edges. Port verified.\n"
    else "=> INTERNAL-EDGE MISMATCH — investigate before trusting the kernel flags.\n")
