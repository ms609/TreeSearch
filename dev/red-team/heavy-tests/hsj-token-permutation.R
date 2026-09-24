# T-376 repro: the HSJ score is not a function of the data.
#
# Red-team area 10, 2026-07-28 (opus verifier). Ported into the repo from the
# session scratchpad because findings.md cited it by a temp path that gets
# swept, and the skill's own rule is that artifacts backing an OPEN finding are
# live. Standalone; NOT wired into the test suite.
#
#   Rscript dev/red-team/heavy-tests/hsj-token-permutation.R [path/to/library]
#
# WHAT MAKES THIS CONFOUND-FREE, and why the obvious version is not:
# permuting matrix ROWS does not work -- `pick_state`'s tie-break uses
# `tb_mintip` (smallest supporting tip index, ts_hsj.cpp:113-136), so row order
# changes tip numbering and any score change is ambiguous. This permutes only
# the ARBITRARY CONTRAST ROW ORDER: taxa, tip numbering, tree, `levels` and
# every token's state set are held fixed, and `PhyDatToMatrix()` is asserted
# byte-identical on every iteration. So the dataset really is the same dataset.
#
# Sweep 1 has ZERO secondaries (m == 0), so `fitch_label_char` is never
# entered -- any variation there isolates ts_hsj.cpp:220 (the token-index vs
# levels-index comparison in `primary_present`) ALONE, with T-375 provably out
# of the loop. BGS is the control and must return a single value.
#
# Result at 1a94403b / 5cffb18d: sweep 1 gives HSJ 2 OR 3; sweep 2 gives 1, 2
# OR 4; BGS constant in both. The canonical `-01?` ordering scores the maximum,
# so the documented construction is the over-counting one.
#
# See also hsj-paper-oracle.R, which gives the absolute correctness gate from
# Hopkins & St John (2021) and settles that rooting-invariance is REQUIRED.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) .libPaths(c(args[[1]], .libPaths()))
suppressMessages({library(TreeSearch); library(TreeTools)})

# Relabel the (arbitrary) contrast row order ONLY. Taxa, tip numbering, tree,
# levels and every token's state set are untouched -> the dataset is identical.
PermuteTokens <- function(d, perm) {
  at <- attributes(d)
  inv <- order(perm)
  out <- lapply(unclass(d), function(x) inv[x])
  at$allLevels <- at$allLevels[perm]
  at$contrast  <- at$contrast[perm, , drop = FALSE]
  attributes(out) <- at
  out
}

tips <- paste0("t", 1:4)
tree <- Preorder(RootTree(BalancedTree(tips), "t1"))

Sweep <- function(base, hierarchy, tag) {
  d0  <- MatrixToPhyDat(base)
  ref <- PhyDatToMatrix(d0)[tips, , drop = FALSE]
  nTok  <- length(attr(d0, "allLevels"))
  perms <- as.matrix(expand.grid(rep(list(seq_len(nTok)), nTok)))
  perms <- perms[apply(perms, 1, function(r) !anyDuplicated(r)), , drop = FALSE]
  res <- NULL
  for (i in seq_len(nrow(perms))) {
    d <- PermuteTokens(d0, perms[i, ])
    # The load-bearing assertion: same data, different token order.
    stopifnot(identical(PhyDatToMatrix(d)[tips, , drop = FALSE], ref))
    res <- rbind(res, data.frame(
      allLevels = paste(attr(d, "allLevels"), collapse = ""),
      HSJ = TreeLength(d, tree = tree, hierarchy = hierarchy,
                       inapplicable = "hsj", hsj_alpha = 1),
      BGS = TreeLength(tree, d, inapplicable = "bgs")))
  }
  cat("\n===", tag, "===\n")
  cat(sprintf("  %d token orderings of ONE identical dataset\n", nrow(res)))
  cat("  HSJ distinct scores:", paste(sort(unique(res$HSJ)), collapse = ", "),
      if (length(unique(res$HSJ)) > 1L) "  <- BUG: score depends on token order\n" else "\n")
  cat("  BGS distinct scores:", paste(sort(unique(res$BGS)), collapse = ", "),
      if (length(unique(res$BGS)) == 1L) "  <- control OK (single value)\n" else
        "  <- CONTROL BROKEN, disregard this sweep\n")
  print(head(res[order(res$HSJ), ], 4), row.names = FALSE)
  print(tail(res[order(res$HSJ), ], 3), row.names = FALSE)
  invisible(res)
}

# (1) ZERO secondaries: m == 0, so fitch_label_char is never called.
#     Any variation isolates ts_hsj.cpp:220 alone.
b1 <- rbind(t1 = c("-", "1"), t2 = c("0", "0"),
            t3 = c("1", "?"), t4 = c("?", "-"))
s1 <- Sweep(b1, CharacterHierarchy(`2` = integer(0)),
            "1: primary only (isolates ts_hsj.cpp:220)")

# (2) With one secondary: adds fitch_label_char's token bit-encoding (T-375).
b2 <- rbind(t1 = c("-", "1", "0"), t2 = c("0", "1", "1"),
            t3 = c("1", "0", "-"), t4 = c("?", "1", "?"))
s2 <- Sweep(b2, CharacterHierarchy(`2` = 3L),
            "2: primary + secondary (adds ts_hsj.cpp:50-57)")

cat("\n---- verdict ----\n")
bad <- length(unique(s1$HSJ)) > 1L || length(unique(s2$HSJ)) > 1L
if (bad) {
  cat("T-376 REPRODUCES: HSJ is not a function of the data.\n")
  quit(status = 1L)
}
cat("HSJ invariant to token order in both sweeps (T-376 would be fixed).\n")
