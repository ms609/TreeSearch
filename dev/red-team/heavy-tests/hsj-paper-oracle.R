# HSJ correctness oracle, derived from Hopkins & St John (2021),
# "Incorporating Hierarchical Characters into Phylogenetic Analysis",
# Syst. Biol. 70(6), doi:10.1093/sysbio/syab005.
#
# Standalone red-team artifact (NOT wired into the test suite) backing findings
# T-374/T-375/T-376. Run it against a candidate fix for those; every check here
# is a property the PAPER requires, not a property of the current code.
#
#   Rscript dev/red-team/heavy-tests/hsj-paper-oracle.R [path/to/library]
#
# Provenance of each check is cited inline. Three of the four FAIL at
# 1a94403b / 5cffb18d, which is what T-376 records.
#
# ---------------------------------------------------------------------------
# WHY THE PAPER SETTLES THE ROOTING QUESTION (the reason this file exists)
#
# The area-10 round of 2026-07-28 left open whether HSJ's rooting-dependence
# was intended. It is not, and the paper is explicit enough to close it:
#
#   p.3  "Following the approach of Fitch (1971), for a fixed phylogenetic
#         tree, we extend the character labelings of the leaves to the internal
#         nodes of the tree and compute the MINIMAL score."
#   p.6  "If a tree already has labelings on all the leaves and the internal
#         nodes, then the score of the tree is the sum of the HSJ
#         dissimilarities across all the branches. [...] we keep track of both
#         the possible score when the controlling primary character is present
#         and when it is absent. [...] We continue until we reach the root of
#         the tree and return the MINIMAL score."
#
# The objective is therefore a MINIMUM over internal-node labelings of a sum of
# dissimilarities over branches. The HSJ dissimilarity is SYMMETRIC in its two
# endpoints (d = "the number of nonmatching secondary characters", p.5), the
# branch set of an unrooted tree does not depend on the rooting, and a minimum
# over labelings introduces no orientation. Hence the objective is
# ROOTING-INVARIANT by construction, and any rooting-dependence is a bug.
#
# Corollary for the implementation: the paper prescribes a two-state DP that
# carries both present/absent possibilities and minimises. A single directional
# DELTRAN-style uppass that commits to ONE resolution (what `fitch_label_char`
# does) is neither the paper's algorithm nor guaranteed minimal, and it is the
# source of the alpha-term rooting-dependence measured in T-374.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) .libPaths(c(args[[1]], .libPaths()))
suppressMessages({library(TreeSearch); library(TreeTools)})

pass <- 0L
fail <- 0L
Check <- function(label, ok, detail = "") {
  if (isTRUE(ok)) {
    pass <<- pass + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    fail <<- fail + 1L
    cat(sprintf("  FAIL  %s%s\n", label, if (nzchar(detail)) paste0("  -- ", detail) else ""))
  }
}

# --- Figure 1 matrix ------------------------------------------------------
# "Example matrix has four taxa [...] four primary characters (character 1,
# 6-8), and four secondary characters (2-5, gray italics) that describe
# variation in the first character." (Fig. 1 caption)
figOne <- rbind(
  t1 = c("1", "1", "1", "1", "1", "0", "0", "0"),
  t2 = c("1", "1", "1", "1", "1", "1", "1", "1"),
  t3 = c("1", "0", "0", "0", "0", "1", "1", "1"),
  t4 = c("1", "0", "0", "0", "0", "0", "0", "0")
)
figOneLeft  <- ape::read.tree(text = "((t1,t2),(t3,t4));")
figOneRight <- ape::read.tree(text = "((t1,t4),(t2,t3));")

# TreeSearch's ValidateHierarchy requires the dataset to contain a "-" token,
# and Fig. 1 has no absent taxa ("coded as '0', none shown here"), so the
# absolute-score checks below extend the matrix with two absent taxa. The
# alpha = 0 identity is used for the extended matrix because it needs no
# hand-derived expected value -- see check 2.
extended <- rbind(
  figOne,
  t5 = c("0", "-", "-", "-", "-", "0", "0", "0"),
  t6 = c("0", "-", "-", "-", "-", "1", "1", "1")
)

hierarchy <- CharacterHierarchy("1" = 2:5)

# =========================================================================
cat("\n[1] Figure 1 absolute scores (primaries only)\n")
cat("    Paper Fig. 1 'Prim. only': left = 6, right = 3.\n")
# This check exercises plain Fitch, not HSJ, and is here as a control: if it
# fails, the matrix has been transcribed wrongly and nothing else is meaningful.
primOnly <- MatrixToPhyDat(figOne[, c(1, 6, 7, 8)])
Check("left  primaries-only == 6", isTRUE(all.equal(TreeLength(figOneLeft,  primOnly), 6)))
Check("right primaries-only == 3", isTRUE(all.equal(TreeLength(figOneRight, primOnly), 3)))

# =========================================================================
cat("\n[2] The alpha = 0 identity  (p.6)\n")
cat("    \"the HSJ approach is equivalent to the Fitch approaches when\n")
cat("     alpha = 0.0, as the contributions of the secondary characters are\n")
cat("     ignored\" -- so HSJ(alpha = 0) MUST equal Fitch over the primary\n")
cat("     characters ALONE (controlling primary INCLUDED: at alpha = 0 it\n")
cat("     still contributes 1 per gain/loss, exactly like any primary).\n")
extDs   <- MatrixToPhyDat(extended)
extPrim <- MatrixToPhyDat(extended[, c(1, 6, 7, 8)])
extNoCtl <- MatrixToPhyDat(extended[, c(6, 7, 8)])
extCtl   <- MatrixToPhyDat(extended[, 1, drop = FALSE])

violations <- 0L
droppedControlling <- 0L
for (i in 1:20) {
  tr <- Preorder(as.phylo(i, 6, tipLabels = rownames(extended)))
  hsj0 <- TreeLength(tr, extDs, hierarchy = hierarchy,
                     inapplicable = "hsj", hsj_alpha = 0)
  if (!isTRUE(all.equal(hsj0, TreeLength(tr, extPrim)))) violations <- violations + 1L
  # Discriminator: does the shortfall equal the controlling primary's own
  # Fitch length? If so the controlling primary is contributing NOTHING.
  if (isTRUE(all.equal(hsj0, TreeLength(tr, extNoCtl)))) {
    droppedControlling <- droppedControlling + 1L
  }
}
Check(sprintf("HSJ(alpha=0) == Fitch(primaries) on 20 trees (violations: %d)", violations),
      violations == 0L)
Check(sprintf("controlling primary contributes (dropped on %d/20 trees)", droppedControlling),
      droppedControlling == 0L,
      "HSJ(alpha=0) == Fitch(NON-controlling primaries only), i.e. the controlling primary's gains/losses are silently ignored -- T-376")

# =========================================================================
cat("\n[3] Rooting invariance  (p.3 / p.6, see header derivation)\n")
cat("    The objective is a minimum over internal-node labelings of a sum of\n")
cat("    SYMMETRIC dissimilarities over the branches of an UNROOTED tree, so\n")
cat("    it cannot depend on where the tree is rooted.\n")
base <- Preorder(as.phylo(7, 6, tipLabels = rownames(extended)))

# GUARD AGAINST A VACUOUS PASS. T-374's rooting-dependence lives ENTIRELY in
# the alpha*d/m secondary term (at alpha = 0 the DP contribution is invariant --
# measured 20 x 10 rootings). So an invariance check is only meaningful on a
# matrix where that term actually contributes something. On THIS matrix, while
# T-376 is unfixed, every tip is misclassified absent, the secondaries are
# inapplicable everywhere, and the alpha term is inert -- so the invariance
# check below would pass while testing nothing at all. Refuse to report that
# as a pass.
alphaLive <- TreeLength(base, extDs, hierarchy = hierarchy,
                        inapplicable = "hsj", hsj_alpha = 1) -
             TreeLength(base, extDs, hierarchy = hierarchy,
                        inapplicable = "hsj", hsj_alpha = 0)
if (isTRUE(all.equal(alphaLive, 0))) {
  fail <- fail + 1L
  cat("  INCONCLUSIVE (counted as FAIL)  the alpha term is INERT on this matrix\n")
  cat("        (HSJ(alpha=1) - HSJ(alpha=0) == 0), so rooting invariance cannot\n")
  cat("        be tested here and a green result would be vacuous. Expected while\n")
  cat("        T-376 is unfixed: every tip is misclassified absent, so the\n")
  cat("        secondaries are inapplicable everywhere. Once check [2] passes,\n")
  cat("        this becomes a live test. Independent evidence that it FAILS when\n")
  cat("        the term IS live: a 10-tip matrix gives HSJ alpha=1 -> 24.5/25\n")
  cat("        across its 10 tip-rootings while alpha=0 gives 20 ten times, and\n")
  cat("        XFORM is rooting-dependent on 165/300 random topologies (T-374).\n")
} else {
  for (alpha in c(0, 0.5, 1)) {
    scores <- vapply(rownames(extended), function(tip) {
      TreeLength(RootTree(base, tip), extDs, hierarchy = hierarchy,
                 inapplicable = "hsj", hsj_alpha = alpha)
    }, double(1))
    Check(sprintf("alpha = %.1f invariant across 6 tip-rootings", alpha),
          length(unique(round(scores, 10))) == 1L,
          paste("scores:", paste(unique(round(scores, 6)), collapse = " / ")))
  }
}

# =========================================================================
cat("\n[4] Per-branch contribution bound  (p.5)\n")
cat("    \"a single controlling primary and its associated secondary\n")
cat("     characters can contribute at most 1 per branch (the same as any\n")
cat("     noncontrolling primary)\" -- with 5 branches' worth of controlling\n")
cat("     primary here, HSJ - Fitch(non-controlling) <= n_branches.\n")
nBranch <- nrow(base$edge)
worst <- -Inf
for (i in 1:20) {
  tr <- Preorder(as.phylo(i, 6, tipLabels = rownames(extended)))
  gap <- TreeLength(tr, extDs, hierarchy = hierarchy,
                    inapplicable = "hsj", hsj_alpha = 1) -
         TreeLength(tr, extNoCtl)
  worst <- max(worst, gap)
}
# Same vacuity caveat as [3]: while the controlling primary contributes nothing,
# the gap is 0 and the bound holds trivially. Report that honestly rather than
# banking it as evidence.
if (isTRUE(all.equal(worst, 0))) {
  cat(sprintf("  VACUOUS (not counted)  gap is 0 on all 20 trees, so the <= %d\n", nBranch))
  cat("        bound holds trivially. Becomes meaningful once check [2] passes.\n")
} else {
  Check(sprintf("controlling-primary contribution <= n_branches (%d); worst observed %.4f",
                nBranch, worst),
        worst <= nBranch + 1e-9)
}

cat(sprintf("\n==== %d passed, %d FAILED ====\n", pass, fail))
cat("A candidate fix for T-374/T-375/T-376 must turn every FAIL into a PASS\n")
cat("without changing any PASS.\n")
if (fail > 0L) quit(status = 1L)
