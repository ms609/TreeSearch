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
# CORRECTION, 2026-08-02, from reading Algorithm 1 (p.15) against the code.
# The corollary this header used to draw -- "the paper prescribes a two-state DP
# and the code does not implement it" -- was wrong, and cost a round of work.
# `score_hierarchy_block()`'s a(n)/p(n) recurrences ARE Algorithm 1 lines 6-7,
# term for term. The two-state DP was already there and was already
# rooting-invariant: its branch costs are symmetric and it minimises over the
# root's own state. The defect was never in the DP; it was that `d(u, v)` was
# read off `fitch_label_char()`'s directional resolution.
#
# Note also that copying Algorithm 1 literally would NOT have fixed this: its
# line 2 sets L(n) to the first-pass Fitch labelling, which is a downpass and so
# root-dependent, and its line 8 updates L(n) in postorder. Algorithm 1 takes
# "Tree, T, with root r" as input. Theorem 2 claims it returns the minimal
# score, which would make it rooting-invariant; the measurements below are
# evidence against that claim as stated.
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

# Fig. 1's ABSOLUTE HSJ scores: 7 (left) and 5 (right) at alpha = 1.  Re-derived
# as 6 + alpha and 3 + 2 * alpha -- two equations satisfied by one consistent
# alpha, obtained without choosing a root.  Column 9 is the inapplicable-token
# carrier ValidateHierarchy demands; its only non-"1" cell is a lone "-", which
# the Fitch pass scores as costing nothing, so it shifts neither expectation.
figOneExt <- cbind(figOne, c("-", "1", "1", "1"))
figOneDs <- MatrixToPhyDat(figOneExt)
for (alpha in c(0, 0.5, 1)) {
  l <- TreeLength(figOneLeft, figOneDs, hierarchy = hierarchy,
                  inapplicable = "hsj", hsj_alpha = alpha)
  r <- TreeLength(figOneRight, figOneDs, hierarchy = hierarchy,
                  inapplicable = "hsj", hsj_alpha = alpha)
  Check(sprintf("alpha = %.1f: left == %.1f and right == %.1f",
                alpha, 6 + alpha, 3 + 2 * alpha),
        isTRUE(all.equal(l, 6 + alpha)) && isTRUE(all.equal(r, 3 + 2 * alpha)),
        sprintf("got left %s, right %s", format(l), format(r)))
}

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

# Root on the edge above every non-root node; that covers all 2n-3 edges.
AllRootings <- function(tr) {
  out <- list()
  for (v in seq_len(max(tr[["edge"]]))) {
    r <- try(RootOnNode(tr, v, resolveRoot = TRUE), silent = TRUE)
    if (!inherits(r, "try-error")) out[[length(out) + 1L]] <- Preorder(r)
  }
  out
}

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
  cat("  NOTE  the Fig.1-derived matrix below is small and its block is ALL-PRESENT,\n")
  cat("        which is the one regime that was ALREADY invariant before T-374 was\n")
  cat("        fixed (see check [5]). Passing here is therefore necessary but not\n")
  cat("        sufficient -- section [3b] carries the discriminating cases, on\n")
  cat("        blocks with MIXED present/absent primaries, which is where T-374\n")
  cat("        actually lived.\n")
  for (alpha in c(0, 0.5, 1)) {
    scores <- vapply(AllRootings(base), function(rt) {
      TreeLength(rt, extDs, hierarchy = hierarchy,
                 inapplicable = "hsj", hsj_alpha = alpha)
    }, double(1))
    Check(sprintf("alpha = %.1f invariant across all %d edge-rootings",
                  alpha, length(AllRootings(base))),
          length(unique(round(scores, 10))) == 1L,
          paste("scores:", paste(unique(round(scores, 6)), collapse = " / ")))
  }
}

# =========================================================================
cat("\n[3b] Rooting invariance on MIXED blocks -- the discriminating check\n")
cat("    T-374's dependence was confined to blocks with mixed present/absent\n")
cat("    primaries. Two mechanisms, both in the SECONDARY labelling and not in\n")
cat("    the a(n)/p(n) DP (which is already Algorithm 1 lines 6-7, and which\n")
cat("    check [5] shows was already invariant):\n")
cat("      1. '-' was admitted as an ordinary state of a secondary, so a node\n")
cat("         INSIDE the present region could be resolved to it, where it is\n")
cat("         disjoint from every present neighbour in every secondary at once\n")
cat("         and the branch was charged d = m, the full alpha.\n")
cat("      2. The residual resolution was a DELTRAN uppass whose direction, and\n")
cat("         tie-break counts whose subtrees, were properties of the rooting.\n")
cat("    Measured at 93c81a9a over every edge-rooting of 30 random 9-tip trees,\n")
cat("    alpha = 1: 21/30 dependent at m = 2 and 26/30 at m = 4, spread 1.00.\n")
cat("    Crossed over an ambiguous ('?') primary, 18/30 dependent pre-fix, and\n")
cat("    over multistate secondaries, which exercise pick_state()'s tie-break.\n")
set.seed(374)
for (secStates in list(c("0", "1"), c("0", "1", "2"))) {
  for (ambiguous in c(FALSE, TRUE)) {
    mixedDep <- 0L
    mixedWorst <- 0
    for (rep in 1:8) {
      nTip <- 9L
      pri <- rep("1", nTip)
      pri[sample.int(nTip, 3L)] <- "0"
      if (ambiguous) pri[sample(which(pri == "1"), 2L)] <- "?"
      live <- pri != "0"
      sec <- matrix("-", nTip, 3L)
      for (j in 1:3) sec[live, j] <- sample(secStates, sum(live), TRUE)
      nh <- matrix(sample(c("0", "1"), nTip * 3L, TRUE), nTip, 3L)
      mm <- cbind(pri, sec, nh)
      rownames(mm) <- paste0("t", seq_len(nTip))
      colnames(mm) <- NULL
      mDs <- MatrixToPhyDat(mm)
      mH <- CharacterHierarchy("1" = 2:4)
      tr <- Preorder(as.phylo(rep, nTip, tipLabels = rownames(mm)))
      sc <- vapply(AllRootings(tr), function(rt)
        TreeLength(rt, mDs, hierarchy = mH, inapplicable = "hsj", hsj_alpha = 1),
        double(1))
      if (diff(range(sc)) > 1e-9) {
        mixedDep <- mixedDep + 1L
        mixedWorst <- max(mixedWorst, diff(range(sc)))
      }
    }
    Check(sprintf("mixed, %d secondary states, '?' primary %-5s: %d/8 dependent, worst %.4f",
                  length(secStates), ambiguous, mixedDep, mixedWorst),
          mixedDep == 0L)
  }
}

# A secondary is freed only where the primary CANNOT be present, not merely
# where it MAY be absent: observing a secondary is evidence the structure is
# present, and freeing it at every "?" primary would discard that and leave a
# block of all-"?" primaries with an empty domain and a silently zero alpha
# term.  Mirrors recode_hierarchy.R's tipSecKnown path (T-379).
qBase <- matrix(c(
  "1",  "0",  "0",  "1",  "1",
  "1",  "0",  "0",  "1",  "0",
  "?",  "0",  "0",  "1",  "0",
  "0",  "-",  "-",  "0",  "0",
  "1",  "1",  "1",  "1",  "1",
  "1",  "1",  "1",  "0",  "1"
), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
qH <- CharacterHierarchy("1" = 2:3)
qTr <- Preorder(ape::read.tree(text = "(t1,((t2,t3),(t4,(t5,t6))));"))
qSc <- vapply(c("0", "1"), function(v) {
  mq <- qBase; mq[3, 2:3] <- v
  TreeLength(qTr, MatrixToPhyDat(mq), hierarchy = qH, inapplicable = "hsj",
             hsj_alpha = 1)
}, double(1))
Check(sprintf("an OBSERVED secondary at a '?' primary still counts (%s vs %s)",
              format(qSc[[1]]), format(qSc[[2]])),
      !isTRUE(all.equal(qSc[[1]], qSc[[2]])))

# =========================================================================
cat("\n[5] All-present closed form -- the regression FLOOR\n")
cat("    Under ANY most-parsimonious reconstruction of secondary j, the number\n")
cat("    of branches on which j changes is FitchLen_j. So when every node in\n")
cat("    the block is present the alpha term is (alpha/m) * sum_j FitchLen_j\n")
cat("    EXACTLY, whichever labelling the uppass picks -- which is why that\n")
cat("    regime was already rooting-invariant, and why T-374 could only ever\n")
cat("    surface on mixed blocks. A fix that breaks this is wrong regardless\n")
cat("    of what it does to the rooting spread.\n")
set.seed(3741)
cfBad <- 0L
cfWorst <- 0
for (rep in 1:8) {
  nTip <- 9L
  nSec <- 4L
  sec <- matrix(sample(c("0", "1"), nTip * nSec, TRUE), nTip, nSec)
  nh <- matrix(sample(c("0", "1"), nTip * 3L, TRUE), nTip, 3L)
  nh[1, 1] <- "-"   # inapplicable-token carrier ValidateHierarchy demands
  pm <- cbind(rep("1", nTip), sec, nh)
  rownames(pm) <- paste0("t", seq_len(nTip))
  colnames(pm) <- NULL
  pDs <- MatrixToPhyDat(pm)
  pH <- CharacterHierarchy("1" = 2:(nSec + 1L))
  tr <- Preorder(as.phylo(rep, nTip, tipLabels = rownames(pm)))
  fitchPri <- TreeLength(tr, MatrixToPhyDat(
    pm[, c(1L, (nSec + 2L):ncol(pm)), drop = FALSE]))
  secLen <- sum(vapply(2:(nSec + 1L), function(j)
    TreeLength(tr, MatrixToPhyDat(pm[, j, drop = FALSE])), double(1)))
  for (alpha in c(0.5, 1)) {
    got <- TreeLength(tr, pDs, hierarchy = pH, inapplicable = "hsj",
                      hsj_alpha = alpha)
    dev <- abs(got - (fitchPri + alpha * secLen / nSec))
    if (dev > 1e-9) { cfBad <- cfBad + 1L; cfWorst <- max(cfWorst, dev) }
  }
}
Check(sprintf("closed form holds on all-present blocks (%d/16 violations, worst %.4f)",
              cfBad, cfWorst),
      cfBad == 0L)

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
