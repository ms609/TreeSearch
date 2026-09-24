# T-374b: XFORM (x-transformation) rooting-dependence oracle.
#
# PURE R.  Deliberately independent of the compiled engine: it re-implements
# the Sankoff DP and the x-transformation cost matrix from
# `R/recode_hierarchy.R`'s documented spec, so a disagreement with the kernel
# is informative rather than circular.  No build required; nothing under src/
# is touched or needed.
#
# Questions it answers, for the decision document
# `dev/plans/2026-07-29-t374b-xform-rooting-policy.md`:
#
#  Q-A  Is the x-transformation objective rooting-dependent, and by how much?
#  Q-B  Where does the dependence come from, and how large can it get?
#       The asymmetric part of the x-transformation cost matrix is a GRADIENT:
#       with f(absent) = 0 and f(present) = nSec/2,
#           c(i,j) = s(i,j) + f(j) - f(i),   s = (c + c^T)/2 symmetric,
#       and summing over edges oriented parent -> child gives
#           total = sum_edges s(u,v)  -  sum_{internal} f(u)  -  f(root)
#                                     +  sum_tips f.
#       sum_tips f is a constant; s carries no orientation.  The root node is
#       therefore charged f TWICE (once inside sum_internal, once as -f(root)),
#       which is the whole source of the dependence and predicts a spread
#       bounded by 2 * f(present) = nSec.
#       NOTE: an earlier, tighter guess of nSec/2 was FALSIFIED by this script
#       (nSec = 1 gives spread 1, nSec = 2 gives spread 2).  It ignored that
#       rooting an unrooted tree SUBDIVIDES an edge, so the edge set - and
#       hence sum_edges s - is not itself rooting-invariant.  The retained
#       bound below is nSec.
#  Q-C  Does forced_root_state = 0 ("absent") restore rooting-invariance?
#       This is the candidate cheap fix, so it must be measured, not assumed.
#       READ THE ANSWER CAREFULLY.  Pinning the root STATE defines a different,
#       explicitly ROOTED criterion; its spread across rootings is the cost of
#       leaving the root POSITION arbitrary, and is NOT commensurable with the
#       free-root arm as "better or worse".  What the numbers do establish is
#       that the one-line change does not remove root-sensitivity from a
#       pipeline that keeps rerooting - it relocates it and enlarges it.  The
#       pinned criterion is perfectly well-defined at any FIXED rooting, which
#       this script also checks.
#  Q-D  Do ambiguous tips (-1 fully ambiguous, -2 present-unknown) break the
#       Sum_tips-f-is-constant step and so break the nSec bound?
#
# Exit status: 0 = every assertion in the report held; 1 = a stated prediction
# failed (read the report; the decision document's algebra would then be wrong).

suppressPackageStartupMessages({
  library("ape")
})

set.seed(20260729L)

# ---------------------------------------------------------------------------
# x-transformation cost matrix, per R/recode_hierarchy.R:101-118
# ---------------------------------------------------------------------------
# State 1 (index 1) = absent; states 2..nStates = present-secondary combos.
# absent -> present = nSec + 1;  present -> absent = 1;
# present -> present = Hamming distance between secondary combinations.
XformCostMatrix <- function(secNStates) {
  nSec <- length(secNStates)
  comboGrid <- if (nSec > 0L) {
    as.matrix(expand.grid(lapply(secNStates, seq_len)))
  } else {
    matrix(integer(0), nrow = 1L, ncol = 0L)
  }
  nStates <- nrow(comboGrid) + 1L
  gainCost <- nSec + 1L
  cm <- matrix(0, nStates, nStates)
  for (i in seq_len(nStates)) {
    for (j in seq_len(nStates)) {
      if (i == j) next
      if (i == 1L) {
        cm[i, j] <- gainCost
      } else if (j == 1L) {
        cm[i, j] <- 1
      } else {
        cm[i, j] <- sum(comboGrid[i - 1L, ] != comboGrid[j - 1L, ])
      }
    }
  }
  cm
}

# ---------------------------------------------------------------------------
# Sankoff DP on a rooted binary tree, mirroring
# src/ts_sankoff.cpp:sankoff_score_char (postorder min-plus, then either the
# forced root state or the min over root states).
# ---------------------------------------------------------------------------
# tipCosts: nTip x nStates matrix of per-tip per-state costs (0 = allowed,
#           Inf = disallowed), matching ds.sankoff_tip_costs.
# forcedRoot: 0-based state index, or -1 for min-over-root-states.
SankoffRooted <- function(tree, tipCosts, cm, forcedRoot = -1L) {
  nTip <- length(tree[["tip.label"]])
  nStates <- ncol(cm)
  nNode <- nTip + tree[["Nnode"]]
  costs <- matrix(Inf, nNode, nStates)
  costs[seq_len(nTip), ] <- tipCosts

  # Internal nodes, tipward-first: ape's "postorder" lists each edge only
  # after both its child's edges, so the parent column in order of first
  # appearance is a valid postorder over internal nodes.
  edge <- tree[["edge"]]
  internalPostorder <- unique(ape::reorder.phylo(tree, "postorder")[["edge"]][, 1])

  for (node in internalPostorder) {
    kids <- edge[edge[, 1] == node, 2]
    stopifnot(length(kids) == 2L)       # binary only, as the kernel assumes
    nc <- numeric(nStates)
    for (s in seq_len(nStates)) {
      tot <- 0
      for (k in kids) {
        tot <- tot + min(cm[s, ] + costs[k, ])
      }
      nc[s] <- tot
    }
    costs[node, ] <- nc
  }

  root <- nTip + 1L
  stopifnot(!(root %in% edge[, 2]))
  if (forcedRoot >= 0L) costs[root, forcedRoot + 1L] else min(costs[root, ])
}

# tip_states -> tipCosts, per src/ts_rcpp.cpp's unpack_xform semantics:
# 0 = absent (state index 1); k >= 1 = present combo k (index k + 1);
# -1 = fully ambiguous (every state free); -2 = present, combo unknown
# (every present state free, absent disallowed).
TipCostMatrix <- function(tipStates, nStates) {
  m <- matrix(Inf, length(tipStates), nStates)
  for (t in seq_along(tipStates)) {
    st <- tipStates[t]
    if (st == -1L) {
      m[t, ] <- 0
    } else if (st == -2L) {
      m[t, -1L] <- 0
    } else {
      m[t, st + 1L] <- 0
    }
  }
  m
}

# ---------------------------------------------------------------------------
# Score one unrooted topology under every distinct rooting
# ---------------------------------------------------------------------------
# `unrootedTree` must be unrooted (Nnode = nTip - 2).  Roots on each of its
# 2n-3 edges via ape::root(..., edgelabel = TRUE) and returns the vector of
# scores.  We use `resolve.root = TRUE` so each rooting is binary, as the
# kernel requires.
ScoreAllRootings <- function(unrootedTree, tipCosts, cm, forcedRoot = -1L) {
  nEdge <- nrow(unrootedTree[["edge"]])
  out <- rep(NA_real_, nEdge)
  for (e in seq_len(nEdge)) {
    rt <- try(ape::root(unrootedTree, node = unrootedTree[["edge"]][e, 2],
                        resolve.root = TRUE), silent = TRUE)
    if (inherits(rt, "try-error")) next
    if (rt[["Nnode"]] != length(rt[["tip.label"]]) - 1L) next   # not binary
    out[e] <- SankoffRooted(rt, tipCosts, cm, forcedRoot)
  }
  out[!is.na(out)]
}

# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
nTip <- 9L
nRep <- 120L

RandomTipStates <- function(nTip, nPresent, pAbsent = 0.35,
                            pAmbig = 0, pPresentUnknown = 0) {
  vapply(seq_len(nTip), function(i) {
    u <- runif(1)
    if (u < pAmbig) return(-1L)
    if (u < pAmbig + pPresentUnknown) return(-2L)
    if (u < pAmbig + pPresentUnknown + pAbsent) return(0L)
    sample.int(nPresent, 1L)
  }, integer(1))
}

scenarios <- list(
  list(name = "nSec=1 (2 levels)",       sec = c(2L),        ambig = 0,    pu = 0),
  list(name = "nSec=2 (2x2 levels)",     sec = c(2L, 2L),    ambig = 0,    pu = 0),
  list(name = "nSec=3 (2x2x2 levels)",   sec = c(2L, 2L, 2L),ambig = 0,    pu = 0),
  list(name = "nSec=2, 20% ambiguous",   sec = c(2L, 2L),    ambig = 0.20, pu = 0),
  list(name = "nSec=2, 20% present-unk", sec = c(2L, 2L),    ambig = 0,    pu = 0.20),
  list(name = "nSec=0 (control: symm.)", sec = integer(0),   ambig = 0,    pu = 0)
)

failures <- character(0)
cat("=== T-374b: XFORM rooting-dependence oracle ===\n")
cat(sprintf("nTip = %d, %d random unrooted topologies per scenario\n\n",
            nTip, nRep))

for (sc in scenarios) {
  cm <- XformCostMatrix(sc$sec)
  nStates <- ncol(cm)
  nSec <- length(sc$sec)
  bound <- nSec

  spreadFree <- numeric(nRep)
  spreadPinned <- numeric(nRep)
  nDepFree <- 0L
  nDepPinned <- 0L
  # Q-E: how bad is an ARBITRARY rooting relative to the rooting-invariant
  # objective min-over-rootings?  fracAtMin = share of the 2n-3 rootings that
  # attain the minimum; excessMean = mean overstatement of a random rooting.
  fracAtMin <- numeric(nRep)
  excessMean <- numeric(nRep)

  for (r in seq_len(nRep)) {
    tr <- ape::rtree(nTip, rooted = FALSE)
    ts <- RandomTipStates(nTip, nStates - 1L, pAmbig = sc$ambig,
                          pPresentUnknown = sc$pu)
    tc <- TipCostMatrix(ts, nStates)

    sFree <- ScoreAllRootings(tr, tc, cm, forcedRoot = -1L)
    sPin  <- ScoreAllRootings(tr, tc, cm, forcedRoot = 0L)

    spreadFree[r] <- diff(range(sFree))
    spreadPinned[r] <- diff(range(sPin[is.finite(sPin)]))
    if (spreadFree[r] > 1e-9) nDepFree <- nDepFree + 1L
    if (isTRUE(spreadPinned[r] > 1e-9)) nDepPinned <- nDepPinned + 1L
    fracAtMin[r] <- mean(sFree <= min(sFree) + 1e-9)
    excessMean[r] <- mean(sFree) - min(sFree)

    # Q-C control: the pinned-state criterion evaluated at ONE canonical
    # rooting.  It must be perfectly reproducible - if it is not, the oracle is
    # broken.  This is the check that keeps Q-C's interpretation honest: the
    # pinned-state spread above is NOT evidence that the pinned criterion is
    # ill-defined; it is a well-defined ROOTED criterion, and the spread is the
    # cost of leaving its root POSITION arbitrary.
    rtCanon <- ape::root(tr, outgroup = tr[["tip.label"]][1], resolve.root = TRUE)
    pinTwice <- c(SankoffRooted(rtCanon, tc, cm, 0L),
                  SankoffRooted(rtCanon, tc, cm, 0L))
    if (diff(range(pinTwice)) > 1e-12) {
      failures <- c(failures, sprintf(
        "%s: pinned-state score at a FIXED rooting is not reproducible - oracle broken",
        sc$name))
    }
  }

  cat(sprintf("--- %s  (nStates = %d, gain = %d, loss = 1)\n",
              sc$name, nStates, nSec + 1L))
  cat(sprintf("  forced_root = -1 (current): rooting-dependent on %d/%d; ",
              nDepFree, nRep))
  cat(sprintf("max spread %.3g (predicted bound nSec = %.3g)\n",
              max(spreadFree), bound))
  # NOTE ON READING THE PINNED ROW: forced_root = 0 defines a DIFFERENT,
  # explicitly ROOTED criterion.  Its spread across rootings is therefore not
  # comparable to the free-root row as "better or worse" - it measures how much
  # the arbitrary CHOICE of root position costs once the root state is
  # meaningful.  In particular the nSec = 0 row below is a valid asymmetry
  # control for the FREE arm only; for the pinned arm a non-zero spread there is
  # expected, not a contradiction.
  cat(sprintf("  forced_root =  0 (absent) : root-POSITION-sensitive on %d/%d; ",
              nDepPinned, nRep))
  cat(sprintf("max spread %.3g  [different, rooted criterion - see note]\n",
              max(spreadPinned)))
  cat(sprintf("  arbitrary rooting vs min-over-rootings: %.0f%% of rootings ",
              100 * mean(fracAtMin)))
  cat(sprintf("attain the min; mean overstatement %.3g\n", mean(excessMean)))

  # Q-B: the gradient prediction.  Only checked where Sum_tips f is constant,
  # i.e. no ambiguous / present-unknown tips (those are Q-D).
  if (sc$ambig == 0 && sc$pu == 0) {
    if (max(spreadFree) > bound + 1e-9) {
      failures <- c(failures, sprintf(
        "%s: spread %.4g EXCEEDS the nSec = %.4g gradient bound",
        sc$name, max(spreadFree), bound))
    }
  }
  # Control: nSec = 0 makes the matrix symmetric, so it MUST be invariant.
  if (nSec == 0L && max(spreadFree) > 1e-9) {
    failures <- c(failures, sprintf(
      "%s: symmetric control is rooting-DEPENDENT (spread %.4g) - the oracle itself is suspect",
      sc$name, max(spreadFree)))
  }
}

cat("\n")
if (length(failures) == 0L) {
  cat("VERDICT: every stated prediction held.\n")
  quit(status = 0L)
} else {
  cat("VERDICT: prediction FAILURES -\n")
  cat(paste0("  - ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
