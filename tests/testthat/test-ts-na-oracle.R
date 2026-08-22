# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Independent-reference oracle for the Brazeau, Guillerme and Smith (2019)
# three-pass inapplicable algorithm (src/ts_fitch_na.h).
#
# Every other NA test in the suite compares one TreeSearch code path against
# another (incremental vs full rescore, cached vs bounded candidate scorers).
# Such relative checks cannot see a defect shared by BOTH sides.  This file
# closes that gap: it reimplements the CRITERION from first principles in R,
# with no code path in common with src/ts_fitch*, and asserts that the engine
# reproduces it exactly.
#
# Criterion, as BGS define it operationally:
#   length = min, over applicability reconstructions that are FITCH-OPTIMAL
#            for the binary applicable/inapplicable character, of
#              sum(within-region Fitch length) + (tip-bearing regions - 1)
#
# The reference brute-forces every applicability labelling of the internal
# nodes, keeps those with the fewest applicability changes, and scores each by
# an explicit dynamic program over the resulting regions.  Tips are either "-"
# (inapplicable) or an applicable state, so their applicability is fixed: this
# is the regime in which the criterion is tie-break-free.  With partially
# ambiguous {-,X} tips BGS additionally prefers the APPLICABLE resolution on a
# tie (maximise homology; see test-ts-na-ambig.R), so the engine may
# legitimately exceed this minimum there -- hence the second test asserts only
# the one-sided bound that holds in every regime.

.NaRefPostorder <- function(kids, root) {
  out <- integer(0)
  stack <- root
  while (length(stack)) {
    nd <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    out <- c(nd, out)
    stack <- c(stack, kids[[nd]])
  }
  out
}

.NaRefLabelCost <- function(app, postorder, kids, parentOf, tipStates, nTip, k) {
  nNode <- length(app)
  cost <- matrix(Inf, nNode, k)
  hasTip <- logical(nNode)
  for (nd in postorder) {
    if (!app[[nd]]) next
    if (nd <= nTip) {
      cost[nd, tipStates[[nd]]] <- 0
      hasTip[[nd]] <- TRUE
    } else {
      accum <- rep(0, k)
      anyTip <- FALSE
      for (kd in kids[[nd]]) {
        if (!app[[kd]]) next
        child <- cost[kd, ]
        accum <- accum + vapply(seq_len(k),
                                function(s) min(child + (seq_len(k) != s)), 0)
        anyTip <- anyTip || hasTip[[kd]]
      }
      cost[nd, ] <- accum
      hasTip[[nd]] <- anyTip
    }
  }
  total <- 0
  nRegion <- 0L
  for (nd in seq_len(nNode)) {
    if (!app[[nd]]) next
    pa <- parentOf[[nd]]
    if (is.na(pa) || !app[[pa]]) {
      total <- total + min(cost[nd, ])
      if (hasTip[[nd]]) nRegion <- nRegion + 1L
    }
  }
  total + max(0L, nRegion - 1L)
}

.NaRefLength <- function(edge, nTip, tipStates, k) {
  nNode <- max(edge)
  parentOf <- rep(NA_integer_, nNode)
  parentOf[edge[, 2]] <- edge[, 1]
  kids <- lapply(seq_len(nNode), function(nd) edge[edge[, 1] == nd, 2])
  root <- setdiff(edge[, 1], edge[, 2])[[1]]
  postorder <- .NaRefPostorder(kids, root)
  tipApp <- vapply(tipStates, function(s) length(s) > 0L, TRUE)

  internals <- seq.int(nTip + 1L, nNode)
  grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(internals))))
  changes <- integer(nrow(grid))
  costs <- numeric(nrow(grid))
  for (r in seq_len(nrow(grid))) {
    app <- logical(nNode)
    app[seq_len(nTip)] <- tipApp
    app[internals] <- grid[r, ]
    changes[[r]] <- sum(app[edge[, 1]] != app[edge[, 2]])
    costs[[r]] <- .NaRefLabelCost(app, postorder, kids, parentOf,
                                  tipStates, nTip, k)
  }
  min(costs[changes == min(changes)])
}

test_that("NA three-pass matches an independent brute-force reference", {
  library("TreeTools", quietly = TRUE)
  lvls <- c("-", "1", "2", "3")
  set.seed(20260805)
  for (case in seq_len(40L)) {
    nTip <- 6L
    mat <- matrix(sample(lvls, nTip, TRUE, c(0.4, 0.22, 0.2, 0.18)),
                  nTip, 1L, dimnames = list(paste0("t", seq_len(nTip)), NULL))
    dataset <- phangorn::phyDat(mat, type = "USER", levels = lvls)
    at <- attributes(dataset)
    tree <- RenumberTips(RandomTree(rownames(mat), root = TRUE), names(dataset))
    edge <- Preorder(tree)[["edge"]]
    tipData <- matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE)
    engine <- TreeSearch:::ts_fitch_score(edge, at[["contrast"]], tipData,
                                          as.integer(at[["weight"]]),
                                          at[["levels"]], concavity = -1)
    tipStates <- lapply(seq_len(nTip), function(i) {
      allowed <- which(at[["contrast"]][tipData[i, 1], ] > 0)
      as.integer(allowed[allowed > 1L] - 1L)
    })
    reference <- .NaRefLength(edge, nTip, tipStates, length(lvls) - 1L)
    expect_equal(as.numeric(engine), as.numeric(reference),
                 info = paste0("case ", case, ": ",
                               paste(mat[, 1], collapse = " ")))
  }
})

test_that("BGS length never falls below the best possible reconstruction", {
  # One-sided invariant that holds however ambiguity ties are broken: the
  # reported length is the cost of SOME reconstruction, so it can never be
  # cheaper than the unconstrained minimum over all node labellings.
  library("TreeTools", quietly = TRUE)
  lvls <- c("-", "1", "2")
  ambig <- list("{-1}" = c("-", "1"), "{-2}" = c("-", "2"),
                "{12}" = c("1", "2"))
  allTok <- c(lvls, names(ambig))
  contrast <- matrix(0, length(allTok), length(lvls),
                     dimnames = list(allTok, lvls))
  for (s in lvls) contrast[s, s] <- 1
  for (nm in names(ambig)) contrast[nm, ambig[[nm]]] <- 1

  set.seed(1234L)
  for (case in seq_len(25L)) {
    nTip <- 6L
    mat <- matrix(sample(allTok, nTip, TRUE), nTip, 1L,
                  dimnames = list(paste0("t", seq_len(nTip)), NULL))
    dataset <- phangorn::phyDat(mat, type = "USER", levels = lvls,
                                ambiguity = names(ambig), contrast = contrast)
    at <- attributes(dataset)
    tree <- RenumberTips(RandomTree(rownames(mat), root = TRUE), names(dataset))
    edge <- Preorder(tree)[["edge"]]
    tipData <- matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE)
    engine <- TreeSearch:::ts_fitch_score(edge, at[["contrast"]], tipData,
                                          as.integer(at[["weight"]]),
                                          at[["levels"]], concavity = -1)
    nNode <- max(edge)
    opts <- lapply(seq_len(nNode), function(nd) {
      if (nd > nTip) return(0:2)
      as.integer(which(at[["contrast"]][tipData[nd, 1], ] > 0) - 1L)
    })
    g <- as.matrix(do.call(expand.grid, opts))
    pa <- g[, edge[, 1], drop = FALSE]
    ch <- g[, edge[, 2], drop = FALSE]
    bothApp <- (pa > 0) & (ch > 0)
    lower <- min(rowSums(bothApp & (pa != ch)) +
                   pmax(0, rowSums(g > 0) - rowSums(bothApp) - 1))
    expect_gte(as.numeric(engine), lower)
  }
})
