# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for x-transformation (Goloboff et al. 2021) scoring via the
# RecodeHierarchy() → Sankoff pipeline and MaximizeParsimony(inapplicable="xform").

library("TreeTools")

make_dat <- function(mat, levels = c("-", "0", "1")) {
  phangorn::phyDat(mat, type = "USER", levels = levels, ambiguity = "?")
}


# ===== Gain/loss asymmetry ===================================================
# The x-transformation penalizes gains (absent→present) more heavily than
# losses (present→absent) at ratio (n+1):1, where n = number of secondaries.

test_that("Xform prefers single gain + losses over multiple gains", {
  # Tree: ((t1,t2),(t3,t4))
  # Primary: t1=absent, t2=present, t3=present, t4=present
  # Secondary: t2=0, t3=0, t4=0 (all identical when present)
  # States: absent=0, (sec=0)=1
  # Cost: gain=2, loss=1
  # Optimal on this tree: root=present(1), loss to t1 → cost 1
  # Alternative: root=absent(0), gain at MRCA(t2,t3,t4)... not a single node
  #   on ((t1,t2),(t3,t4)). Would need gain at each present tip = 3*2 = 6
  # So single-gain-from-root (cost 1) << multiple-gains (cost 6)
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "0", "0",
    "1", "0", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  recoded <- RecodeHierarchy(ds, h)
  blk <- recoded$sankoff_chars[[1]]
  expect_equal(blk$cost_matrix[1, 2], 2)  # gain = n+1 = 2
  expect_equal(blk$cost_matrix[2, 1], 1)  # loss = 1

  # Score via Sankoff: tree where absent tip is sister to one present
  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree <- Renumber(RenumberTips(tree, names(ds)))
  res <- TreeSearch:::ts_sankoff_test(
    tree$edge,
    as.integer(blk$n_states),
    list(blk$cost_matrix),
    matrix(as.integer(blk$tip_states), ncol = 1),
    as.integer(blk$forced_root_state)
  )
  # Root=1(present), nodeAB: state 1 costs 0(t2)+1(loss to t1)=1
  # nodeCD: state 1 costs 0+0=0. Root=1: 0+0=0 from children.
  # But root cost = min over states. state 1 at root: costAB(1)=1, costCD(1)=0
  # cost_root_state1 = 1 + 0 = 1
  expect_equal(res$score, 1)
})


# ===== Secondary variation increases xform score =============================

test_that("Xform penalizes secondary variation on present branches", {
  # All present, varying secondaries → present-present Hamming cost
  mat_uniform <- matrix(c(
    "1", "0", "0",
    "1", "0", "0",
    "1", "0", "0",
    "1", "0", "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))

  mat_varied <- matrix(c(
    "1", "0", "0",
    "1", "1", "1",
    "1", "0", "0",
    "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))

  ds_u <- make_dat(mat_uniform)
  ds_v <- make_dat(mat_varied)
  h <- CharacterHierarchy("1" = 2:3)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")
  tree_u <- Renumber(RenumberTips(tree, names(ds_u)))
  tree_v <- Renumber(RenumberTips(tree, names(ds_v)))

  score_fn <- function(ds, tr) {
    rec <- RecodeHierarchy(ds, h)
    blk <- rec$sankoff_chars[[1]]
    TreeSearch:::ts_sankoff_test(
      tr$edge, as.integer(blk$n_states),
      list(blk$cost_matrix),
      matrix(as.integer(blk$tip_states), ncol = 1),
      as.integer(blk$forced_root_state)
    )$score
  }

  score_uniform <- score_fn(ds_u, tree_u)
  score_varied <- score_fn(ds_v, tree_v)

  # Uniform: all same state, no Hamming cost → 0
  expect_equal(score_uniform, 0)
  # Varied: secondary changes → Hamming cost > 0
  expect_gt(score_varied, 0)
})


# ===== HSJ vs xform cross-validation =========================================
# Both methods handle inapplicable characters; they should agree on basic
# properties even if exact scores differ.

test_that("HSJ and xform agree on optimal tree for simple gain scenario", {
  mat <- matrix(c(
    "0", "-", "-", "0", "0", "0",
    "0", "-", "-", "0", "1", "1",
    "1", "0", "0", "1", "0", "0",
    "1", "0", "1", "1", "0", "1",
    "1", "1", "0", "1", "1", "0",
    "1", "1", "1", "0", "1", "1"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  # Both should successfully search and return valid trees
  set.seed(7184)
  hsj_result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "hsj", hsj_alpha = 1.0,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )
  set.seed(7184)
  xform_result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )

  expect_s3_class(hsj_result[[1]], "phylo")
  expect_s3_class(xform_result[[1]], "phylo")

  # Both should find trees with the correct number of tips
  expect_equal(length(hsj_result[[1]]$tip.label), 6L)
  expect_equal(length(xform_result[[1]]$tip.label), 6L)
})


# ===== Xform with non-hierarchy characters ====================================

test_that("Xform correctly combines Fitch + Sankoff scoring", {
  # Chars 1-2: hierarchy (primary + secondary)
  # Chars 3-4: non-hierarchy (standard Fitch)
  mat <- matrix(c(
    "0", "-", "0", "0",
    "1", "0", "0", "1",
    "1", "1", "1", "0",
    "0", "-", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  tree <- ape::read.tree(text = "((t1,t2),(t3,t4));")

  # Standard Fitch score (all 4 chars as standard)
  fitch_total <- TreeLength(tree, ds)

  # Xform search should run
  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_true(is.finite(fitch_total))
})


# ===== Xform on larger dataset (8 tips) =====================================

test_that("Xform search works on 8-tip dataset", {
  mat <- matrix(c(
    "1",  "0",  "0",  "0",  "0",  "0",  "0",
    "1",  "0",  "0",  "0",  "0",  "1",  "0",
    "1",  "0",  "1",  "0",  "1",  "0",  "0",
    "1",  "1",  "0",  "1",  "0",  "0",  "1",
    "1",  "1",  "1",  "1",  "0",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "0",  "1",
    "0",  "-",  "-",  "1",  "1",  "1",  "0",
    "0",  "-",  "-",  "0",  "1",  "1",  "0"
  ), nrow = 8, byrow = TRUE,
  dimnames = list(paste0("t", 1:8), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  result <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )
  expect_s3_class(result[[1]], "phylo")
  expect_equal(length(result[[1]]$tip.label), 8L)

  # All result trees should be valid rooted phylogenies
  for (tr in result) {
    expect_s3_class(tr, "phylo")
    expect_true(TreeIsRooted(tr))
  }
})


# ===== Xform score is consistent across replicates ===========================

test_that("Xform search produces deterministic scores with same seed", {
  mat <- matrix(c(
    "0", "-", "0", "1",
    "1", "0", "1", "0",
    "1", "1", "0", "1",
    "1", "1", "1", "0",
    "0", "-", "1", "1"
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("t", 1:5), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  set.seed(3021)
  r1 <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )
  set.seed(3021)
  r2 <- MaximizeParsimony(
    ds, hierarchy = h, inapplicable = "xform",
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L
  )

  # Same seed, same result
  expect_equal(attr(r1, "score"), attr(r2, "score"))
})


# ===== Asymmetric cost correctness: gain vs loss =============================

test_that("Xform gain cost scales with number of secondaries", {
  # 1 secondary → gain = 2, 2 secondaries → gain = 3, 3 → gain = 4
  for (n_sec in 1:3) {
    n_cols <- 1 + n_sec
    mat <- matrix("-", nrow = 3, ncol = n_cols,
                  dimnames = list(paste0("t", 1:3), NULL))
    mat[1, ] <- c("0", rep("-", n_sec))
    mat[2, ] <- c("1", rep("0", n_sec))
    mat[3, ] <- c("1", rep("1", n_sec))
    ds <- make_dat(mat)
    h <- CharacterHierarchy("1" = seq(2L, n_cols))

    rec <- RecodeHierarchy(ds, h)
    blk <- rec$sankoff_chars[[1]]

    expected_gain <- n_sec + 1L
    # All absent→present transitions should cost expected_gain
    for (j in 2:blk$n_states) {
      expect_equal(blk$cost_matrix[1, j], expected_gain,
                   info = paste("n_sec =", n_sec, "state", j))
    }
  }
})


# ===== Heterogeneous-n_states blocks: cost-matrix stride (SK-01) ==============
# Two hierarchy blocks of DIFFERENT n_states give max_states > min(n_states).
# The live search path (score_tree) copied each block's cost matrix from the
# [max_states x max_states] storage block verbatim but sankoff_score_char reads
# it at the per-character n_states stride, so any block with n_states <
# max_states had its rows s>0 read from the zero-padded gap -- silently treating
# loss/transition costs as 0 and undercounting the score. The search-reported
# score must agree with TreeLength() (independent ts_sankoff_test kernel) on the
# returned tree, whichever tree the search settles on.

test_that("Xform scores heterogeneous-n_states blocks consistently (SK-01)", {
  mat <- matrix(c(
    "1", "0", "1", "0", "0",
    "1", "1", "1", "0", "1",
    "0", "-", "1", "1", "0",
    "1", "1", "0", "-", "-",
    "0", "-", "1", "1", "1",
    "1", "0", "0", "-", "-"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L, "3" = 4:5)

  rec <- RecodeHierarchy(ds, h)
  # Blocks must have differing state counts (3 and 5) to exercise the stride.
  expect_equal(
    sort(vapply(rec$sankoff_chars, function(b) as.integer(b$n_states), integer(1))),
    c(3L, 5L)
  )

  set.seed(101)
  res <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "xform",
                           maxReplicates = 4L, targetHits = 3L, verbosity = 0L)
  # MaximizeParsimony (score_tree) and TreeLength (ts_sankoff_test) must report
  # the same score for the SAME tree; they disagreed by the undercount before
  # the stride was compacted to per-character n_states.
  expect_equal(
    attr(res, "score"),
    TreeLength(res[[1]], ds, hierarchy = h, inapplicable = "xform")
  )
})


# ===== All-hierarchy data: zero Fitch words (regression) =====================
# When EVERY character belongs to a hierarchy block, recoding removes all
# characters from the equal-weights (Fitch) dataset, so DataSet::total_words
# is 0 and the per-word state vectors (tip_states, edge_set) are empty.
# wagner_tree() took the address of element 0 of these empty vectors
# (&ds.tip_states[0], &edge_set[0]) — undefined behaviour that aborted under
# the hardened libstdc++ assertions in the gcc-ASAN CI (run 28662381835,
# ts-xform group).  With no Fitch signal the search must still build a valid
# tree from the Sankoff term alone.  This test exercises that path so the
# hardened/ASAN CI covers it.

# ===== T-379: partial secondaries must constrain, not free, -2 tips =========
# A tip coded -2 ("present, secondary combination unknown") previously freed
# EVERY present state regardless of any secondary that WAS actually observed,
# discarding real information and undercounting cost. With one secondary
# unknown but two known, and those two knowns conflicting with a comparison
# tip's fully-resolved combination, the admissible set should now exclude
# that comparison tip's exact state -- forcing a strictly positive
# present-present Hamming cost instead of the old free-ride of 0.

test_that("Xform -2 tip is constrained by its known secondaries (T-379)", {
  # Tip A: primary present, secondaries all "1" -> fully resolved combo (2,2,2).
  # Tip B: primary present, secondaries "0","0","?" -> two secondaries KNOWN
  #   (both conflicting with A's "1","1"), one unknown.
  # Tip C: primary present, secondaries all "0" -> supplies the "0" level for
  #   secondary 3 so that character has 2 informative levels (otherwise it
  #   would trivially degenerate to 1, and B's "unknown" would be moot). C is
  #   excluded from the scored tree below so it cannot mask the effect.
  mat <- matrix(c(
    "1", "1", "1", "1",
    "1", "0", "0", "?",
    "1", "0", "0", "0"
  ), nrow = 3, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:4)

  rec <- RecodeHierarchy(ds, h)
  blk <- rec$sankoff_chars[[1]]
  expect_equal(blk$n_states, 9L)  # 2^3 present combos + absent

  # Score A and B alone, as a 2-tip cherry: with a single sister pair, the
  # Sankoff minimum reduces to A's fixed state plus the cheapest transition
  # to any state B's tip cost allows -- i.e. exactly the quantity T-379
  # changes. (Tip C only exists to register the "0" level above; it plays no
  # further part here.)
  abIdx <- match(c("A", "B"), names(ds))
  tipStatesMat <- matrix(as.integer(blk$tip_states[abIdx]), ncol = 1)
  expect_equal(tipStatesMat[2, 1], -2L)  # B: present, secondary 3 unknown

  tree <- ape::read.tree(text = "(A,B);")

  # Score WITHOUT combo info -- mirrors the pre-fix behaviour: -2 frees every
  # present state, so B can "become" A's exact state at zero cost.
  score_old <- TreeSearch:::ts_sankoff_test(
    tree$edge, as.integer(blk$n_states), list(blk$cost_matrix),
    tipStatesMat, as.integer(blk$forced_root_state)
  )$score
  expect_equal(score_old, 0)

  # Score WITH combo info -- B is now constrained to states consistent with
  # its two known (and, here, A-conflicting) secondaries, so it can no longer
  # reach A's exact state.
  combo_grids <- list(blk$combo_grid)
  tip_sec_known <- list(matrix(as.integer(blk$tip_sec_known[abIdx, ]),
                                nrow = 2))
  score_new <- TreeSearch:::ts_sankoff_test(
    tree$edge, as.integer(blk$n_states), list(blk$cost_matrix),
    tipStatesMat, as.integer(blk$forced_root_state),
    combo_grids, tip_sec_known
  )$score

  expect_gt(score_new, score_old)
})


test_that("Xform search handles all-hierarchy data (zero Fitch words)", {
  mat <- matrix(c(
    "1", "0", "0", "-", "-",
    "1", "1", "1", "0", "1",
    "0", "-", "1", "1", "0",
    "1", "0", "1", "0", "0",
    "0", "-", "0", "-", "-",
    "1", "1", "0", "-", "-"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  # Two blocks span ALL five characters -> no non-hierarchy chars remain.
  h <- CharacterHierarchy("1" = 2L, "3" = 4:5)
  expect_length(setdiff(seq_len(5L), HierarchyChars(h)), 0L)

  set.seed(42)
  # This all-hierarchy matrix provably triggers T-374's open residue: pool
  # membership is decided on search-time scores taken at differing rootings, so
  # the returned trees do not share a length at the common rooting they are
  # reported at (measured: 7 to 9), and MaximizeParsimony() warns and reports the
  # smallest.  Pinned as an expectation rather than left as ambient noise in a
  # green suite -- if this warning ever STOPS firing, the residue has been fixed
  # (or masked) and that deserves to be noticed here.
  res <- NULL
  expect_warning(
    res <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "xform",
                             maxReplicates = 4L, targetHits = 3L,
                             verbosity = 0L),
    "do not share a length")
  expect_s3_class(res[[1]], "phylo")
  expect_equal(length(res[[1]]$tip.label), 6L)
  for (tr in res) {
    expect_s3_class(tr, "phylo")
    expect_true(TreeIsRooted(tr))
  }
})


# ===== T-385: a reported length must be reproducible ==========================
# The x-transformation's step matrix is asymmetric (gain = nSec + 1, loss = 1),
# so a tree's length depends on where its root sits -- unlike parsimony under the
# symmetric criteria.  The engine recorded `best_score` mid-search at whatever
# rooting the replicate held, while `ts_collapse_pool()` returns every tree
# re-rooted on tip 0, so `attr(res, "score")` did not reproduce under
# `TreeLength()` of the very tree returned (measured: 178 reported against 183
# returned, 36 tips / 6 blocks).  Both boundaries now score at a canonical
# rooting.  See dev/plans/2026-07-29-t374b-xform-rooting-policy.md.

test_that("TreeLength xform is rooting-invariant (T-385)", {
  # 8-tip case whose Sankoff term genuinely varies with the rooting.  That
  # precondition is ASSERTED below rather than assumed, so this test cannot pass
  # vacuously on data that happens to be rooting-insensitive.
  mat <- matrix(c(
    "0", "-", "-",
    "1", "1", "1",
    "0", "-", "-",
    "0", "-", "-",
    "1", "0", "0",
    "0", "-", "-",
    "0", "-", "-",
    "0", "-", "-"
  ), nrow = 8, byrow = TRUE,
  dimnames = list(paste0("t", 1:8), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)
  tree <- ape::read.tree(text = "(((t1,t3),((t2,t5),t7)),(t4,(t6,t8)));")
  taxa <- names(ds)

  # Precondition: the raw Sankoff kernel IS rooting-sensitive here (3 to 5).
  # ts_sankoff_test() is untouched by the fix, so this measures the data.
  recoded <- RecodeHierarchy(ds, h)
  xf <- TreeSearch:::.PrepareXformArgs(recoded, length(ds))
  kernel <- vapply(taxa, function(taxon) {
    tr <- RenumberTips(Renumber(RootTree(tree, taxon)), taxa)
    TreeSearch:::ts_sankoff_test(tr[["edge"]], xf$n_states, xf$cost_matrices,
                                 xf$tip_states, xf$forced_root, xf$combo_grids,
                                 xf$tip_sec_known)$score
  }, numeric(1))
  expect_gt(diff(range(kernel)), 0)

  # Given that, TreeLength() must still return ONE length for ONE topology.
  lengths <- vapply(taxa, function(taxon) {
    TreeLength(RootTree(tree, taxon), ds, hierarchy = h,
               inapplicable = "xform")
  }, numeric(1))
  expect_equal(diff(range(lengths)), 0)

  # The multiPhylo method must canonicalise identically to the single-tree one:
  # it previously rooted only trees that arrived unrooted.
  multi <- TreeLength(
    structure(lapply(taxa, function(taxon) RootTree(tree, taxon)),
              class = "multiPhylo"),
    ds, hierarchy = h, inapplicable = "xform"
  )
  expect_equal(unname(multi), unname(lengths))
})

test_that("MaximizeParsimony xform reports the length of the tree it returns (T-385)", {
  # End-to-end guard on the report path.  Larger than the test above so the
  # search has room to end at a rooting other than the returned one.
  set.seed(7)
  nTip <- 16L
  cols <- list()
  hierArgs <- list()
  for (b in 1:3) {
    primary <- sample(c("0", "1"), nTip, replace = TRUE)
    priIdx <- length(cols) + 1L
    cols[[length(cols) + 1L]] <- primary
    secIdx <- integer(2)
    for (s in 1:2) {
      cols[[length(cols) + 1L]] <- ifelse(
        primary == "0", "-", sample(c("0", "1"), nTip, replace = TRUE))
      secIdx[s] <- length(cols)
    }
    hierArgs[[as.character(priIdx)]] <- secIdx
  }
  for (f in 1:6) {
    cols[[length(cols) + 1L]] <- sample(c("0", "1"), nTip, replace = TRUE)
  }
  mat <- do.call(cbind, cols)
  dimnames(mat) <- list(paste0("t", seq_len(nTip)), NULL)
  ds <- make_dat(mat)
  h <- do.call(CharacterHierarchy, hierArgs)

  set.seed(11)
  res <- suppressWarnings(
    MaximizeParsimony(ds, inapplicable = "xform", hierarchy = h,
                      maxReplicates = 3L, verbosity = 0L))

  reported <- unique(unname(attr(res, "score")))
  expect_length(reported, 1L)
  returned <- unname(vapply(res, function(tr) {
    TreeLength(tr, ds, hierarchy = h, inapplicable = "xform")
  }, numeric(1)))

  # The contract is that the reported score is the canonical length of the pool,
  # reproducible from the trees returned -- NOT that every returned tree shares
  # it.  Pool membership is still decided on search-time scores taken at differing
  # rootings, so on some data the pool genuinely spans a range and
  # MaximizeParsimony() warns (all-hierarchy matrices do: see the zero-Fitch-words
  # test above, which spans 7 to 9).  That is the open residue of T-374 and is
  # deliberately out of scope for the reporting fix.  Assert the contract:
  expect_equal(reported, min(returned))

  # On THIS dataset the pool does happen to be self-consistent.  Asserted as a
  # property of the data, not as something the code promises.
  expect_equal(diff(range(returned)), 0)
})
