library("TreeTools", quietly = TRUE)

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

# --- Input validation ---

test_that("MaximizeParsimony stops with message when dataset is NULL", {
  expect_error(
    MaximizeParsimony(NULL, maxReplicates = 1L, targetHits = 1L,
                      verbosity = 0L),
    "`dataset` cannot be NULL."
  )
})

# --- Strategy presets ---

test_that("strategy = 'sprint' runs and returns valid result", {
  set.seed(3418)
  result <- MaximizeParsimony(ds, strategy = "sprint",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_true(attr(result, "score") > 0)
  expect_equal(NTip(result[[1]]), NTip(ds))
})

test_that("candidates_evaluated attribute is reported for serial search", {
  # Diagnostic counter (TNT "rearrangements examined" analogue): a positive,
  # finite scalar for a single-threaded search. See MaximizeParsimony @return.
  set.seed(3418)
  result <- MaximizeParsimony(ds, strategy = "sprint",
                               maxReplicates = 2L, targetHits = 1L,
                               nThreads = 1L, verbosity = 0L)
  ce <- attr(result, "candidates_evaluated")
  expect_type(ce, "double")
  expect_length(ce, 1L)
  expect_true(is.finite(ce) && ce > 0)
})

test_that("strategy = 'intensive' (opt-in) runs and returns valid result", {
  set.seed(5726)
  result <- MaximizeParsimony(ds, strategy = "intensive",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_equal(attr(result, "score"), TreeLength(result[[1]], ds),
               tolerance = 0.01)
})

test_that("strategy = 'default' runs and returns valid result", {
  set.seed(5726)
  result <- MaximizeParsimony(ds, strategy = "default",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("strategy = 'thorough' runs and returns valid result", {
  set.seed(8103)
  result <- MaximizeParsimony(ds, strategy = "thorough",
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that(".AutoStrategy selects on size and signal density", {
  AS <- TreeSearch:::.AutoStrategy

  # Small datasets: always sprint
  expect_equal(AS(20, 100), "sprint")
  expect_equal(AS(30, 500), "sprint")

  # Few chars (< 100 patterns) -> flat landscape -> always default
  expect_equal(AS(50, 25),   "default")  # small, very few chars
  expect_equal(AS(65, 80),   "default")  # large enough tip count, but nChar < 100
  expect_equal(AS(200, 99),  "default")  # large, but still nChar < 100

  # Mid-size (31-64 tips) with enough chars -> default (not large enough)
  expect_equal(AS(60, 300), "default")  # nChar >= 100 but nTip < 65
  expect_equal(AS(64, 200), "default")  # nChar >= 100 but nTip < 65

  # Large (>= 65 tips) with enough chars -> thorough
  # Signal density does NOT gate thorough: more chars = more benefit (T-068 benchmark)
  expect_equal(AS(65, 100),   "thorough")  # boundary case: 65 tips, 100 chars
  expect_equal(AS(74, 200),   "thorough")  # 74 tips, ratio 2.7
  expect_equal(AS(75, 250),   "thorough")  # ratio 3.3
  expect_equal(AS(100, 200),  "thorough")  # ratio 2.0
  expect_equal(AS(75, 400),   "thorough")  # ratio 5.3 — high ratio still benefits
  expect_equal(AS(119, 2800), "thorough")  # just below large threshold
  expect_equal(AS(125, 2800), "large")    # >= 120 tips -> large
  expect_equal(AS(200, 100),  "large")    # >= 120 tips -> large
  expect_equal(AS(200, 1200), "large")    # >= 120 tips -> large
})

test_that("strategy = 'auto' selects based on dataset size", {
  set.seed(2944)
  # Vinther2008 has 23 tips -> should auto-select "sprint"
  result <- MaximizeParsimony(ds, strategy = "auto",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("strategy = 'none' uses raw parameter defaults", {
  set.seed(6017)
  result <- MaximizeParsimony(ds, strategy = "none",
                               maxReplicates = 2L, targetHits = 1L,
                               ratchetCycles = 1L, driftCycles = 0L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("explicit params override strategy preset", {
  set.seed(1589)
  # Sprint has driftCycles=0; override to 1
  result <- MaximizeParsimony(ds, strategy = "sprint",
                               driftCycles = 1L,
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("unknown strategy gives warning", {
  set.seed(4821)
  expect_warning(
    MaximizeParsimony(ds, strategy = "nonexistent",
                      maxReplicates = 1L, targetHits = 1L,
                      verbosity = 0L),
    "Unknown strategy"
  )
})

# --- maxSeconds timeout ---

test_that("maxSeconds stops search before maxReplicates", {
  set.seed(7392)
  # Use a near-zero timeout so even a small dataset triggers it reliably.
  # The timeout is checked between replicates, so the first may complete

  # but subsequent ones won't start.
  result <- MaximizeParsimony(ds, maxReplicates = 1000L, targetHits = 1000L,
                               maxSeconds = 0.001, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_lt(attr(result, "replicates"), 1000L)
})

test_that("maxSeconds = 0 means no timeout", {
  set.seed(8456)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               maxSeconds = 0, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_false(attr(result, "timed_out"))
})

test_that("perturbStopFactor fires and sets perturb_stop attribute", {
  # perturbStopFactor=1 on Vinther2008 (23 tips) means limit = 23 reps.
  set.seed(4618)
  result <- MaximizeParsimony(ds, maxReplicates = 500L, targetHits = 500L,
                               control = SearchControl(
                                 perturbStopFactor = 1L,
                                 ratchetCycles = 1L),
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_lt(attr(result, "replicates"), 500L)
  expect_true(attr(result, "perturb_stop"))
  expect_false(attr(result, "timed_out"))
})

test_that("verbosity = 1 prints 'Search complete' summary to console", {
  set.seed(3071)
  # MaximizeParsimony emits two streams at verbosity = 1: cli messages
  # via message() ("Strategy: ...", "Search complete: ...") and C++
  # Rprintf progress via stdout ("Replicate N/M", "Converged: ...").
  # Capture both so they don't leak into testthat output, then assert
  # that the expected lines were produced.
  msg_lines <- character()
  stdout_lines <- capture.output(
    msg_lines <- capture.output(
      MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                        verbosity = 1L),
      type = "message"
    )
  )
  expect_true(any(grepl("Search complete", msg_lines)))
  expect_true(any(grepl("Replicate", stdout_lines)))
  expect_true(any(grepl("Converged|score", stdout_lines)))
})

# --- nThreads ---

test_that("nThreads = 1 (serial) runs correctly", {
  set.seed(5193)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               nThreads = 1L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("nThreads = 2 (parallel) runs correctly", {
  set.seed(6274)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               nThreads = 2L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  # Score should be reasonable (not garbage from parallel corruption)
  expect_true(attr(result, "score") < 200)
})

# --- User-supplied starting tree (warm-start) ---

test_that("user tree is used as warm start", {
  # Build a known tree
  set.seed(9847)
  user_tree <- RandomTree(ds, root = TRUE)
  user_tree <- Preorder(user_tree)

  result <- MaximizeParsimony(ds, tree = user_tree,
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  result_score <- attr(result, "score")

  # Result should be at least as good as the input tree
  input_score <- TreeLength(user_tree, ds)
  expect_true(result_score <= input_score)
})

test_that("multiPhylo input uses first tree as warm start", {
  set.seed(3571)
  trees <- list(
    RandomTree(ds, root = TRUE),
    RandomTree(ds, root = TRUE)
  )
  class(trees) <- "multiPhylo"
  trees <- Preorder(trees)

  result <- MaximizeParsimony(ds, tree = trees,
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

# --- timings attribute ---

test_that("timings attribute is returned", {
  set.seed(2689)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  timings <- attr(result, "timings")
  expect_false(is.null(timings))
  expect_true(is.numeric(timings))
  expect_true(all(timings >= 0))
  expect_true("wagner_ms" %in% names(timings))
  expect_true("ratchet_ms" %in% names(timings))
})

# --- IW with strategy ---

test_that("IW mode works with strategy presets", {
  set.seed(4012)
  result <- MaximizeParsimony(ds, concavity = 10, strategy = "sprint",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  cpp_score <- attr(result, "score")
  tl_score <- TreeLength(result[[1]], ds, concavity = 10)
  expect_equal(cpp_score, tl_score, tolerance = 0.01)
})

# --- Output tree validity ---

test_that("output trees have valid preorder numbering", {
  set.seed(8734)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  for (tree in result) {
    # Trees should have correct number of tips and edges
    expect_equal(NTip(tree), NTip(ds))
    expect_equal(nrow(tree$edge), 2 * (NTip(ds) - 1))
    # Tip labels should match
    expect_true(all(TipLabels(tree) %in% names(ds)))
  }
})

# --- T-039 regression: constraint on small fully-resolved trees ---

test_that("Fully-resolving constraint on 5-tip tree does not crash", {
  ds5 <- phangorn::phyDat(
    matrix(c("0","0","0","1","1","0","1","0","1","0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(4172)
  result <- MaximizeParsimony(ds5, constraint = cons,
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_equal(NTip(result[[1]]), 5L)
})

test_that("AdditionTree with fully-resolving constraint works", {
  ds5 <- phangorn::phyDat(
    matrix(c("0","0","0","1","1","0","1","0","1","0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(6091)
  wt <- AdditionTree(ds5, constraint = cons)
  expect_s3_class(wt, "phylo")
  expect_equal(NTip(wt), 5L)
  expect_equal(nrow(wt$edge), 8L)
})

test_that("Constrained Wagner tree works with multiple seeds", {
  ds5 <- phangorn::phyDat(
    matrix(c("0","0","0","1","1","0","1","0","1","0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")
  at <- attributes(ds5)
  consArgs <- TreeSearch:::.PrepareConstraint(cons, ds5)

  for (s in c(1, 6, 42)) {
    set.seed(s)
    result <- do.call(TreeSearch:::ts_random_wagner_tree, c(
      list(contrast = at$contrast,
           tip_data = matrix(unlist(ds5, use.names = FALSE), nrow = 5, byrow = TRUE),
           weight = at$weight, levels = at$levels),
      consArgs))
    expect_true(is.finite(result$score), info = paste("seed", s))
    expect_equal(nrow(result$edge), 8L, info = paste("seed", s))
  }
})

# --- Intra-replicate fusing (T-258) ---

test_that("intraFuse runs without error", {
  set.seed(8517)
  result <- MaximizeParsimony(ds, strategy = "sprint",
                              maxReplicates = 5L, targetHits = 2L,
                              maxSeconds = 3, intraFuse = TRUE,
                              verbosity = 0L, nThreads = 1L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_lte(attr(result, "score"), 100)  # should find reasonable score
})

test_that("intraFuse with dataset size change does not crash", {
  ds_large <- inapplicable.phyData[["Agnarsson2004"]]  # 62 tips
  ds_small <- inapplicable.phyData[["Vinther2008"]]     # 23 tips

  # Run on larger dataset first with intra-fuse
  set.seed(9014)
  r1 <- MaximizeParsimony(ds_large, strategy = "sprint",
                          maxReplicates = 3L, targetHits = 2L,
                          maxSeconds = 3, intraFuse = TRUE,
                          verbosity = 0L, nThreads = 1L)
  expect_true(is.finite(attr(r1, "score")))

  # Then run on smaller dataset with intra-fuse (regression test for segfault)
  set.seed(9015)
  r2 <- MaximizeParsimony(ds_small, strategy = "sprint",
                          maxReplicates = 3L, targetHits = 2L,
                          maxSeconds = 3, intraFuse = TRUE,
                          verbosity = 0L, nThreads = 1L)
  expect_true(is.finite(attr(r2, "score")))
})

# --- collapse = TRUE (TNT-style zero-length-branch contraction) ---

# A soft polytomy: (b1,b2) & (b3,b4) supported cherries plus a fan of identical
# taxa that has no internal support.  Every resolution of the fan is an MPT, so
# collapse = FALSE returns many fully-resolved trees while collapse = TRUE must
# contract the fan to a single polytomy -> exactly one distinct collapsed tree.
.SoftPolytomyData <- function(nFan = 5L) {
  backbone <- c("b1", "b2", "b3", "b4")
  fan <- paste0("f", seq_len(nFan))
  taxa <- c(backbone, fan)
  mat <- matrix("0", nrow = length(taxa), ncol = 0,
                dimnames = list(taxa, NULL))
  addChar <- function(m, ones) {
    cbind(m, ifelse(rownames(m) %in% ones, "1", "0"))
  }
  for (i in 1:3) mat <- addChar(mat, fan)            # fan synapomorphy x3
  for (i in 1:2) mat <- addChar(mat, c("b1", "b2"))  # backbone cherry x2
  for (i in 1:2) mat <- addChar(mat, c("b3", "b4"))  # backbone cherry x2
  phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))
}

test_that("collapse = TRUE contracts a soft polytomy to one collapsed tree", {
  phy <- .SoftPolytomyData(5L)
  set.seed(1L)
  resolved <- MaximizeParsimony(
    phy, concavity = Inf, maxReplicates = 12L, strategy = "intensive",
    control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
    verbosity = 0L, collapse = FALSE)
  set.seed(1L)
  collapsed <- MaximizeParsimony(
    phy, concavity = Inf, maxReplicates = 12L, strategy = "intensive",
    control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
    verbosity = 0L, collapse = TRUE)

  # Many resolved variants of one collapsed topology.
  expect_gt(length(resolved), 1L)
  # All collapse to a single distinct topology.
  expect_equal(length(collapsed), 1L)
  expect_equal(attr(collapsed, "n_topologies"), 1L)
  # The single tree is a genuine polytomy (fewer internal nodes than binary).
  expect_lt(collapsed[[1]][["Nnode"]], NTip(phy) - 1L)
  # The supported backbone cherries survive collapse (the fan, not the
  # backbone, is contracted).  Root on a fan tip so neither cherry sits at the
  # root, then each cherry is the exclusive descendant set of its MRCA.
  rooted <- TreeTools::RootTree(collapsed[[1]], "f1")
  isClade <- function(members) {
    mrca <- ape::getMRCA(rooted, members)
    !is.null(mrca) &&
      setequal(rooted[["tip.label"]][phangorn::Descendants(rooted, mrca,
                                                           "tips")[[1]]],
               members)
  }
  expect_true(isClade(c("b1", "b2")))
  expect_true(isClade(c("b3", "b4")))
})

test_that("collapse = TRUE is a no-op when no branch is unsupported", {
  # Vinther2008 MPTs are fully resolved (no zero-length branches): collapse must
  # leave every tree binary and topologically unchanged.
  set.seed(3418)
  resolved <- MaximizeParsimony(ds, strategy = "sprint",
                                maxReplicates = 3L, targetHits = 1L,
                                verbosity = 0L, collapse = FALSE)
  set.seed(3418)
  collapsed <- MaximizeParsimony(ds, strategy = "sprint",
                                 maxReplicates = 3L, targetHits = 1L,
                                 verbosity = 0L, collapse = TRUE)
  nTip <- NTip(ds)
  # No spurious collapse: every returned tree stays fully binary.
  expect_true(all(vapply(collapsed, function(t) t[["Nnode"]] == nTip - 1L,
                         logical(1))))
  # Same set of topologies (each collapsed tree matches a resolved tree, RF 0).
  expect_true(all(apply(
    as.matrix(TreeDist::RobinsonFoulds(collapsed, resolved)), 1L, min) == 0))
})

test_that("collapse = TRUE shows an enforced clade and collapses the rest", {
  # "Show the enforced clade": a constraint encodes external evidence for a
  # grouping the matrix does not support, so it must stay visible even though it
  # sits on a zero-length branch.  Forcing (f1, f2) inside the otherwise
  # unsupported fan gives it a zero-length branch; collapse must keep it while
  # still contracting the unconstrained fan resolutions (the philosophy that
  # makes a constrained collapse meaningful, vs the old skip-under-constraint).
  phy <- .SoftPolytomyData(5L)
  constraint <- ape::read.tree(
    text = "((f1, f2), (b1, b2, b3, b4, f3, f4, f5));")
  set.seed(1L)
  collapsed <- MaximizeParsimony(
    phy, constraint = constraint, concavity = Inf, maxReplicates = 12L,
    strategy = "intensive",
    control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
    verbosity = 0L, collapse = TRUE)

  isClade <- function(tr, members) {
    tr <- TreeTools::RootTree(tr, "b1")
    mrca <- ape::getMRCA(tr, members)
    !is.null(mrca) &&
      setequal(tr[["tip.label"]][phangorn::Descendants(tr, mrca, "tips")[[1]]],
               members)
  }
  # The enforced (f1, f2) clade is shown in every returned tree ...
  expect_true(all(vapply(collapsed, isClade, logical(1),
                         members = c("f1", "f2"))))
  # ... yet the rest of the fan is still contracted (tree is a polytomy), so the
  # collapse acted everywhere except on the protected constraint split.
  expect_true(all(vapply(collapsed, function(t)
    t[["Nnode"]] < NTip(phy) - 1L, logical(1))))
})

test_that("collapse = TRUE contracts a fully-uninformative dataset to a star (T-331)", {
  # When every character is parsimony-uninformative (all-constant here), every
  # binary resolution ties at the same (zero) score, so the collapsed answer is
  # a single star: one polytomy, one topology -- never the inflated binary
  # count that a total_words == 0 no-op would produce.
  taxa <- paste0("t", 1:5)
  mat <- matrix("0", nrow = length(taxa), ncol = 2,
                dimnames = list(taxa, NULL))
  phy <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))

  collapsed <- MaximizeParsimony(phy, verbosity = 0L, collapse = TRUE)

  expect_equal(length(collapsed), 1L)
  expect_equal(attr(collapsed, "n_topologies"), 1L)
  expect_equal(collapsed[[1]][["Nnode"]], 1L)
})
