# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for `inapplicable = "missing"` (pure-Fitch mode), which recodes every
# gap-bearing token as missing data via .GapsAsMissing().  The distinguishing
# behaviour from a naive "-"->"?" text substitution is that {state, -}
# polymorphisms must also become fully ambiguous ("?"), rather than surviving
# as genuine inapplicable tokens (which the engine strips to a pure gap).

library(TreeSearch)
library(TreeTools)

# Build a phyDat carrying {state, -} ambiguity tokens.
make_ambig_phyDat <- function(tip_states, tip_labels, levels, ambig_tokens) {
  all_tokens <- c(levels, names(ambig_tokens))
  n_states <- length(levels)
  contrast <- matrix(0, nrow = length(all_tokens), ncol = n_states,
                     dimnames = list(all_tokens, levels))
  for (i in seq_along(levels)) contrast[levels[i], levels[i]] <- 1
  for (nm in names(ambig_tokens)) contrast[nm, ambig_tokens[[nm]]] <- 1
  phangorn::phyDat(
    matrix(tip_states, ncol = 1, dimnames = list(tip_labels, NULL)),
    type = "USER", levels = levels,
    ambiguity = names(ambig_tokens), contrast = contrast
  )
}


test_that(".GapsAsMissing leaves gap-free data untouched", {
  m <- matrix(c("0", "0", "1", "1"), ncol = 1,
              dimnames = list(paste0("t", 1:4), NULL))
  ds <- MatrixToPhyDat(m)
  expect_false("-" %in% attr(ds, "levels"))
  expect_identical(TreeSearch:::.GapsAsMissing(ds), ds)
})


test_that(".GapsAsMissing recodes every gap-bearing token to fully ambiguous", {
  # White-box: this is the behaviour that separates "missing" from the old
  # text-substitution hack.  A {state, -} polymorphism must become "?" (all
  # states); the old hack left it as a genuine inapplicable token, which the
  # engine then stripped to a pure gap.
  ds <- make_ambig_phyDat(c("1", "1", "{-2}", "1"),
                          paste0("t", 1:4), c("-", "1", "2"),
                          list("{-2}" = c("-", "2")))
  rec <- TreeSearch:::.GapsAsMissing(ds)
  tok <- attr(rec, "allLevels")
  ct <- attr(rec, "contrast")
  # Pure gap "-" and the {-,2} polymorphism are now all-ones ("?").
  expect_true(all(ct[match("-", tok), ] == 1))
  expect_true(all(ct[match("{-2}", tok), ] == 1))
  # Applicable-only tokens are unchanged.
  expect_equal(unname(ct[match("1", tok), ]), c(0, 1, 0))
  expect_equal(unname(ct[match("2", tok), ]), c(0, 0, 1))
})


test_that("gaps-as-missing reduces to Fitch (exact hand-computed score)", {
  # Character 0,0,1,- on ((t1,t2),(t3,t4)); "-" -> "?" gives 0,0,1,{0,1}.
  # Fitch: cherry (t1,t2) = {0}; cherry (t3,t4) = {1}; root union = 1 step.
  tr4 <- ape::read.tree(text = "((t1, t2), (t3, t4));")
  m <- matrix(c("0", "0", "1", "-"), ncol = 1,
              dimnames = list(paste0("t", 1:4), NULL))
  ds <- MatrixToPhyDat(m)
  expect_equal(TreeLength(tr4, TreeSearch:::.GapsAsMissing(ds)), 1)
})


test_that("contrast recode matches raw '-'->'?' substitution on pure gaps", {
  # On pure-gap data the contrast-surgery path and a raw text substitution
  # (a wholly independent construction) must agree.
  m <- matrix(c("0", "0", "1", "-",
                "1", "-", "1", "0",
                "-", "-", "0", "1"), ncol = 3,
              dimnames = list(paste0("t", 1:4), NULL))
  ds <- MatrixToPhyDat(m)
  helperPath <- TreeSearch:::.GapsAsMissing(ds)
  sub <- m
  sub[sub == "-"] <- "?"
  subPath <- MatrixToPhyDat(sub)
  for (txt in c("((t1, t2), (t3, t4));", "((t1, t3), (t2, t4));")) {
    tr <- ape::read.tree(text = txt)
    expect_equal(TreeLength(tr, helperPath), TreeLength(tr, subPath))
  }
})


test_that("'missing' treats {state, -} as missing, not inapplicable", {
  # t3 = {-,1}; under "missing" this is "?" = {1,2}, so it scores identically
  # to coding t3 explicitly as "?".
  tr4 <- ape::read.tree(text = "((t1, t2), (t3, t4));")
  ds_amb <- make_ambig_phyDat(c("2", "2", "{-1}", "1"),
                              paste0("t", 1:4), c("-", "1", "2"),
                              list("{-1}" = c("-", "1")))
  mq <- matrix(c("2", "2", "?", "1"), ncol = 1,
               dimnames = list(paste0("t", 1:4), NULL))
  expect_equal(TreeLength(tr4, TreeSearch:::.GapsAsMissing(ds_amb)),
               TreeLength(tr4, MatrixToPhyDat(mq)))
})


test_that("'missing' differs from BGS on real inapplicable data", {
  # The mode must actually bite: on a dataset with genuine inapplicable
  # regions, gaps-as-missing drops the region cost that BGS pays.
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  tr <- RootTree(PectinateTree(ds), 1)
  bgs <- TreeLength(tr, ds, inapplicable = "bgs")
  missing <- TreeLength(tr, TreeSearch:::.GapsAsMissing(ds))
  expect_equal(bgs, 139)
  expect_equal(missing, 138)
  expect_lt(missing, bgs)
})


test_that("MaximizeParsimony(inapplicable = 'missing') wires through the recode", {
  # "missing" mode must be identical to searching a pre-recoded dataset under
  # the standard (BGS == Fitch on gap-free data) engine.  With a fixed seed
  # and a single thread the two searches are reproducible and must match.
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  set.seed(1)
  r1 <- MaximizeParsimony(ds, inapplicable = "missing",
                          maxReplicates = 2L, targetHits = 1L,
                          nThreads = 1L, verbosity = 0L)
  set.seed(1)
  r2 <- MaximizeParsimony(TreeSearch:::.GapsAsMissing(ds), inapplicable = "bgs",
                          maxReplicates = 2L, targetHits = 1L,
                          nThreads = 1L, verbosity = 0L)
  expect_true(is.finite(attr(r1, "score")))
  expect_equal(attr(r1, "score"), attr(r2, "score"))
})


# --- TreeLength support for inapplicable = "missing" ------------------------

test_that("TreeLength(inapplicable = 'missing') matches the explicit recode", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  tr <- RootTree(PectinateTree(ds), 1)
  expect_equal(TreeLength(tr, ds, inapplicable = "missing"),
               TreeLength(tr, TreeSearch:::.GapsAsMissing(ds)))
  expect_equal(TreeLength(tr, ds, inapplicable = "missing"), 138)
  expect_lt(TreeLength(tr, ds, inapplicable = "missing"),
            TreeLength(tr, ds, inapplicable = "bgs"))
})


test_that("TreeLength.list recodes after the taxon-subset round-trip", {
  # A tree on a taxon subset triggers .Recompress() (a phyDat round-trip)
  # inside TreeLength.list.  The recode must be applied AFTER that, or the
  # "-" level is restored and gaps revert to inapplicable scoring.
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  keep <- names(ds)[1:18]
  tr <- RootTree(PectinateTree(keep), 1)
  trees <- structure(list(tr, tr), class = "multiPhylo")
  bgs <- TreeLength(trees, ds, inapplicable = "bgs")
  missing <- TreeLength(trees, ds, inapplicable = "missing")
  expect_length(missing, 2L)
  expect_true(all(missing < bgs))
})


test_that("TreeLength.numeric threads 'missing' through to scoring", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  set.seed(7)
  viaMode <- TreeLength(3L, ds, inapplicable = "missing")
  set.seed(7)
  viaRecode <- TreeLength(3L, TreeSearch:::.GapsAsMissing(ds))
  expect_equal(viaMode, viaRecode)
})
