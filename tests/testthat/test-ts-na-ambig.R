# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# {inapplicable, state} = {-,X} partial ambiguity under BGS (Brazeau, Guillerme
# & Smith 2019 three-pass inapplicable parsimony).
#
# A {-,X} tip is ambiguous on the binary APPLICABILITY character (0 = inapplic-
# able, 1 = applicable): it is a "?" there.  BGS reconstructs applicability by
# Fitch, resolving each {-,X} tip in tree context — applicable-PREFERRED on a
# tie — and only then counts state steps.  This MAXIMISES HOMOLOGY (De Laet)
# rather than minimising homoplasy: an observed state is a homology statement
# that must be paid for, so a {-,X} tip is NOT free to become inapplicable
# merely to dodge a state step.

# Helper: build a one-character phyDat with {-,S} ambiguity tokens.
make_ambig_phyDat <- function(tip_states, tip_labels, levels, ambig_tokens) {
  all_tokens <- c(levels, names(ambig_tokens))
  n_states <- length(levels)
  contrast <- matrix(0, nrow = length(all_tokens), ncol = n_states,
                     dimnames = list(all_tokens, levels))
  for (i in seq_along(levels)) contrast[levels[i], levels[i]] <- 1
  for (nm in names(ambig_tokens)) {
    contrast[nm, ambig_tokens[[nm]]] <- 1
  }
  phangorn::phyDat(
    matrix(tip_states, ncol = 1, dimnames = list(tip_labels, NULL)),
    type = "USER", levels = levels,
    ambiguity = names(ambig_tokens), contrast = contrast
  )
}


test_that("{-,X} resolves APPLICABLE when applicability supports it", {
  library("TreeTools")
  library("TreeSearch")
  lvls2 <- c("-", "1", "2")
  ambig2 <- list("{-2}" = c("-", "2"))

  # ((t1, t2), (t3, t4)) with 1 1 {-,2} 1: t3 sits among applicable "1" tips, so
  # the applicability character forces it applicable -> state 2 -> one step.
  # Stripping t3 to "-" would give 0; the homology of the observed "2" costs.
  tree4 <- read.tree(text = "((t1, t2), (t3, t4));")
  dat_11m21 <- make_ambig_phyDat(c("1", "1", "{-2}", "1"),
                                 tree4$tip.label, lvls2, ambig2)
  dat_11_1  <- make_ambig_phyDat(c("1", "1", "-", "1"),
                                 tree4$tip.label, lvls2, ambig2)
  expect_equal(TreeLength(tree4, dat_11m21), 1L)
  expect_equal(TreeLength(tree4, dat_11_1), 0L)
  expect_gt(TreeLength(tree4, dat_11m21), TreeLength(tree4, dat_11_1))

  # Larger balanced tree, same mechanism: 1 1 {-,2} 1 3 3 -> {-2} applicable.
  tree6 <- read.tree(text = "((t1, t2), ((t3, t4), (t5, t6)));")
  lvls3 <- c("-", "1", "2", "3")
  dat_6tip    <- make_ambig_phyDat(c("1", "1", "{-2}", "1", "3", "3"),
                                   tree6$tip.label, lvls3, ambig2)
  dat_6tip_na <- make_ambig_phyDat(c("1", "1", "-", "1", "3", "3"),
                                   tree6$tip.label, lvls3, ambig2)
  expect_equal(TreeLength(tree6, dat_6tip), 2L)
  expect_equal(TreeLength(tree6, dat_6tip_na), 1L)
})


test_that("forced-applicable {-,X} (cherry with an applicable sibling)", {
  library(TreeTools)
  # ((t1, t2), ((t3, t5), t4)); t1 = {-,2} is in a cherry with the applicable
  # t2 = 1, so applicability is FORCED to applicable at (t1,t2) -> t1 = state 2.
  tree <- read.tree(text = "((t1, t2), ((t3, t5), t4));")
  dat <- make_ambig_phyDat(c("{-2}", "1", "-", "2", "1"),
                           c("t1", "t2", "t3", "t4", "t5"),
                           c("-", "1", "2"), list("{-2}" = c("-", "2")))
  expect_equal(TreeLength(tree, dat), 2L)
})


test_that("applicability tie breaks toward applicable, not the cheaper -", {
  library(TreeTools)
  # (t1,(t2,(t3,(t4,t5)))); 1 2 1 {-,2} - : making t4 applicable vs inapplicable
  # is EQUALLY parsimonious on the applicability character (a genuine tie), so
  # BGS takes the applicable resolution (maximise homology) -> t4 = 2 is counted
  # -> length 2, even though resolving t4 = "-" would score 1.  This is the
  # deliberate consequence of Fitch-optimising applicability first.
  tree <- read.tree(text = "(t1, (t2, (t3, (t4, t5))));")
  dat <- make_ambig_phyDat(c("1", "2", "1", "{-2}", "-"),
                           c("t1", "t2", "t3", "t4", "t5"),
                           c("-", "1", "2"), list("{-2}" = c("-", "2")))
  expect_equal(TreeLength(tree, dat), 2L)
})


test_that("{-,X} resolves INAPPLICABLE when its neighbourhood is inapplicable", {
  library(TreeTools)
  tree <- PectinateTree(7)
  tree$tip.label <- paste0("t", 1:7)
  lvls <- c("-", "1", "2", "3")
  ambig <- list("{-1}" = c("-", "1"), "{-2}" = c("-", "2"), "{-3}" = c("-", "3"))

  # {-1} - - 2 2 3 3: the {-1} tip sits beside two inapplicable tips, so the
  # applicability character resolves it inapplicable -> same as - - - 2 2 3 3.
  dat_one     <- make_ambig_phyDat(c("{-1}", "-", "-", "2", "2", "3", "3"),
                                   tree$tip.label, lvls, ambig)
  dat_pure_na <- make_ambig_phyDat(c("-", "-", "-", "2", "2", "3", "3"),
                                   tree$tip.label, lvls, ambig)
  expect_equal(TreeLength(tree, dat_one), 1L)
  expect_equal(TreeLength(tree, dat_one), TreeLength(tree, dat_pure_na))
})


test_that("same-state {-,X} tips join the applicable region for free", {
  library(TreeTools)
  tree <- PectinateTree(7)
  tree$tip.label <- paste0("t", 1:7)
  lvls <- c("-", "1", "2", "3")
  ambig <- list("{-1}" = c("-", "1"),
                "{-2}" = c("-", "2"),
                "{-3}" = c("-", "3"))

  # Three {-,2} tips beside 2 2: all resolve applicable AS state 2, homologous
  # with the observed 2s -> no extra step (length 1, = - - - 2 2 3 3).
  dat_same <- make_ambig_phyDat(c("{-2}", "{-2}", "{-2}", "2", "2", "3", "3"),
                                tree$tip.label, lvls, ambig)
  expect_equal(TreeLength(tree, dat_same), 1L)

  # But {-,1}{-,2}{-,3} resolve to DISTINCT states 1,2,3 -> homoplasy is real
  # (length 3), NOT the 1 that stripping them all to "-" once gave.
  dat_distinct <- make_ambig_phyDat(c("{-1}", "{-2}", "{-3}", "2", "2", "3", "3"),
                                    tree$tip.label, lvls, ambig)
  expect_equal(TreeLength(tree, dat_distinct), 3L)
})


test_that("missing data (all states + NA) is not over-counted", {
  library(TreeTools)
  tree <- PectinateTree(7)
  tree$tip.label <- paste0("t", 1:7)
  lvls <- c("-", "1", "2", "3")
  ambig <- list(
    "{-1}" = c("-", "1"), "{-2}" = c("-", "2"), "{-3}" = c("-", "3"),
    "miss" = c("-", "1", "2", "3")
  )

  dat_miss <- make_ambig_phyDat(
    c("miss", "miss", "miss", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  dat_na <- make_ambig_phyDat(
    c("-", "-", "-", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  expect_lte(TreeLength(tree, dat_miss), TreeLength(tree, dat_na))
})


test_that("pure applicable characters unaffected", {
  library(TreeTools)
  tree <- PectinateTree(7)
  tree$tip.label <- paste0("t", 1:7)
  lvls <- c("-", "1", "2", "3")
  ambig <- list("{-1}" = c("-", "1"), "{-2}" = c("-", "2"),
                "{-3}" = c("-", "3"))

  dat <- make_ambig_phyDat(
    c("1", "2", "3", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  expect_equal(TreeLength(tree, dat), 3L)
})


test_that("pure inapplicable characters unaffected by {-,X} handling", {
  library(TreeTools)

  # Regression: standard NA cases (no {-,X} ambiguity) are untouched by the
  # faithful {-,X} encoding.
  tree <- PectinateTree(7)
  tree$tip.label <- paste0("t", 1:7)
  lvls <- c("-", "1", "2", "3")
  ambig <- list("{-1}" = c("-", "1"))

  cases <- list(
    list(data = c("-", "-", "-", "2", "2", "3", "3"), score = 1L,
         label = "---2233"),
    list(data = c("-", "1", "-", "2", "2", "2", "-"), score = 1L,
         label = "-1-222-"),
    list(data = c("1", "1", "1", "2", "2", "2", "2"), score = 1L,
         label = "1112222"),
    list(data = c("-", "-", "-", "-", "-", "-", "-"), score = 0L,
         label = "all inapplicable")
  )
  for (case in cases) {
    dat <- make_ambig_phyDat(case$data, tree$tip.label, lvls, ambig)
    expect_equal(TreeLength(tree, dat), case$score,
                 label = case$label)
  }
})
