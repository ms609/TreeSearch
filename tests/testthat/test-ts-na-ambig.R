# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Test that {inapplicable, state} ambiguity tokens are handled correctly
# by stripping applicable bits when the inapplicable bit is present.
#
# When a tip has both the inapplicable bit and an applicable-state bit
# (e.g. {-,2}), the applicable bits are stripped in build_dataset(),
# treating the tip as pure inapplicable. This matches MorphyLib's
# SingleCharMorphy behaviour. The three-pass NA algorithm cannot
# correctly resolve {-,X} ambiguity in tree context.

# Helper: build a phyDat with {-,S} ambiguity tokens
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


test_that("{-,S} resolves as inapplicable when that is more parsimonious", {
  library(TreeSearch)
  library(TreeTools)
  tree <- PectinateTree(7)
  tree$tip.label <- paste0("t", 1:7)
  lvls <- c("-", "1", "2", "3")
  ambig <- list("{-1}" = c("-", "1"), "{-2}" = c("-", "2"),
                "{-3}" = c("-", "3"))

  # {-1}{-2}{-3}2233: all three ambig tips share inapplicable; resolving
  # as inapplicable gives score 1 (= score of ---2233).
  dat_ambig <- make_ambig_phyDat(
    c("{-1}", "{-2}", "{-3}", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  dat_pure_na <- make_ambig_phyDat(
    c("-", "-", "-", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  expect_equal(TreeLength(tree, dat_ambig), 1L)
  expect_equal(TreeLength(tree, dat_ambig), TreeLength(tree, dat_pure_na))

  # Variants: same-state ambig tokens, single ambig token
  dat_same <- make_ambig_phyDat(
    c("{-2}", "{-2}", "{-2}", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  expect_equal(TreeLength(tree, dat_same), TreeLength(tree, dat_pure_na))

  dat_one <- make_ambig_phyDat(
    c("{-1}", "-", "-", "2", "2", "3", "3"),
    tree$tip.label, lvls, ambig
  )
  expect_equal(TreeLength(tree, dat_one), TreeLength(tree, dat_pure_na))
})


test_that("{-,S} is stripped to pure inapplicable", {
  library(TreeTools)

  # With the strip approach, {-,X} tokens collapse to pure inapplicable.
  # The three-pass NA algorithm cannot correctly resolve {-,X} ambiguity
  # in tree context, so we conservatively treat them as inapplicable.
  tree4 <- read.tree(text = "((t1, t2), (t3, t4));")
  lvls2 <- c("-", "1", "2")
  ambig2 <- list("{-2}" = c("-", "2"))

  dat_11m21 <- make_ambig_phyDat(
    c("1", "1", "{-2}", "1"),
    tree4$tip.label, lvls2, ambig2
  )
  dat_11_1 <- make_ambig_phyDat(
    c("1", "1", "-", "1"),
    tree4$tip.label, lvls2, ambig2
  )
  # {-2} stripped to pure inapplicable: same score as 11-1
  expect_equal(TreeLength(tree4, dat_11m21), TreeLength(tree4, dat_11_1))

  # Also test on a balanced tree with more tips
  tree6 <- read.tree(text = "((t1, t2), ((t3, t4), (t5, t6)));")
  lvls3 <- c("-", "1", "2", "3")
  ambig3 <- list("{-2}" = c("-", "2"))
  dat_6tip <- make_ambig_phyDat(
    c("1", "1", "{-2}", "1", "3", "3"),
    tree6$tip.label, lvls3, ambig3
  )
  dat_6tip_na <- make_ambig_phyDat(
    c("1", "1", "-", "1", "3", "3"),
    tree6$tip.label, lvls3, ambig3
  )
  # {-2} stripped to pure inapplicable: same score as 11-133
  expect_equal(TreeLength(tree6, dat_6tip), TreeLength(tree6, dat_6tip_na))
})


test_that("missing data (all states + NA) is not collapsed", {
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


test_that("pure inapplicable characters unaffected by {-,S} fix", {
  library(TreeTools)

  # Regression: ensure the NA algorithm still works for standard cases
  # (no {-,X} ambiguity tokens) after the fix.
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
