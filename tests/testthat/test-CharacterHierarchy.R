test_that("CharacterHierarchy constructor works", {
  # Simple

  h <- CharacterHierarchy("1" = 2:5)
  expect_s3_class(h, "CharacterHierarchy")
  expect_length(h, 1)
  expect_equal(h[[1]]$controlling, 1L)
  expect_equal(h[[1]]$dependents, 2:5)
  expect_length(h[[1]]$children, 0)

  # Multiple blocks
  h2 <- CharacterHierarchy("1" = 2:5, "10" = 11:12)
  expect_length(h2, 2)
  expect_equal(h2[[2]]$controlling, 10L)
  expect_equal(h2[[2]]$dependents, 11:12)

  # Nested hierarchy
  h3 <- CharacterHierarchy("1" = list(2, 4, 5, "3" = 9:10))
  expect_length(h3, 1)
  expect_equal(h3[[1]]$controlling, 1L)
  expect_true(3L %in% h3[[1]]$dependents)
  expect_length(h3[[1]]$children, 1)
  expect_equal(h3[[1]]$children[[1]]$controlling, 3L)
  expect_equal(h3[[1]]$children[[1]]$dependents, 9:10)

  # A sub-controller listed explicitly AND as a named sub-block (e.g. the
  # documented `list(2, 3, 4, 5, "3" = 9:10)`) must appear once, not twice,
  # in `dependents` (RTS-003).
  h4 <- CharacterHierarchy("1" = list(2, 3, 4, 5, "3" = 9:10))
  expect_equal(h4[[1]]$dependents, c(2L, 3L, 4L, 5L))
  expect_false(anyDuplicated(h4[[1]]$dependents) > 0)
})

test_that("CharacterHierarchy rejects bad input", {
  expect_error(CharacterHierarchy(), "At least one")
  expect_error(CharacterHierarchy(2:5), "must be named")
  expect_error(CharacterHierarchy("abc" = 2:3), "integer indices")
})

test_that("print.CharacterHierarchy runs", {
  h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
  expect_output(print(h), "Char 1 controls")
  expect_output(print(h), "Char 6 controls")
})

test_that("HierarchyChars extracts all indices", {
  h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
  chars <- HierarchyChars(h)
  expect_setequal(chars, 1:8)

  # Nested
  h2 <- CharacterHierarchy("1" = list(2, 3, "3" = 9:10))
  chars2 <- HierarchyChars(h2)
  expect_setequal(chars2, c(1, 2, 3, 9, 10))
})

test_that("HierarchyControlling returns top-level controllers", {
  h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
  expect_equal(HierarchyControlling(h), c(1L, 6L))
})

test_that("HierarchyFromNames parses TNT-style names", {
  nms <- c("sup_tail", "sub_tail_colour", "sub_tail_shape",
           "sup_wing", "sub_wing_venation", "eyes")
  h <- HierarchyFromNames(nms)
  expect_s3_class(h, "CharacterHierarchy")
  expect_setequal(HierarchyControlling(h), c(1L, 4L))
  expect_setequal(HierarchyChars(h), c(1, 2, 3, 4, 5))
})

test_that("HierarchyFromNames returns NULL with no hierarchy", {
  nms <- c("eyes", "legs", "wings")
  expect_null(HierarchyFromNames(nms))
})

test_that("HierarchyFromNames warns on orphan sub_ tags", {
  nms <- c("sup_tail", "sub_tail_colour", "sub_arm_length")
  expect_warning(HierarchyFromNames(nms), "no corresponding sup_")
})

test_that("ValidateHierarchy passes on well-formed data", {
  mat <- matrix(c(
    "0", "-", "-", "1",
    "0", "-", "-", "0",
    "1", "0", "1", "1",
    "1", "1", "0", "0"
  ), nrow = 4, byrow = TRUE)
  rownames(mat) <- LETTERS[1:4]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2:3)
  expect_silent(ValidateHierarchy(h, ds))
})

test_that("ValidateHierarchy catches non-inapplicable secondaries", {
  mat <- matrix(c(
    "0", "1", "1",
    "1", "0", "1",
    "1", "1", "0"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2:3)
  expect_error(ValidateHierarchy(h, ds), "non-inapplicable")
})

test_that("ValidateHierarchy catches non-binary controlling character", {
  mat <- matrix(c(
    "0", "-",
    "1", "0",
    "2", "1"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1", "2"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2L)
  expect_error(ValidateHierarchy(h, ds), "binary")
})

test_that("ValidateHierarchy catches out-of-range indices", {
  mat <- matrix(c("0", "-", "1", "0"), nrow = 2, byrow = TRUE)
  rownames(mat) <- c("A", "B")
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 99L)
  expect_error(ValidateHierarchy(h, ds), "out of range")
})

test_that("ValidateHierarchy catches double-claimed characters", {
  mat <- matrix(c(
    "0", "-", "-", "0", "-",
    "1", "0", "1", "1", "0",
    "1", "1", "0", "0", "-"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2:3, "4" = c(3L, 5L))
  expect_error(ValidateHierarchy(h, ds), "multiple hierarchy blocks")
})

test_that(".NonHierarchyWeights subtracts hierarchy chars", {
  mat <- matrix(c(
    "0", "-", "0", "1",
    "1", "0", "1", "0",
    "1", "1", "0", "0"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")

  h <- CharacterHierarchy("1" = 2L)
  w_orig <- attr(ds, "weight")
  w_adj <- TreeSearch:::.NonHierarchyWeights(ds, h)

  # Adjusted weights should be non-negative

  expect_true(all(w_adj >= 0L))
  # Total weight should decrease by the number of hierarchy chars
  expect_equal(sum(w_adj), sum(w_orig) - 2L)
})

test_that(".BuildTipLabels creates correct matrix", {
  mat <- matrix(c(
    "0", "-", "1",
    "1", "0", "0",
    "0", "1", "1"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")

  tl <- TreeSearch:::.BuildTipLabels(ds)
  expect_equal(nrow(tl), 3L)
  expect_equal(ncol(tl), 3L)
  # Values should be 0-based token indices
  expect_true(all(tl >= 0L))
})

test_that(".HierarchyToBlocks converts to 0-based flat list", {
  h <- CharacterHierarchy("1" = 2:4, "5" = 6:7)
  blocks <- TreeSearch:::.HierarchyToBlocks(h)
  expect_length(blocks, 2)
  expect_equal(blocks[[1]]$primary, 0L)
  expect_equal(blocks[[1]]$secondaries, 1:3)
  expect_equal(blocks[[2]]$primary, 4L)
  expect_equal(blocks[[2]]$secondaries, 5:6)
})

test_that(".HierarchyToBlocks flattens nested hierarchies", {
  h <- CharacterHierarchy("1" = list(2, 4, "3" = 9:10))
  blocks <- TreeSearch:::.HierarchyToBlocks(h)
  expect_gte(length(blocks), 2)
  # First block: primary=0, secondaries should include 1 and 3 (chars 2 and 4)
  expect_equal(blocks[[1]]$primary, 0L)
  # Nested block: primary=2 (char 3), secondaries=8:9 (chars 9, 10)
  nested <- blocks[[2]]
  expect_equal(nested$primary, 2L)
  expect_equal(nested$secondaries, c(8L, 9L))
})

test_that(".NonHierarchyWeights preserves non-hierarchy patterns", {
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "0", "1", "0"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")

  h <- CharacterHierarchy("1" = 2L)
  idx <- attr(ds, "index")
  w_orig <- attr(ds, "weight")
  w_adj <- TreeSearch:::.NonHierarchyWeights(ds, h)

  # Character 3 is not in the hierarchy; its pattern should keep its weight
  # unless it shares a pattern with a hierarchy character
  non_h_chars <- setdiff(seq_along(idx), HierarchyChars(h))
  for (ci in non_h_chars) {
    pat <- idx[ci]
    # Pattern weight should be at least 1 for non-hierarchy chars
    # (could be reduced if shared with a hierarchy char)
    expect_gte(w_adj[pat], 0L)
  }
})


# =========================================================================
# Nested hierarchies: the documented example must actually validate (T-395)
# =========================================================================
# .ParseOneBlock() deliberately records a sub-controller BOTH as a dependent of
# its parent and as the controlling character of its own block -- that dual role
# is the entire content of "nested".  .ValidateBlock() used to count the second
# occurrence as a double claim, so every nested hierarchy that could be written
# was rejected, including this file's own roxygen example.  R CMD check could not
# see it because CharacterHierarchy() itself never validates: the example
# constructs fine and is never scored.

test_that("a nested hierarchy validates, and the documented example works", {
  # Char 1 controls {2, 3, 4, 5}; char 3 additionally controls {9, 10}.
  nested <- CharacterHierarchy("1" = list(2, 3, 4, 5, "3" = 9:10))

  expect_equal(nested[[1]]$controlling, 1L)
  expect_equal(sort(nested[[1]]$dependents), c(2L, 3L, 4L, 5L))
  expect_equal(length(nested[[1]]$children), 1L)
  expect_equal(nested[[1]]$children[[1]]$controlling, 3L)
  expect_equal(nested[[1]]$children[[1]]$dependents, 9:10)

  # Coding invariants: secondaries "-" where their controller codes absence;
  # char 3 binary where it applies; chars 9-10 "-" wherever char 3 is not "1".
  mat <- rbind(
    t1 = c("1", "0", "1", "0", "1", "0", "1", "0", "0", "1"),
    t2 = c("1", "1", "1", "1", "0", "1", "0", "1", "1", "1"),
    t3 = c("1", "0", "0", "1", "1", "0", "1", "1", "-", "-"),
    t4 = c("0", "-", "-", "-", "-", "1", "0", "0", "-", "-"),
    t5 = c("0", "-", "-", "-", "-", "0", "1", "1", "-", "-"),
    t6 = c("0", "-", "-", "-", "-", "1", "1", "0", "-", "-"))
  ds <- phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                         ambiguity = "?")

  expect_silent(ValidateHierarchy(nested, ds))

  # A genuine double claim must STILL be rejected -- the fix must not have
  # simply disabled the check.  Char 2 is claimed by both blocks here.
  expect_error(
    ValidateHierarchy(CharacterHierarchy("1" = 2:3, "5" = c(2L, 6L)), ds),
    "multiple hierarchy blocks"
  )
  # And a nested block's own dependents are still checked against other blocks.
  expect_error(
    ValidateHierarchy(CharacterHierarchy("1" = list(2, "3" = 9:10), "5" = 9L),
                      ds),
    "multiple hierarchy blocks"
  )
})

test_that("a nested hierarchy scores, and stays rooting-invariant", {
  nested <- CharacterHierarchy("1" = list(2, 3, 4, 5, "3" = 9:10))
  mat <- rbind(
    t1 = c("1", "0", "1", "0", "1", "0", "1", "0", "0", "1"),
    t2 = c("1", "1", "1", "1", "0", "1", "0", "1", "1", "1"),
    t3 = c("1", "0", "0", "1", "1", "0", "1", "1", "-", "-"),
    t4 = c("0", "-", "-", "-", "-", "1", "0", "0", "-", "-"),
    t5 = c("0", "-", "-", "-", "-", "0", "1", "1", "-", "-"),
    t6 = c("0", "-", "-", "-", "-", "1", "1", "0", "-", "-"))
  ds <- phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                         ambiguity = "?")
  tree <- Renumber(RenumberTips(
    ape::read.tree(text = "(((t1,t2),t3),(t4,(t5,t6)));"), rownames(mat)))

  score <- TreeLength(tree, ds, hierarchy = nested, inapplicable = "hsj",
                      hsj_alpha = 1)
  expect_true(is.finite(score))
  expect_gt(score, 0)

  # The HSJ score is a minimum over labellings of a sum of SYMMETRIC
  # dissimilarities on an unrooted tree, so it cannot depend on the rooting.
  rooted <- vapply(rownames(mat), function(tip) {
    TreeLength(Renumber(RootTree(tree, tip)), ds, hierarchy = nested,
               inapplicable = "hsj", hsj_alpha = 1)
  }, double(1))
  expect_equal(diff(range(rooted)), 0)

  # The x-transformation genuinely does not implement nesting; its error must
  # be the informative one, not a double-claim complaint from the validator.
  expect_error(RecodeHierarchy(ds, nested), "[Nn]ested")
})

test_that("HierarchyFromNames detects nesting by tag extension", {
  # Flat, documented case: unchanged.
  flat <- HierarchyFromNames(c("sup_tail", "sub_tail_colour", "sub_tail_shape",
                               "sup_wing", "sub_wing_venation", "eyes"))
  expect_equal(length(flat), 2L)
  expect_equal(sort(flat[[1]]$dependents), c(2L, 3L))
  expect_equal(length(flat[[1]]$children), 0L)

  # Nested: `tail_tip` extends `tail`, so char 3 is both a dependent of char 1
  # and the controller of char 4.  The old first-component tag match collapsed
  # every depth onto the outermost tag, so this returned two flat blocks -- one
  # of them empty -- and the nesting was silently lost.
  nested <- HierarchyFromNames(c("sup_tail", "sub_tail_colour",
                                 "sup_tail_tip", "sub_tail_tip_gloss"))
  expect_equal(length(nested), 1L)
  expect_equal(nested[[1]]$controlling, 1L)
  expect_equal(sort(nested[[1]]$dependents), c(2L, 3L))
  expect_equal(length(nested[[1]]$children), 1L)
  expect_equal(nested[[1]]$children[[1]]$controlling, 3L)
  expect_equal(nested[[1]]$children[[1]]$dependents, 4L)

  # Longest-match, not first-match: the deeper dependent must not be captured
  # by the shallower tag.
  expect_false(4L %in% nested[[1]]$dependents)

  # A shared prefix without an underscore boundary is NOT nesting.
  fin <- HierarchyFromNames(c("sup_tail", "sup_tailfin", "sub_tailfin_x"))
  expect_equal(length(fin), 2L)

  # Three levels deep.
  deep <- HierarchyFromNames(c("sup_a", "sub_a_x", "sup_a_b", "sub_a_b_y",
                               "sup_a_b_c", "sub_a_b_c_z"))
  expect_equal(length(deep), 1L)
  expect_equal(deep[[1]]$children[[1]]$controlling, 3L)
  expect_equal(deep[[1]]$children[[1]]$children[[1]]$controlling, 5L)
  expect_equal(deep[[1]]$children[[1]]$children[[1]]$dependents, 6L)

  # An orphan sub_ still warns.
  expect_warning(HierarchyFromNames(c("sup_tail", "sub_nose_shape")),
                 "no corresponding sup_")
})
