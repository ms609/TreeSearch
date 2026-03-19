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

test_that("hierarchy_chars extracts all indices", {
  h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
  chars <- hierarchy_chars(h)
  expect_setequal(chars, 1:8)

  # Nested
  h2 <- CharacterHierarchy("1" = list(2, 3, "3" = 9:10))
  chars2 <- hierarchy_chars(h2)
  expect_setequal(chars2, c(1, 2, 3, 9, 10))
})

test_that("hierarchy_controlling returns top-level controllers", {
  h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
  expect_equal(hierarchy_controlling(h), c(1L, 6L))
})

test_that("hierarchy_from_names parses TNT-style names", {
  nms <- c("sup_tail", "sub_tail_colour", "sub_tail_shape",
           "sup_wing", "sub_wing_venation", "eyes")
  h <- hierarchy_from_names(nms)
  expect_s3_class(h, "CharacterHierarchy")
  expect_setequal(hierarchy_controlling(h), c(1L, 4L))
  expect_setequal(hierarchy_chars(h), c(1, 2, 3, 4, 5))
})

test_that("hierarchy_from_names returns NULL with no hierarchy", {
  nms <- c("eyes", "legs", "wings")
  expect_null(hierarchy_from_names(nms))
})

test_that("hierarchy_from_names warns on orphan sub_ tags", {
  nms <- c("sup_tail", "sub_tail_colour", "sub_arm_length")
  expect_warning(hierarchy_from_names(nms), "no corresponding sup_")
})

test_that("validate_hierarchy passes on well-formed data", {
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
  expect_silent(validate_hierarchy(h, ds))
})

test_that("validate_hierarchy catches non-inapplicable secondaries", {
  mat <- matrix(c(
    "0", "1", "1",
    "1", "0", "1",
    "1", "1", "0"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2:3)
  expect_error(validate_hierarchy(h, ds), "non-inapplicable")
})

test_that("validate_hierarchy catches non-binary controlling character", {
  mat <- matrix(c(
    "0", "-",
    "1", "0",
    "2", "1"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1", "2"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2L)
  expect_error(validate_hierarchy(h, ds), "binary")
})

test_that("validate_hierarchy catches out-of-range indices", {
  mat <- matrix(c("0", "-", "1", "0"), nrow = 2, byrow = TRUE)
  rownames(mat) <- c("A", "B")
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 99L)
  expect_error(validate_hierarchy(h, ds), "out of range")
})

test_that("validate_hierarchy catches double-claimed characters", {
  mat <- matrix(c(
    "0", "-", "-", "0", "-",
    "1", "0", "1", "1", "0",
    "1", "1", "0", "0", "-"
  ), nrow = 3, byrow = TRUE)
  rownames(mat) <- LETTERS[1:3]
  ds <- phangorn::phyDat(mat, type = "USER",
                         levels = c("-", "0", "1"), ambiguity = "?")
  h <- CharacterHierarchy("1" = 2:3, "4" = c(3L, 5L))
  expect_error(validate_hierarchy(h, ds), "multiple hierarchy blocks")
})
