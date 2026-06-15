# `RandomMorphyTree()` uses C's MWC RNG (static global state in
# build_postorder.h), which is NOT seeded by set.seed().
# `RandomTreeScore()` builds a tree with TreeTools::RandomTree() (R's RNG) and
# scores it with the native Fitch kernel, so it IS reproducible via set.seed().
# The heavy distribution-of-scores tests below are statistical and slow; they
# run only under TREESEARCH_EXTENDED_TESTS.  Routine correctness of
# RandomTreeScore() is covered by test-RandomTreeScore.R.

ScoreObj <- function (char) {
  # `char` is a Morphy token string without a trailing ";"
  SingleCharMorphy(char)
}

test_that("four-tip trees are randomly distributed", {
  nTrees <- 36000
  stringency <- 1e-6
  nTip <- 4
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees, 1/(nTip - 1))
  rTrees <- vapply(logical(nTrees), function (XX)
    unlist(RandomMorphyTree(nTip)), integer((nTip * 4) - 3))
  expect_true(all(rTrees[1 + (seq_len(nTip - 1)), ] %in% nTip + seq_len(nTip - 2)))
  expect_lt(expectedBounds[1], sum(rTrees[2, ] == 5))
  expect_gt(expectedBounds[2], sum(rTrees[2, ] == 5))
  expect_lt(expectedBounds[1], sum(rTrees[3, ] == 5))
  expect_gt(expectedBounds[2], sum(rTrees[3, ] == 5))
  expect_lt(expectedBounds[1], sum(rTrees[4, ] == 5))
  expect_gt(expectedBounds[2], sum(rTrees[4, ] == 5))

  expect_true(all(table(rTrees[c(9, 12), ])[seq_len(nTip - 1)] > expectedBounds[1]))
  expect_true(all(table(rTrees[c(9, 12), ])[seq_len(nTip - 1)] < expectedBounds[2]))

  expect_true(all(table(rTrees[c(10, 13), ])[seq_len(nTip - 1)] < nTrees - expectedBounds[1]))
  expect_true(all(table(rTrees[c(10, 13), ])[seq_len(nTip - 1)] > nTrees - expectedBounds[2]))
})

test_that("four-tip trees are randomly scored", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("TREESEARCH_EXTENDED_TESTS"))),
              "Set TREESEARCH_EXTENDED_TESTS=true to run extended tests")
  set.seed(0)

  nTrees <- 12000
  stringency <- 1e-6
  nTip <- 4

  morphyObj <- ScoreObj("0011")

  expectedBounds <- qbinom(c(stringency, 1 - stringency), nTrees,
                           NUnrooted(nTip - 1L) / NUnrooted(nTip))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))
  expect_lt(expectedBounds[1], sum(scores==1))
  expect_gt(expectedBounds[2], sum(scores==1))
})

test_that("five-tip trees are randomly scored", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("TREESEARCH_EXTENDED_TESTS"))),
              "Set TREESEARCH_EXTENDED_TESTS=true to run extended tests")
  set.seed(0)
  nTrees <- 12000
  stringency <- 1e-6
  nTip <- 5
  morphyObj <- ScoreObj("00011")
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees,
                           NUnrooted(nTip - 1) / NUnrooted(nTip))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))
  expect_equal(2L, max(scores))
  expect_lt(expectedBounds[1], sum(scores == 1))
  expect_gt(expectedBounds[2], sum(scores == 1))
})


test_that("six-tip trees are randomly scored", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("TREESEARCH_EXTENDED_TESTS"))),
              "Set TREESEARCH_EXTENDED_TESTS=true to run extended tests")
  set.seed(0)

  nTrees <- 12000
  stringency <- 1e-6
  nTip <- 6

  morphyObj <- ScoreObj("000011")
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees,
                           NUnrooted(5) / NUnrooted(6))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))

  expect_true(max(scores) == 2)
  expect_lt(expectedBounds[1], sum(scores==1))
  expect_gt(expectedBounds[2], sum(scores==1))

  morphyObj <- ScoreObj("001122")
  expectedBounds <- qbinom(c(stringency, 1 - stringency), nTrees,
                           7 / NUnrooted(nTip))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj),
                   integer(1))

  expect_true(all(scores %in% 2:4))
  expect_lt(expectedBounds[1], sum(scores == 2))
  expect_gt(expectedBounds[2], sum(scores == 2))

  morphyObj <- ScoreObj("000111")
  expectedBounds <- qbinom(c(stringency, 1-stringency), nTrees,
                           3 * 3 / NUnrooted(nTip))
  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj), integer(1))

  expect_true(max(scores) == 3)
  expect_lt(expectedBounds[1], sum(scores == 1))
  expect_gt(expectedBounds[2], sum(scores == 1))

})

test_that("twelve-tip trees are randomly scored", {
  skip_if_not(isTRUE(as.logical(Sys.getenv("TREESEARCH_EXTENDED_TESTS"))),
              "Set TREESEARCH_EXTENDED_TESTS=true to run extended tests")
  set.seed(0)
  nTrees <- 24000
  stringency <- 1e-6
  nTip <- 12
  morphyObj <- ScoreObj("000000011111")
  expectedBounds <- qbinom(c(stringency, 1 - stringency), nTrees,
                           NUnrooted(7) * (2 * 7 - 3) *
                           NUnrooted(5) * (2 * 5 - 3) / NUnrooted(nTip))

  scores <- vapply(logical(nTrees),
                   function (XX) RandomTreeScore(morphyObj),
                   integer(1L))

  expect_equal(5L, max(scores))
  nScoring1 <- sum(scores == 1)
  expect_lt(expectedBounds[1], nScoring1)
  expect_gt(expectedBounds[2], nScoring1)
})
