test_that("All labellings counted (3 tokens)", {
  
  labels <- c(2, 3, 2)
  nTip <- sum(labels)
  tree <- TreeTools::BalancedTree(nTip)
  first <- labels[[1]]
  leftAfterFirst <- nTip - first
  second <- labels[[2]]
  nComb <- choose(nTip, first) * choose(leftAfterFirst, second)
  allComb <- apply(combn(nTip, first), 2, function(i) {
    apply(combn(leftAfterFirst, second), 2, function(j) {
      ret <- rep(1, nTip)
      ret2 <- rep(3, leftAfterFirst)
      ret2[j] <- 2
      ret[-i] <- ret2
      ret
    }, simplify = FALSE)
  }) |>
    unlist() |>
    paste0(collapse = "") |>
    StringToPhyDat(tips = nTip, byTaxon = FALSE)
  expect <- table(CharacterLength(tree, allComb, compress = TRUE))
  
  expect_equal(FixedTreeCount(tree, c(nTip), 2),
               log(c("0" = 1, "1" = 0, "2" = 0)))
  expect_equal(FixedTreeCount(tree, c(nTip - 1, 1), 2),
               log(c("0" = 0, "1" = nTip, "2" = 0)))
               
  expect_equal(FixedTreeCount(tree, c(2, 3, 2), 5),
               log(c("0" = 0, "1" = 0, expect, "5" = 0)))
  
  expect_equal(
    FixedTreeCountBatch(tree, cbind(c(nTip, 0, 0), c(nTip - 1, 0, 1), labels)) |>
      exp(),
    cbind(c(1, rep(0, 5)),
          c(0, 7, rep(0, 4)),
          c("0" = 0, "1" = 0, expect, "5" = 0))
  )
})

test_that("All labellings counted - 4 tokens", {
  
  labels <- c(2, 4, 2, 2)
  nTip <- sum(labels)
  tree <- TreeTools::BalancedTree(nTip)
  expect <- if (recalculateExpectation <- FALSE) {
    first <- labels[[1]]
    leftAfter1 <- nTip - first
    second <- labels[[2]]
    leftAfter2 <- leftAfter1 - second
    third <- labels[[3]]
    nComb <- prod(choose(nTip, first),
                  choose(leftAfter1, second),
                  choose(leftAfter2, third))
    allComb <- apply(combn(nTip, first), 2, function(i) {
      apply(combn(leftAfter1, second), 2, function(j) {
        apply(combn(leftAfter2, third), 2, function(k) {
          ret <- rep(1, nTip)
          ret2 <- rep(2, leftAfter1)
          ret3 <- rep(3, leftAfter2)
          ret3[-k] <- 4
          ret2[-j] <- ret3
          ret[-i] <- ret2
          ret
        }, simplify = FALSE)
      })
    }) |>
      unlist() |>
      paste0(collapse = "") |>
      StringToPhyDat(tips = nTip, byTaxon = FALSE)
    table(CharacterLength(tree, allComb, compress = TRUE))
  } else {
    c(`3` = 24L, `4` = 972L, `5` = 6096L, `6` = 11808L)
  }
  
  expect_equal(exp(FixedTreeCount(tree, c(nTip), 2)),
               c("0" = 1, "1" = 0, "2" = 0))
  expect_equal(exp(FixedTreeCount(tree, c(nTip - 1, 1), 2)),
               c("0" = 0, "1" = nTip, "2" = 0))
               
  expect_equal(exp(FixedTreeCount(tree, labels)[as.character(3:6)]), expect)
})
