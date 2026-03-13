# Helper to prepare phyDat for C++ engine
prep_pd <- function(pd) {
  list(
    contrast = attr(pd, "contrast"),
    tip_data = t(vapply(pd, I, pd[[1]])),
    weight = attr(pd, "weight"),
    levels = attr(pd, "levels")
  )
}

test_that("ts_wagner_tree produces valid tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  tip_data <- t(vapply(pd, I, pd[[1]]))

  result <- ts_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  expect_true(is.list(result))
  expect_equal(ncol(result$edge), 2L)
  expect_equal(nrow(result$edge), 2L * (n_tip - 1L))
  expect_true(result$score > 0)
})

test_that("Wagner score matches ts_fitch_score", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  tip_data <- t(vapply(pd, I, pd[[1]]))

  result <- ts_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  fitch_check <- ts_fitch_score(
    result$edge,
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  expect_equal(result$score, fitch_check)
})

test_that("Wagner score matches TreeLength", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  tip_data <- t(vapply(pd, I, pd[[1]]))

  result <- ts_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  result_tree <- structure(list(
    edge = result$edge,
    tip.label = names(pd),
    Nnode = n_tip - 1L
  ), class = "phylo")

  expect_equal(result$score, TreeLength(result_tree, pd))
})

test_that("All tips present exactly once", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  tip_data <- t(vapply(pd, I, pd[[1]]))

  result <- ts_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  # Tips are nodes 1..n_tip in the edge matrix
  child_nodes <- result$edge[, 2]
  tips_found <- sort(child_nodes[child_nodes <= n_tip])
  expect_equal(tips_found, seq_len(n_tip))
})

test_that("Same addition order gives same tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  tip_data <- t(vapply(pd, I, pd[[1]]))

  order_seq <- seq_len(n_tip)
  r1 <- ts_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels"),
    order_seq
  )
  r2 <- ts_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels"),
    order_seq
  )

  expect_identical(r1$edge, r2$edge)
  expect_identical(r1$score, r2$score)
})

test_that("Random Wagner trees vary", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  tip_data <- t(vapply(pd, I, pd[[1]]))

  set.seed(7263)
  r1 <- ts_random_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )
  set.seed(1984)
  r2 <- ts_random_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  expect_false(identical(r1$edge, r2$edge))
})

test_that("Random Wagner score verified", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  tip_data <- t(vapply(pd, I, pd[[1]]))

  set.seed(5511)
  result <- ts_random_wagner_tree(
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  fitch_check <- ts_fitch_score(
    result$edge,
    attr(pd, "contrast"), tip_data,
    attr(pd, "weight"), attr(pd, "levels")
  )

  expect_equal(result$score, fitch_check)
})

test_that("Small tree (5 tips) is correct", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd5 <- congreveLamsdellMatrices[[1]][1:5]
  tip_data <- t(vapply(pd5, I, pd5[[1]]))

  result <- ts_wagner_tree(
    attr(pd5, "contrast"), tip_data,
    attr(pd5, "weight"), attr(pd5, "levels")
  )

  expect_equal(nrow(result$edge), 8L)
  expect_true(result$score > 0)

  fitch_check <- ts_fitch_score(
    result$edge,
    attr(pd5, "contrast"), tip_data,
    attr(pd5, "weight"), attr(pd5, "levels")
  )
  expect_equal(result$score, fitch_check)
})

test_that("Medium tree (20 tips) completes without error", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd20 <- congreveLamsdellMatrices[[1]][1:20]
  tip_data <- t(vapply(pd20, I, pd20[[1]]))

  expect_no_error({
    result <- ts_wagner_tree(
      attr(pd20, "contrast"), tip_data,
      attr(pd20, "weight"), attr(pd20, "levels")
    )
  })

  fitch_check <- ts_fitch_score(
    result$edge,
    attr(pd20, "contrast"), tip_data,
    attr(pd20, "weight"), attr(pd20, "levels")
  )
  expect_equal(result$score, fitch_check)
})

test_that("Multiple datasets produce verified scores", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  for (i in 1:3) {
    pd <- congreveLamsdellMatrices[[i]]
    tip_data <- t(vapply(pd, I, pd[[1]]))

    set.seed(3000 + i)
    result <- ts_random_wagner_tree(
      attr(pd, "contrast"), tip_data,
      attr(pd, "weight"), attr(pd, "levels")
    )

    fitch_check <- ts_fitch_score(
      result$edge,
      attr(pd, "contrast"), tip_data,
      attr(pd, "weight"), attr(pd, "levels")
    )
    expect_equal(result$score, fitch_check,
                 info = paste("Dataset", i))
  }
})
