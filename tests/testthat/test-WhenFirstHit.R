test_that("WhenFirstHit()", {
  library("TreeTools", quietly = TRUE)
  trees <- list(
     seed_00 = as.phylo(1, 8),
     ratch1_01 = as.phylo(2, 8),
     ratch1_02 = as.phylo(3, 8),
     ratch4_44 = as.phylo(4, 8),
     final_99 = as.phylo(5, 8)
  )
  wfhTrees <- WhenFirstHit(trees)
  expect_equal(attr(wfhTrees, "firstHit"),
               structure(c(seed = 1L, ratch1 = 2L, ratch4 = 1L, final = 1L),
                         dim = 4L,
                         dimnames = list(
                           whenHit = c("seed", "ratch1", "ratch4", "final")
                         ), class = "table")
  )
  
  expect_equal(attr(WhenFirstHit(trees), "firstHit"),
               attr(WhenFirstHit(wfhTrees), "firstHit"))
  
  noInfo <- as.phylo(1:10, 8)
  expect_equal(WhenFirstHit(noInfo), noInfo)
})
