test_that("PrepareDataProfile() handles empty matrices", {
  dat <- TreeTools::MatrixToPhyDat(matrix(c(0, 1, rep("?", 5)),
                                          dimnames = list(letters[1:7], NULL)))
  result <- suppressMessages(PrepareDataProfile(dat))
  expect_equal(attr(result, "info.amounts"), numeric(0))
  expect_equal(attr(result, "nr"), 0L)
})

Dehash <- function (x) {
  lapply(x, function (xi) {
    attr(xi, ".match.hash") <- NULL
    if (!is.null(dimnames(xi))) {
      dimnames(xi) <- Dehash(dimnames(xi))
    }
    xi
  })
}

test_that("PrepareDataProfile()", {
  
  # Easy one
  mtx <- cbind(c("0", "0", 1,1,1,1),
               c(0,0,1,1,1,1),# again
               c(0,0,0,1,1,"?"))
  rownames(mtx) <- letters[seq_len(nrow(mtx))]
  phy1 <- TreeTools::MatrixToPhyDat(mtx)
  expect_equivalent(phy1, PrepareDataProfile(phy1))
  # PrepareDataProfile renormalizes token labels to 1..k; check structural
  # attributes that should be preserved, not levels/allLevels/contrast
  pp1 <- PrepareDataProfile(phy1)
  expect_equal(attr(pp1, "weight"), attr(phy1, "weight"))
  expect_equal(attr(pp1, "nr"), attr(phy1, "nr"))
  expect_equal(attr(pp1, "nc"), attr(phy1, "nc"))
  expect_equal(attr(pp1, "index"), attr(phy1, "index"))
  
  # Flipped binary char: PrepareDataProfile does not flip-normalize, so phy2
  # produces 3 unique patterns (not deduplicated with phy1's 2)
  mtx <- cbind(c("0", "0", 1,1,1,1),
               c(1,1,0,0,0,0),# flipped
               c(0,0,0,1,1,"{012}"))
  rownames(mtx) <- letters[seq_len(nrow(mtx))]
  phy2 <- TreeTools::MatrixToPhyDat(mtx)
  pp2 <- PrepareDataProfile(phy2)
  expect_equal(attr(pp2, "nr"), 3L)
  expect_equal(attr(pp2, "nc"), attr(pp1, "nc"))
  # Both informative binary patterns have the same information content
  expect_equal(attr(pp2, "info.amounts")[, 1], attr(pp1, "info.amounts")[, 1])
  
  
  mtx <- cbind(c("0", "0", 1,1,1, "2", "2", 3,3,3,3),
               c("?", "?", 1,1,1, "?", "?", 0,0,0,0),
               c(0,0,1,1,1,2,2,3,3,3,3),# again
               c(rep("?", 5), "2", "2", 0,0,0,0),
               c("?", "?", 1,1,1, 1,1, 0,0,0,0),
               c("0", "1", rep("?", 9))
               )
  rownames(mtx) <- letters[seq_len(nrow(mtx))]
  dataset <- TreeTools::MatrixToPhyDat(mtx)
  
  # After T-107/T-144: 4-state chars with 11 tips are within the
  # MaddisonSlatkin feasibility threshold (k=4, max=18 tips), so they are
  # preserved as 4-state — no binary reduction, no warning.
  pd <- PrepareDataProfile(dataset)
  expect_equal(4L, length(attr(pd, "levels")))
  expect_equal(dim(PhyDatToMatrix(pd)), c(11L, 6L))
  expect_equal(c(1L, 2L, 1L, 3L, 4L, 5L), attr(pd, "index"))
  expect_equal(c(2L, 1L, 1L, 1L, 1L), attr(pd, "weight"))
  
  dataset2 <- TreeTools::MatrixToPhyDat(mtx[!mtx[, 1] %in% c(0, 2), ])
  # The first informative pattern in dataset2 matches the informative pattern
  # in pd (both are the same 2-state split with 3 tips vs 4 tips)
  expect_equal(attr(PrepareDataProfile(dataset2), "info.amounts")[, 1, drop = FALSE],
               attr(pd, "info.amounts")[1:3, 2, drop = FALSE])
  
  
  data("Lobo", package = "TreeTools")
  expect_warning(prep <- PrepareDataProfile(Lobo.phy))
  expect_equal(c(17, attr(prep, "nr")),
               dim(attr(prep, "info.amounts")))
  
})
