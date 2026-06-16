# The native scorer (PrepareData + TreeScore/EdgeListScore) computes the CORRECT
# Brazeau-Guillerme-Smith inapplicable scores -- NOT bit-equivalence with legacy
# MorphyLib, which mis-scored `{1-}` ambiguous-with-inapplicable tokens
# (Lobo pattern 93).

test_that("PrepareData() builds a ParsimonyData handle", {
  data("Lobo", package = "TreeTools")
  obj <- PrepareData(Lobo.phy)
  expect_true(is.ParsimonyData(obj))
  expect_equal(obj[["nTip"]], length(Lobo.phy))
})

test_that("native TreeScore() equals native TreeLength()", {
  data("Lobo", package = "TreeTools")
  data("inapplicable.datasets")
  datasets <- list(Lobo = Lobo.phy,
                   CL42 = congreveLamsdellMatrices[[42]],
                   CL10 = congreveLamsdellMatrices[[10]])
  for (ds in datasets) {
    obj <- PrepareData(ds)
    for (tr in list(TreeTools::PectinateTree(names(ds)),
                    TreeTools::BalancedTree(names(ds)))) {
      expect_equal(TreeScore(tr, obj), TreeLength(tr, ds))
    }
  }
})

test_that("native scorer fixes the MorphyLib {1-} bug (Lobo pattern 93)", {
  data("Lobo", package = "TreeTools")
  tr <- TreeTools::PectinateTree(names(Lobo.phy))
  obj <- PrepareData(Lobo.phy)
  # Correct BGS length is 273; the legacy MorphyLib path returned 274.
  expect_equal(TreeScore(tr, obj), 273)
})

test_that("SingleCharData() + EdgeListScore() score natively per character", {
  data("Lobo", package = "TreeTools")
  ds <- Lobo.phy
  tr <- TreeTools::RenumberTips(TreeTools::PectinateTree(names(ds)), names(ds))
  tr <- TreeTools::Postorder(tr)
  strs <- TreeTools::PhyToString(ds, byTaxon = FALSE, useIndex = FALSE,
                                 concatenate = FALSE)
  objs <- lapply(strs, SingleCharData)
  e <- tr[["edge"]]
  perChar <- vapply(objs, function(o) EdgeListScore(e[, 1], e[, 2], o),
                    double(1))
  expect_equal(sum(perChar * attr(ds, "weight")), TreeLength(tr, ds))
})
