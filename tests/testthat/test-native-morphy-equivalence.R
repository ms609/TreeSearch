# Stage 1 of dropping MorphyLib: MorphyLength() scores the default (inapplicable)
# gap mode via the native Fitch kernel. We assert the CORRECT native
# (Brazeau-Guillerme-Smith) scores -- NOT bit-equivalence with legacy MorphyLib,
# which mis-scores `{1-}` ambiguous-with-inapplicable tokens (Lobo pattern 93).

test_that("PhyDat2Morphy() attaches native data only for the inapplicable mode", {
  data("Lobo", package = "TreeTools")
  mo <- PhyDat2Morphy(Lobo.phy)                    # default gap = "inapplicable"
  on.exit(UnloadMorphy(mo), add = TRUE)
  expect_false(is.null(attr(mo, "native", exact = TRUE)))

  ma <- PhyDat2Morphy(Lobo.phy, "ambiguous")       # rare mode: MorphyLib for now
  on.exit(UnloadMorphy(ma), add = TRUE)
  expect_null(attr(ma, "native", exact = TRUE))
})

test_that("native MorphyTreeLength() equals native TreeLength()", {
  data("Lobo", package = "TreeTools")
  data("inapplicable.datasets")
  datasets <- list(Lobo = Lobo.phy,
                   CL42 = congreveLamsdellMatrices[[42]],
                   CL10 = congreveLamsdellMatrices[[10]])
  for (ds in datasets) {
    for (tr in list(TreeTools::PectinateTree(names(ds)),
                    TreeTools::BalancedTree(names(ds)))) {
      mo <- PhyDat2Morphy(ds)
      expect_equal(MorphyTreeLength(tr, mo), round(TreeLength(tr, ds)))
      UnloadMorphy(mo)
    }
  }
})

test_that("native MorphyLength fixes the MorphyLib {1-} bug (Lobo pattern 93)", {
  data("Lobo", package = "TreeTools")
  tr <- TreeTools::PectinateTree(names(Lobo.phy))
  mo <- PhyDat2Morphy(Lobo.phy)
  on.exit(UnloadMorphy(mo))
  # Correct BGS length is 273; the legacy MorphyLib path returned 274.
  expect_equal(MorphyTreeLength(tr, mo), 273)
})

test_that("SingleCharMorphy() + MorphyLength() score natively per character", {
  data("Lobo", package = "TreeTools")
  ds <- Lobo.phy
  tr <- TreeTools::PectinateTree(names(ds))
  strs <- TreeTools::PhyToString(ds, byTaxon = FALSE, useIndex = FALSE,
                                 concatenate = FALSE)
  mObjs <- lapply(strs, SingleCharMorphy)
  on.exit(lapply(mObjs, UnloadMorphy))
  e <- tr[["edge"]]
  perChar <- vapply(mObjs, MorphyLength, parent = e[, 1], child = e[, 2],
                    integer(1))
  expect_equal(sum(perChar * attr(ds, "weight")), round(TreeLength(tr, ds)))
})
