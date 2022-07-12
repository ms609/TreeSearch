test_that("QuartetResolution()", {
  expect_equal(
    QuartetResolution(inapplicable.trees[["Vinther2008"]],
                      c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")),
               c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                 3, 3, 3, 2, 1, 3, 1, 1, 3, 2, 1, 1, 1, 3, 1, 3, 3, 3, 1, 1, 3,
                 3, 2, 3, 1, 3, 2, 2, 3, 2, 1, 3, 1, 1, 3, 2, 2, 2, 1, 2, 3, 2,
                 3, 2, 1, 1, 3, 2, 3, 2, 3, 3, 1))
  expect_equal(
    QuartetResolution(inapplicable.trees[["Vinther2008"]],
                      c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")),
    QuartetResolution(inapplicable.trees[["Vinther2008"]],
                      c("Nemertean", "Halkieria", "Wiwaxia", "Acaenoplax"))
  )
               
})
