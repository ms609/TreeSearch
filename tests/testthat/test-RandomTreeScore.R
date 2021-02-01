test_that("Errors are handled", {
  mo <- mpl_new_Morphy()
  expect_warning(RandomTreeScore(1, mo))
  mpl_delete_Morphy(mo)
  
  expect_error(RandomMorphyTree(-1))
  expect_error(RandomMorphyTree(0))
  expect_error(RandomMorphyTree(1))
})
