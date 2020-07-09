context("RMorphy.C")

test_that("NULL pointers don't cause crash", {
  ptr <- mpl_new_Morphy()
  expect_equal(0, mpl_delete_Morphy(ptr))
  expect_true(is.na(mpl_delete_Morphy(ptr)))
})

test_that("Pointers survive garbage collection", {
  ptr <- mpl_new_Morphy()
  gc()
  expect_equal(0, mpl_delete_Morphy(ptr))
})
