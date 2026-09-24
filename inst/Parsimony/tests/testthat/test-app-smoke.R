test_that("App starts and a dataset loads without errors", {
  app <- new_app_driver("Smoke")
  on.exit(app$stop(), add = TRUE)

  # The app boots with no dataset (data-dataSource defaults to "file", the
  # upload placeholder), so nothing is "in memory" until a dataset is chosen.
  # Load a bundled tree distribution and confirm the app ingests it.
  app$set_inputs(`data-dataSource` = "Sun2018")
  wait_stable(app, 30000)

  vals <- app$get_values()
  # dataSource is namespaced in the data module
  expect_equal(vals$input[["data-dataSource"]], "Sun2018")

  # Sidebar should now report the loaded tree count ("N trees in memory: ...")
  results_text <- paste(as.character(vals$output[["search-results"]]),
                        collapse = "")
  expect_true(grepl("trees? in memory", results_text))

  # No search has been run yet
  expect_equal(vals$export$searchCount, 0L)
})
