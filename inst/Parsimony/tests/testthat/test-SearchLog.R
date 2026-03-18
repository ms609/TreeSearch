test_that("Search log workflow produces expected outputs", {
  app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "SearchLog"
  )
  on.exit(app$stop(), add = TRUE)

  # --- EW search ---
  app$set_inputs(`data-dataSource` = "Wills2012", timeout_ = 4000)
  app$click("search-searchConfig")
  app$wait_for_idle(timeout = 5000)
  app$set_inputs(`search-concavity` = 1.1)
  app$set_inputs(`search-epsilon` = 1)
  app$set_inputs(`search-implied.weights` = "off")
  app$set_inputs(`search-strategy` = "sprint")
  app$set_inputs(`search-maxReplicates` = 5)
  app$set_inputs(`search-targetHits` = 3)
  app$click("search-modalGo")
  # ExtendedTask returns immediately; poll exported counter
  app$wait_for_value(export = "searchCount",
                     ignore = list(NULL, 0L),
                     timeout = 120000)
  app$click("search-searchConfig")
  app$wait_for_idle(timeout = 5000)

  zipFile <- app$get_download("dl-saveZip")
  expect_true(file.exists(zipFile))
  nwkFile <- app$get_download("dl-saveNwk")
  expect_true(file.exists(nwkFile))

  # --- IW search ---
  app$set_inputs(`search-implied.weights` = "on")
  app$set_inputs(`search-strategy` = "default")
  app$set_inputs(`search-maxReplicates` = 3)
  app$set_inputs(`search-targetHits` = 2)
  app$set_inputs(`search-epsilon` = 0)
  app$click("search-modalGo")
  app$wait_for_value(export = "searchCount",
                     ignore = list(NULL, 0L, 1L),
                     timeout = 200000)

  zipFile2 <- app$get_download("dl-saveZip")
  expect_true(file.exists(zipFile2))
  nexFile <- app$get_download("dl-saveNex")
  expect_true(file.exists(nexFile))
})
