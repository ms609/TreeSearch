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
  app$set_inputs(dataSource = "Wills2012", timeout_ = 4000)
  app$click("searchConfig")
  app$wait_for_idle(timeout = 5000)
  app$set_inputs(concavity = 1.1)
  app$set_inputs(epsilon = 1)
  app$set_inputs(`implied.weights` = "off")
  app$set_inputs(strategy = "sprint")
  app$set_inputs(maxReplicates = 5)
  app$set_inputs(targetHits = 3)
  app$click("modalGo")
  # ExtendedTask returns immediately; poll exported counter
  app$wait_for_value(export = "searchCount",
                     ignore = list(NULL, 0L),
                     timeout = 120000)
  app$click("searchConfig")
  app$wait_for_idle(timeout = 5000)

  zipFile <- app$get_download("saveZip")
  expect_true(file.exists(zipFile))
  nwkFile <- app$get_download("saveNwk")
  expect_true(file.exists(nwkFile))

  # --- IW search ---
  app$set_inputs(`implied.weights` = "on")
  app$set_inputs(strategy = "default")
  app$set_inputs(maxReplicates = 3)
  app$set_inputs(targetHits = 2)
  app$set_inputs(epsilon = 0)
  app$click("modalGo")
  app$wait_for_value(export = "searchCount",
                     ignore = list(NULL, 0L, 1L),
                     timeout = 200000)

  zipFile2 <- app$get_download("saveZip")
  expect_true(file.exists(zipFile2))
  nexFile <- app$get_download("saveNex")
  expect_true(file.exists(nexFile))
})
