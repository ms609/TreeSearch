test_that("Character viewing and tree manipulation works", {
  app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "ViewChars"
  )
  on.exit(app$stop(), add = TRUE)

  # savePlotZip download content is snapshotted at the same states the legacy
  # tests/shinytest/ViewChars.R did (via snapshotDownload), so retiring the v1
  # script is lossless. normalize_download() (setup.R) scrubs volatile lines;
  # compare_file_text gives a readable diff on mismatch.
  expect_plot_script <- function() {
    app$expect_download("dl-savePlotZip",
                        compare = compare_file_text,
                        transform = normalize_download)
  }

  app$set_inputs(`data-dataSource` = "Agnarsson2004")
  app$wait_for_idle(timeout = 30000)
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-consP` = 0.5)
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-neverDrop` = "Argiope") # Avoid resetting root
  app$set_inputs(`consensus-keepNTips` = 1)
  app$set_inputs(`consensus-keepNTips` = 2)
  app$expect_values() # Check for correct display of invalid input (no script)

  app$set_inputs(`consensus-keepNTips` = 61)
  app$set_inputs(`consensus-outgroup` = "Argiope", timeout_ = 5000)
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-keepNTips` = 59) # Check tips kept legend changes to 17
  app$set_inputs(`consensus-excludedTip` = "Emertonella", timeout_ = 200000)
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-neverDrop` = "Emertonella")
  # QuickRogue triggered; keepNTips will change to 61
  app$set_inputs(`consensus-keepNTips` = 59)
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-outgroup` = character(0))
  app$set_inputs(`consensus-outgroup` = "Thymoites")
  app$expect_values()
  expect_plot_script()

  app$set_inputs(plotFormat = "ind")
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-plottedChar` = 0)
  app$set_inputs(`consensus-mapDisplay` = "tipsRight")
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-mapDisplay` = character(0))
  app$wait_for_idle(timeout = 5000)
  app$set_inputs(`consensus-plottedChar` = 1)
  app$set_inputs(`consensus-plottedChar` = 2)
  app$set_inputs(`consensus-plottedChar` = 3)
  app$set_inputs(`consensus-plottedChar` = 6)
  app$set_inputs(`consensus-plottedChar` = 7)
  app$set_inputs(`consensus-plottedChar` = 8)
  app$set_inputs(`consensus-plottedChar` = 7)
  app$wait_for_idle(timeout = 5000)
  app$expect_values()
  expect_plot_script()

  app$set_inputs(`consensus-plottedChar` = 8)
  app$set_inputs(`consensus-plottedChar` = 11)
  app$set_inputs(`consensus-plottedChar` = 53)
  app$wait_for_idle(timeout = 5000)
  app$set_inputs(`consensus-whichTree` = 7)
  app$wait_for_idle(timeout = 5000)
  app$expect_values()
  expect_plot_script()
})
