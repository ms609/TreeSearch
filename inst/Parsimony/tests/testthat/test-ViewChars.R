test_that("Character viewing and tree manipulation works", {
  app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "ViewChars"
  )
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(dataSource = "Agnarsson2004")
  app$expect_values()
  zipFile <- app$get_download("savePlotZip")
  expect_true(file.exists(zipFile))

  app$set_inputs(consP = 0.5)
  app$expect_values()

  app$set_inputs(neverDrop = "Argiope") # Avoid resetting root
  app$set_inputs(keepNTips = 1)
  app$set_inputs(keepNTips = 2)
  app$expect_values() # Check for correct display of invalid input

  app$set_inputs(keepNTips = 61)
  app$set_inputs(outgroup = "Argiope", timeout_ = 5000)
  app$expect_values()

  app$set_inputs(keepNTips = 59) # Check tips kept legend changes to 17
  app$set_inputs(excludedTip = "Emertonella", timeout_ = 200000)
  app$expect_values()

  app$set_inputs(neverDrop = "Emertonella")
  # QuickRogue triggered; keepNTips will change to 61
  app$set_inputs(keepNTips = 59)
  app$expect_values()

  app$set_inputs(outgroup = character(0))
  app$set_inputs(outgroup = "Thymoites")
  app$expect_values()

  app$set_inputs(plotFormat = "ind")
  app$expect_values()

  app$set_inputs(plottedChar = 0)
  app$set_inputs(mapDisplay = "tipsRight")
  app$expect_values()

  app$set_inputs(mapDisplay = character(0))
  app$set_inputs(plottedChar = 1)
  app$set_inputs(plottedChar = 2)
  app$set_inputs(plottedChar = 3)
  app$set_inputs(plottedChar = 6)
  app$set_inputs(plottedChar = 7)
  app$set_inputs(plottedChar = 8)
  app$set_inputs(plottedChar = 7)
  app$expect_values()

  app$set_inputs(plottedChar = 8)
  app$set_inputs(plottedChar = 11)
  app$set_inputs(plottedChar = 53)
  app$set_inputs(whichTree = 7)
  app$expect_values()
})
