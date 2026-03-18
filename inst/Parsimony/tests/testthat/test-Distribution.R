test_that("Distribution and concordance plots render correctly", {
  app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "Distribution"
  )
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(`data-dataSource` = "Sun2018")
  app$set_inputs(plotFormat = "clus")
  app$set_inputs(`data-treeRange` = c(77, 125))
  app$expect_values()
  zipFile <- app$get_download("dl-savePlotZip")
  expect_true(file.exists(zipFile))

  app$set_inputs(`data-nTree` = 125)
  app$set_inputs(`data-treeRange` = c(1, 125))
  app$expect_values()
  zipFile <- app$get_download("dl-savePlotZip")
  expect_true(file.exists(zipFile))

  app$set_inputs(`clustering-clThresh` = 1)
  app$set_inputs(consP = 0.5)
  app$expect_values()

  app$set_inputs(concordance = "qc")
  app$expect_values()

  app$set_inputs(concordance = "clc", timeout_ = 6000)
  app$set_inputs(plotFormat = "ind", timeout_ = 6000)
  app$expect_values()

  app$set_inputs(plotFormat = "space")
  app$expect_values()

  app$set_inputs(`clustering-clThresh` = 0.5)
  app$expect_values()

  app$set_inputs(plotSize = 400)
  app$set_inputs(`treespace-spaceCol` = "score")
  app$set_inputs(`treespace-spaceDim` = 3)
  app$set_inputs(plotFormat = "space")
  app$set_inputs(`treespace-mapLines` = "seq")
  app$set_inputs(`treespace-relators` = c("Wiwaxia_corrugata", "Tonicella",
                               "Dentalium", "Phoronis"))
  app$expect_values()

  app$set_inputs(`treespace-mapLines` = character(0))
  app$set_inputs(`treespace-spaceCol` = "firstHit")
  app$set_inputs(`treespace-spaceCol` = "score")
  app$set_inputs(distMeth = "pid")
  app$expect_values()

  app$set_inputs(distMeth = "rf")
  app$expect_values()

  app$set_inputs(`data-dataSource` = "Agnarsson2004")
  app$wait_for_idle(timeout = 10000)
  app$set_inputs(distMeth = "qd")
  app$wait_for_idle(timeout = 5000)
  app$expect_values()
})
