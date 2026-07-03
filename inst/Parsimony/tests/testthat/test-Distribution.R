test_that("Distribution and concordance plots render correctly", {
  app <- new_app_driver("Distribution")
  on.exit(app$stop(), add = TRUE)

  # savePlotZip download content is snapshotted at the same states the legacy
  # tests/shinytest/Distribution.R did (via snapshotDownload), so retiring the
  # v1 script is lossless. normalize_download() (setup.R) scrubs volatile lines.
  expect_plot_script <- function() {
    app$expect_download("dl-savePlotZip",
                        compare = compare_file_text,
                        transform = normalize_download)
  }

  app$set_inputs(`data-dataSource` = "Sun2018")
  app$set_inputs(plotFormat = "clus")
  app$set_inputs(`data-treeRange` = c(77, 125))
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(`data-nTree` = 125)
  app$set_inputs(`data-treeRange` = c(1, 125))
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(`clustering-clThresh` = 1)
  app$set_inputs(`consensus-consP` = 0.5)
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(`consensus-concordance` = "qc")
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(`consensus-concordance` = "clc", timeout_ = 6000)
  app$set_inputs(plotFormat = "ind", timeout_ = 6000)
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(plotFormat = "space")
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(`clustering-clThresh` = 0.5)
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(plotSize = 400)
  app$set_inputs(`treespace-spaceCol` = "score")
  app$set_inputs(`treespace-spaceDim` = 3)
  app$set_inputs(plotFormat = "space")
  app$set_inputs(mapLines = "seq")
  app$set_inputs(`treespace-relators` = c("Wiwaxia_corrugata", "Tonicella",
                               "Dentalium", "Phoronis"))
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(mapLines = character(0))
  app$set_inputs(`treespace-spaceCol` = "firstHit")
  app$set_inputs(`treespace-spaceCol` = "score")
  app$set_inputs(distMeth = "pid")
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(distMeth = "rf")
  expect_vals(app)
  expect_plot_script()

  app$set_inputs(`data-dataSource` = "Agnarsson2004")
  wait_stable(app, 10000)
  app$set_inputs(distMeth = "qd")
  wait_stable(app, 5000)
  expect_vals(app)
  expect_plot_script()
})
