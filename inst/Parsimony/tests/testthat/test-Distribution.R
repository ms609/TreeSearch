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
  # Let the dataset load finish before setting anything that depends on it.
  # UpdateAllTrees() resets treeRange to the full 1:nTrees span whenever the tree
  # count changes (mod_data.R), so a `data-treeRange` that lands mid-load is
  # silently clobbered by the reset -- and set_inputs()' default 4 s wait is not
  # long enough for Sun2018's 125 trees, as the "Server did not update any output
  # values within 4 seconds" warning on this line reported. Whether the clobber
  # beat the capture varied by machine, which is what made this snapshot record
  # trees[1:125] on some runs and trees[77:125] on others.
  wait_stable(app)
  app$set_inputs(plotFormat = "clus")
  app$set_inputs(`data-treeRange` = c(77, 125))
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(`data-nTree` = 125)
  app$set_inputs(`data-treeRange` = c(1, 125))
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(`clustering-clThresh` = 1)
  app$set_inputs(`consensus-consP` = 0.5)
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(`consensus-concordance` = "qc")
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(`consensus-concordance` = "clc", timeout_ = 6000)
  app$set_inputs(plotFormat = "ind", timeout_ = 6000)
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(plotFormat = "space")
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(`clustering-clThresh` = 0.5)
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(plotSize = 400)
  app$set_inputs(`treespace-spaceCol` = "score")
  app$set_inputs(`treespace-spaceDim` = 3)
  app$set_inputs(plotFormat = "space")
  app$set_inputs(mapLines = "seq")
  app$set_inputs(`treespace-relators` = c("Wiwaxia_corrugata", "Tonicella",
                               "Dentalium", "Phoronis"))
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(mapLines = character(0))
  app$set_inputs(`treespace-spaceCol` = "firstHit")
  app$set_inputs(`treespace-spaceCol` = "score")
  app$set_inputs(distMeth = "pid")
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(distMeth = "rf")
  wait_stable(app)
  expect_plot_script()

  app$set_inputs(`data-dataSource` = "Agnarsson2004")
  wait_stable(app, 10000)
  app$set_inputs(distMeth = "qd")
  wait_stable(app, 5000)
  wait_stable(app)
  expect_plot_script()
})
