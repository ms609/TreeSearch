test_that("App starts and default dataset loads without errors", {
  app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "Smoke"
  )
  on.exit(app$stop(), add = TRUE)

  # Default dataset (Wills2012) should auto-load
  app$wait_for_idle(timeout = 10000)

  # Verify app is alive and has trees loaded
  vals <- app$get_values()
  expect_true(!is.null(vals$input$dataSource))

  # Sidebar should show tree count
  results <- vals$output$results
  expect_true(is.character(results))
  expect_true(grepl("trees? in memory", results))

  # No error notifications (shiny-notification-panel hidden in test mode,

  # but we can verify the search count export is initialized)
  expect_equal(vals$export$searchCount, 0L)
})
