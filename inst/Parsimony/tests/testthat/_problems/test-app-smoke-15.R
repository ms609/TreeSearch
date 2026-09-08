# Extracted from test-app-smoke.R:15

# test -------------------------------------------------------------------------
app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "Smoke"
  )
on.exit(app$stop(), add = TRUE)
app$wait_for_idle(timeout = 10000)
vals <- app$get_values()
