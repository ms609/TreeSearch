# Extracted from test-app-smoke.R:8

# test -------------------------------------------------------------------------
app <- AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = "Smoke"
  )
