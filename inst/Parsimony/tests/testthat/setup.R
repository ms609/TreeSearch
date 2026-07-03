library(shinytest2)

# ---------------------------------------------------------------------------
# new_app_driver(): shared AppDriver factory for the three integration tests.
#
#  * Disable expect_values() screenshots (expect_values_screenshot_args = FALSE)
#    everywhere. The browser screenshot (.png) is flaky even run-to-run on one
#    machine (anti-aliasing / render timing), so it is not a reliable signal;
#    the value (json) + download-content snapshots carry the real coverage.
#  * Generous load_timeout for the headless-Chrome boot under load.
# ---------------------------------------------------------------------------
new_app_driver <- function(name, ...) {
  AppDriver$new(
    app_dir = "../../",
    seed = 0,
    load_timeout = 200000,
    shiny_args = list(test.mode = TRUE),
    name = name,
    expect_values_screenshot_args = FALSE,
    ...
  )
}

# ---------------------------------------------------------------------------
# expect_vals(): expect_values() wrapper.
#
# expect_values() records the app's rendered plot output as an OS-dependent
# hashed data URI, so the value snapshots are only portable when the runner
# matches the machine that generated the baselines. The shiny CI job therefore
# runs on Windows (the development platform) -- see .github/workflows -- which
# keeps the rendered output consistent with the committed baselines. Kept as a
# thin wrapper so a normalising transform can be slotted in here if a residual
# rendering difference ever surfaces.
# ---------------------------------------------------------------------------
expect_vals <- function(app) {
  app$expect_values()
}

# ---------------------------------------------------------------------------
# wait_stable(): retrying wrapper around app$wait_for_idle().
#
# wait_for_idle() can transiently error ("An error occurred while waiting for
# Shiny to be stable"), most often right after a heavy dataset load / plot.
# Retry a couple of times before propagating, so a transient chromote hiccup
# doesn't fail an otherwise-healthy run. A first-try success is the common path
# and incurs no delay.
# ---------------------------------------------------------------------------
wait_stable <- function(app, timeout = 30000, attempts = 3L) {
  for (i in seq_len(attempts)) {
    ok <- tryCatch(
      { app$wait_for_idle(timeout = timeout); TRUE },
      error = function(e) if (i == attempts) stop(e) else FALSE
    )
    if (isTRUE(ok)) break
    Sys.sleep(1)
  }
  invisible(app)
}

# ---------------------------------------------------------------------------
# normalize_download(): transform for AppDriver$expect_download().
#
# Scrubs volatile content from downloaded R-script / Nexus / log files so the
# snapshot baseline is portable across machines and runs. Ports the ad-hoc
# normalization the legacy `tests/shinytest/*.R` scripts did by hand.
#
# testthat's expect_snapshot_file() calls transform() with a character vector
# of the file's lines. Applied to BOTH the current output and the stored
# baseline, so the substitutions only need to be self-consistent.
#
# The leading sub("\r$", ...) is load-bearing: authored on Windows but the
# sense-check leg and other contributors run Linux, where readLines() keeps a
# trailing \r on CRLF files. Stripping it here makes the baseline byte-identical
# regardless of the OS that generated it (a prior CRLF mismatch corrupted a
# Hamilton deploy; see memory hamilton-deploy-crlf-corruption).
# ---------------------------------------------------------------------------
normalize_download <- function(lines) {
  lines <- sub("\r$", "", lines)
  # Session / plot log timestamps: "# # TreeSearch session log: 2026-.. # # #"
  lines <- sub("(TreeSearch (session|plot) log:).*",
               "\\1 <DATE-AND-TIME>", lines, perl = TRUE)
  # ape's write.nexus header: "[R-package APE, Fri Jul 03 ...]"
  lines <- sub("\\[R-package APE, .*\\]",
               "[R-package APE, <DATE-AND-TIME>]", lines, perl = TRUE)
  # System line (saveZip log; savePlotZip already scrubs its own in testmode)
  lines <- sub("^# System: .*", "# System: <SYS-INFO>", lines, perl = TRUE)
  # Package/R version lines: "# - TreeSearch 2.0.0.9999" -> "# - TreeSearch <VERSION>"
  lines <- sub("^(# - [A-Za-z.]+ ).*", "\\1<VERSION>", lines, perl = TRUE)
  lines
}
