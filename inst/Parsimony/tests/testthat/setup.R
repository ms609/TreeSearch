library(shinytest2)

#' Set options for the duration of `code`, restoring the prior values on
#' exit. Minimal single-purpose stand-in for withr::with_options so the
#' package doesn't need withr just for this.
with_options <- function(new, code) {
  old <- options(new)
  on.exit(options(old))
  code
}

# ---------------------------------------------------------------------------
# new_app_driver(): shared AppDriver factory for the three integration tests.
#
#  * Disable expect_values() screenshots (expect_values_screenshot_args = FALSE)
#    everywhere. The browser screenshot (.png) is flaky even run-to-run on one
#    machine (anti-aliasing / render timing), so it is not a reliable signal.
#  * Generous load_timeout for the headless-Chrome boot under load.
#
# NB: these tests deliberately do NOT use app$expect_values(). Its value
# snapshots embed the app's rendered plot output -- both a data-URI image hash
# AND float plot geometry (usr coordinates, bounding boxes) -- which differ from
# machine to machine even on the same OS (confirmed on GHA windows-latest vs the
# dev box). The portable coverage is instead: the download-script snapshots
# (expect_download, which encode the app state), the interaction path itself (a
# crash / binding error / thrown reactive fails the test), and SearchLog's
# deterministic log markers. wait_stable() replaces the implicit settle-wait
# that expect_values() used to provide before each download capture.
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
# wait_stable(): retrying wrapper around app$wait_for_idle().
#
# wait_for_idle() can transiently error ("An error occurred while waiting for
# Shiny to be stable"), most often right after a heavy dataset load / plot.
# Retry a couple of times before propagating, so a transient chromote hiccup
# doesn't fail an otherwise-healthy run. A first-try success is the common path
# and incurs no delay.
#
# debounceWait is load-bearing, not padding. mod_data's nTree / treeRange
# watchers are debounce()d (aJiffy = 42 ms, typingJiffy = 105 ms), and a pending
# debounce timer does NOT make Shiny busy -- there is nothing to recompute until
# it expires, so wait_for_idle() can return before a debounced watcher has even
# seen the input the test just set. A download captured at that point encodes the
# state BEFORE the last set_inputs(); on a machine where the timer does fire in
# time it encodes the state after. The same test then yields different snapshots
# on different machines. Sleeping past the longest debounce window and waiting
# again lets that work start and finish, which is what makes these baselines
# reproducible rather than timing-dependent.
#
# This is how the Distribution baseline came to record `trees[1:125]` for a
# state its test had set to c(77, 125) -- noticed only once the Coreset
# dependency fix let CI reach the suite at all.
# ---------------------------------------------------------------------------
wait_stable <- function(app, timeout = 30000, attempts = 3L,
                        debounceWait = 0.25) {
  for (i in seq_len(attempts)) {
    ok <- tryCatch(
      {
        app$wait_for_idle(timeout = timeout)
        Sys.sleep(debounceWait)
        app$wait_for_idle(timeout = timeout)
        TRUE
      },
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
  # WideSample() thinning seed (T-360): freshly drawn each run, so scrub the
  # value while still asserting that a set.seed() call precedes WideSample().
  lines <- sub("^set\\.seed\\([0-9]+\\)$", "set.seed(<SEED>)", lines, perl = TRUE)
  lines
}
