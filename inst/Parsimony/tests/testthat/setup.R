library(shinytest2)

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
