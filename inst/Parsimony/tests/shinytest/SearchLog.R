app <- ShinyDriver$new("../../", seed = 0, loadTimeout = 2e+05,
                       shinyOptions = list(test.mode = TRUE))
app$snapshotInit("SearchLog")

# Helper: poll exported searchCount until it exceeds `prev`.
# Needed because MaximizeParsimony runs asynchronously via ExtendedTask;
# setInputs(modalGo = "click") returns immediately after invoke().
waitForSearch <- function(app, prev = 0L, timeout_s = 120, poll_s = 2) {
  deadline <- Sys.time() + timeout_s
  while (Sys.time() < deadline) {
    vals <- app$getAllValues()
    count <- vals$export$searchCount
    if (!is.null(count) && count > prev) return(invisible(count))
    Sys.sleep(poll_s)
  }
  stop("Timed out waiting for search to complete")
}

app$setInputs(dataSource = "Wills2012", timeout_ = 4000)
app$setInputs(searchConfig = "click")
app$setInputs(concavity = 1.1) # Set whilst visible; remembered later?
app$setInputs(epsilon = 1) # Set whilst visible; remembered later?
app$setInputs(`implied.weights` = "off")
app$setInputs(strategy = "sprint")
app$setInputs(maxReplicates = 5)
app$setInputs(targetHits = 3)
app$setInputs(modalGo = "click")
searchesDone <- waitForSearch(app, prev = 0L)
app$setInputs(searchConfig = "click")
app$snapshotDownload("saveZip")
app$snapshotDownload("saveNwk")
app$setInputs(`implied.weights` = "on")
app$setInputs(strategy = "default")
app$setInputs(maxReplicates = 3)
app$setInputs(targetHits = 2)
app$setInputs(epsilon = 0) # No tolerance line here
app$setInputs(modalGo = "click")
searchesDone <- waitForSearch(app, prev = searchesDone)
app$snapshotDownload("saveZip")
app$snapshotDownload("saveNex")

# Replace non-stable elements of file content with tags
downloads <- list.files("SearchLog-current", "*.download", full.names = TRUE)

for (file in downloads) {
  lines <- readLines(file)
  lines <- sub("TreeSearch session log: 2[\\d\\-]{9} [012][\\d:]{7}",
               "TreeSearch session log: <DATE-AND-TIME>", 
               lines, perl = TRUE)
  lines[2] <- sub("\\[R-package APE, .*\\]",
               "[R-package APE, <DATE-AND-TIME>]",
               lines[2], perl = TRUE)
  lines[4] <- sub("# System: .*", "# System: <SYS-INFO>", lines[4])
  lines[5:9] <- sub("^(# \\- \\w+ ).*$", "\\1<VERSION>",
                    lines[5:9], perl = TRUE)
  writeLines(lines, con = file)
}
