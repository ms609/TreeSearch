app <- ShinyDriver$new("../../", seed = 0, loadTimeout = 2e+05,
                       shinyOptions = list(test.mode = TRUE))
app$snapshotInit("SearchLog")

app$setInputs(dataSource = "Wills2012", timeout_ = 4000)
app$setInputs(searchConfig = "click")
app$setInputs(concavity = 1.1) # Set whilst visible; remembered later?
app$setInputs(epsilon = 1) # Set whilst visible; remembered later?
app$setInputs(`implied.weights` = "off")
app$setInputs(finalIter = 1.4)
app$setInputs(maxHits = 1.3)
app$setInputs(startIter = 1.6)
app$setInputs(ratchIter = 4)
app$setInputs(tbrIter = 2)
app$setInputs(modalGo = "click", timeout_ = 1e05)
app$setInputs(searchConfig = "click")
app$snapshotDownload("saveZip")
app$snapshotDownload("saveNwk")
app$setInputs(`implied.weights` = "on")
app$setInputs(ratchIter = 2)
app$setInputs(maxHits = 1)
app$setInputs(tbrIter = 1)
app$setInputs(startIter = 1.2)
app$setInputs(epsilon = 0) # No tolerance line here
app$setInputs(finalIter = 1)
app$setInputs(modalGo = "click", timeout_ = 2e05)
app$snapshotDownload("saveZip")
app$snapshotDownload("saveNex")

# Replace non-stable elements of file content with tags
downloads <- list.files("SearchLog-current", "*.download", full.names = TRUE)

for (file in downloads) {
  lines <- readLines(file)
  lines <- sub("TreeSearch session log: 2[\\d\\-]{9} [012][\\d:\\.]{7,}",
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
