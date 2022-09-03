source("split-support/config.R")

if(!dir.exists("split-support/alignments")) {
  source("split-support/simulate.R")
}

on.exit(unlink("*.tmp.tre"))
for (aln in alns) {
  if (file.exists(TNTFile(aln, "ew"))) {
    message("Results found for ", aln)
  } else {
    system2(
      tntExec, 
      paste0("run split-support/tnt-ew.run", " ",
             DataFile(aln), " ",
             TNTFile(aln, "ew"), " ;")
    )
  }
}

# Validation step
for (aln in alns) {
  thisFile <- TNTFile(aln, "ew")
  if (file.exists(thisFile) && length(readLines(thisFile)) < 60) {
    message("Deleting incomplete analysis: ", aln)
    file.remove(thisFile)
  }
}
