source("split-support/config.R")

if(!dir.exists("split-support/alignments")) {
  source("split-support/simulate.R")
}

for (aln in alns) {
  if (file.exists(paste0("split-support/tnt/", aln, ".ew.log"))) {
    message("Results found for ", aln)
  } else {
    system2(
      tntExec, 
      paste0("run split-support/tnt-ew.run ",
             "split-support/alignments/", aln, ".nex ",
             "split-support/TNT/", aln, ".ew.log ;")
    )
  }
}

warnings()
