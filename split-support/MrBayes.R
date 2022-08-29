if(!dir.exists("split-support/alignments")) {
  source("simulate.R")
}

# Download MrBayes from 
# https://nbisweden.github.io/MrBayes/download.html
# and unzip/install to a convenient local folder.

mbPath <-  "C:/Programs/Phylogeny/MrBayes/bin/mb.3.2.7-win64.exe"

template <- readLines("split-support/mb.nex")

for (i in formatC(1:1000, width = 4, flag = 0)) {
  mbFile <- paste0("split-support/MrBayes/aln", i)
  
  if (file.exists(paste0(mbFile, ".trprobs"))) {
    message("Tree probabilities found for alignment ", i)
  } else {
    on.exit(unlink(mbFile))
    writeLines(
      c(readLines(paste0("split-support/alignments/aln", i, ".nex")), template),
      mbFile
    )
    system2(mbPath, mbFile)
    
    
    # Remove unneeded results files
    keepExt <- c(
      "con\\.tre", # Consensus tree - why not
      "parts", "tstat", # Partitions an dprobabilities
      # "trprobs" # Sampled trees and probabilities
      # "mcmc" # Standard deviations of splits - see tstat
      "pstat" # Convergence diagnostics
    )
    
    outFiles <- list.files(path = "split-support/MrBayes/",
                           pattern = paste0("aln", i),
                           full.names = TRUE)
    
    unlink(outFiles[-grep(paste0("(", paste0(keepExt, collapse = "|"), ")$"),
                          outFiles)])
    
  }
}
