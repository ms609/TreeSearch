source("split-support/config.R")

DataFile <- function(aln) {
  paste0("split-support/alignments/", aln, ".nex")
}

MBFile <- function(aln, suffix = NULL) {
  paste0("split-support/MrBayes/", aln, if(!is.null(suffix)) ".", suffix)
}

for (aln in alns) {
  
  partitions <- as.Splits(read.table(MBFile(aln, "parts"),
                                     skip = 2 + nTip)[, 2])
  pp <- read.table(MBFile(aln, "tstat"), skip = 1,
                   header = TRUE, comment.char = "")
  dataset <- MatrixToPhyDat(matrix(unlist(read.nexus.data(DataFile(aln))), 
                                   nrow = nTip, byrow = TRUE,
                                   dimnames = list(tips, NULL)))
  QuartetConcordance(partitions, dataset)
  
  if (file.exists(paste0(mbFile, ".trprobs"))) {
    message("Tree probabilities found for alignment ", aln)
  } else {
    on.exit(unlink(mbFile))
    writeLines(
      c(readLines(paste0("split-support/alignments/", aln, ".nex")), template),
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
