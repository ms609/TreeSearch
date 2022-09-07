source("split-support/config.R")

if(!dir.exists("split-support/alignments")) {
  source("split-support/simulate.R")
}

for (aln in alns) {
  if (file.exists(IQFile(aln))) {
    message("Results found for ", aln)
  } else {
    seqs <- read.nexus.data(DataFile(aln))
    phyle <- IQFile(aln)
    on.exit(unlink(phyle))
    writeLines(c(
      paste(length(seqs), length(seqs[[1]])),
      "",
      apply(cbind(names(seqs),
            vapply(seqs, paste0, character(1), collapse = "")), 1, paste,
            collapse = "\t")
      ), phyle)
    system2(
      iqExec,
      # http://www.iqtree.org/doc/Command-Reference
      paste0(" -s ", phyle,
             " -st DNA ", # Sequence type: DNA
             " -mset JC ", # Model: Jukes-Cantor
             " -mrate E ", # Equal rates only
             " -nt auto -ntmax 6 ", # Number of threads
             " -seed 1 ", # Set random seed for reproducibility
             #" -b 1000", # Nonparametric bootstrap
             " -bb 1000 ", # Number of ultrafast bootstrap replicates
             " -bnni ", # Avoids branch support overestimates in UFB
             " -lbp 1000 ", # Fast local bootstrap probability
             " -alrt 1000 ", # Approximate likelihood ratio test
             " -abayes ", # Approximate Bayes test
             " --redo-tree ", # Overwrite previous run results
             ""
             )
    )
    
    # Remove unneeded results files
    keepExt <- c(
      "contree", # Consensus tree - why not
      "splits\\.nex" # Convergence diagnostics
    )
    
    outFiles <- list.files(path = iqDir, pattern = aln, full.names = TRUE)
    
    unlink(outFiles[-grep(paste0("(", paste0(keepExt, collapse = "|"), ")$"),
                          outFiles)])
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

