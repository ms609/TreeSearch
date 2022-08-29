# Size of tree and dataset
nTip <- 48
tips <- paste0("t", seq_len(nTip))
nChar <- nTip * 2

# Number of alignments to analyse
nAln <- 1000
alns <- paste0("aln", formatC(seq_len(nAln), width = 4, flag = 0))


# Download command-driven TNT from 
# http://www.lillo.org.ar/phylogeny/tnt/ZIPCHTNT.ZIP
# and unzip to a convenient local folder.

# Command to launch TNT executable
tntExec <- "C:/Programs/Phylogeny/tnt/tnt.exe"

# Download MrBayes from 
# https://nbisweden.github.io/MrBayes/download.html
# and unzip/install to a convenient local folder.

# Command to launch MrBayes executable
mbPath <- "C:/Programs/Phylogeny/MrBayes/bin/mb.3.2.7-win64.exe"


# Patterns to use when creating files
ConcFile <- function(aln) {
  paste0("split-support/concordance/", aln, ".txt")
}

DataFile <- function(aln) {
  paste0("split-support/alignments/", aln, ".nex")
}

MBFile <- function(aln, suffix = NULL) {
  paste0("split-support/MrBayes/", aln, if(!is.null(suffix)) ".", suffix)
}

TNTFile <- function(aln, wt) {
  paste0("split-support/TNT/", aln, ".", wt, ".sym")
}