### Analytical setup ###

# Size of tree and dataset
nTip <- 48
tips <- paste0("t", seq_len(nTip))
nChar <- nTip * 2

# Number of alignments to analyse
nAln <- 1000
alns <- paste0("aln", formatC(seq_len(nAln), width = 4, flag = 0))

### Location of executables ###

# Download command-driven TNT from 
# http://www.lillo.org.ar/phylogeny/tnt/ZIPCHTNT.ZIP
# and unzip to a convenient local folder.

# Command to launch TNT executable
tntExec <- "C:/Programs/Phylogeny/tnt/tnt.exe"

# Download MrBayes from 
# https://nbisweden.github.io/MrBayes/download.html
# and unzip/install to a convenient local folder.

# Command to launch MrBayes executable
mbExec <- "C:/Programs/Phylogeny/MrBayes/bin/mb.3.2.7-win64.exe"

# Download IQ-tree from http://www.iqtree.org/
# and unzip/install to a convenient local folder.

# Command to launch IQ-tree executable
iqExec <- "C:/Programs/Phylogeny/iqtree-2.2.0-Windows/bin/iqtree2.exe"

### Location of output files ###
# Set up directory structure
CreateDir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# Patterns to use when creating files
CreateDir("split-support/concordance")
ConcFile <- function(aln) {
  paste0("split-support/concordance/", aln, ".txt")
}

CreateDir("split-support/alignments")
DataFile <- function(aln, ext = ".nex") {
  paste0("split-support/alignments/", aln, ext)
}

CreateDir("split-support/iqtree")
IQFile <- function(aln) {
  paste0("split-support/iqtree/", aln, ".phy")
}

CreateDir("split-support/MrBayes")
MBFile <- function(aln, suffix = NULL) {
  paste0("split-support/MrBayes/", aln, if(!is.null(suffix)) ".", suffix)
}

CreateDir("split-support/TNT")
TNTFile <- function(aln, wt = "ew") {
  paste0("split-support/TNT/", aln, ".", wt, ".out")
}
