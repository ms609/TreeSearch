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
iqDir <- "split-support/iqtree/"
mbDir <- "split-support/MrBayes/"
tntDir <- "split-support/TNT/"
alnDir <- "split-support/alignments/"
concDir <- "split-support/concordance/"

# Set up directory structure
CreateDir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

# Patterns to use when creating files
CreateDir(concDir)
ConcFile <- function(aln) {
  paste0(concDir, "/", aln, ".txt")
}

CreateDir(alnDir)
DataFile <- function(aln, ext = ".nex") {
  paste0(alnDir, "/", aln, ext)
}

CreateDir(iqDir)
IQFile <- function(aln, suffix = "") {
  paste0(iqDir, "/", aln, ".phy", suffix)
}

CreateDir(mbDir)
MBFile <- function(aln, suffix = NULL) {
  paste0(mbDir, "/", aln, if(!is.null(suffix)) ".", suffix)
}

CreateDir(tntDir)
TNTFile <- function(aln, wt = "ew") {
  paste0(tntDir, "/", aln, ".", wt, ".out")
}
