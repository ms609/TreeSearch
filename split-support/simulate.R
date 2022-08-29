# Set constants
nTip <- 48

# Set up directory structure
CreateDir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}
CreateDir("split-support/alignments")
CreateDir("split-support/tnt")
CreateDir("split-support/MrBayes")

# Create reference tree
set.seed(0)
referenceTree <- ape::rtree(nTip)
referenceTree$tip.label <- paste0("t", seq_len(nTip))
referenceTree <- RootTree(referenceTree, "t1")
plot(referenceTree)
write.tree(referenceTree, file = "split-support/reference.tre")

# Simulate alignments
for (i in 1:1000) {
  write.nexus.data(
    toupper(PhyDatToMatrix(
      phangorn::simSeq(referenceTree, nTip * 10) # Jukes-Cantor model
    )), 
    file = paste0("split-support/alignments/aln",
                  formatC(i, width = 4, flag = "0"),
                  ".nex")
  )
}
