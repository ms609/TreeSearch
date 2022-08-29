# Set constants
nTip <- 48
nChar <- nTip * 2

# Set up directory structure
CreateDir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}
CreateDir("split-support/alignments")
CreateDir("split-support/tnt")
CreateDir("split-support/MrBayes")

# Create reference tree
set.seed(0)
referenceTree <- ape::rtree(nTip, equiprob = TRUE)
referenceTree$tip.label <- paste0("t", seq_len(nTip))
referenceTree <- RootTree(referenceTree, "t1")
treeLength <- sum(referenceTree$edge.length)
rate <- 12 / treeLength
print(signif(rate)) # mb.nex: prset brlenspr=unconstrained:uniform(0,<RATE>);
plot(referenceTree)
write.tree(referenceTree, file = "split-support/reference.tre")

# Simulate alignments
for (i in formatC(1:1000, width = 4, flag = "0")) {
  write.nexus.data(
    toupper(PhyDatToMatrix(
      phangorn::simSeq(referenceTree, nChar,
                       rootseq = rep("a", nChar),
                       rate = rate
                       ) # Jukes-Cantor model
    )),
    file = paste0("split-support/alignments/aln", i, ".nex"),
    interleaved = FALSE 
  )
}
