set.seed(0)
referenceTree <- RootTree(ape::rtree(48), "t1")
plot(referenceTree)
write.tree(referenceTree, file = "split-support/reference.tre")
CreateDir <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}
CreateDir("split-support/alignments")
CreateDir("split-support/tnt")
CreateDir("split-support/MrBayes")

for (i in 1:1000) {
  write.nexus.data(
    phangorn::simSeq(referenceTree, 288), # Jukes-Cantor model
    file = paste0("split-support/alignments/aln",
                  formatC(i, width = 4, flag = "0"),
                  ".nex")
  )
}
