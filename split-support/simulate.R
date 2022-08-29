set.seed(0)
referenceTree <- RootTree(ape::rtree(48), "t1")
plot(referenceTree)
write.tree(referenceTree, file = "split-support/reference.tre")
dir.create("split-support/alignments")
dir.create("split-support/tnt")
dir.create("split-support/MrBayes")

for (i in 1:1000) {
  write.nexus.data(
    phangorn::simSeq(referenceTree, 288), # Jukes-Cantor model
    file = paste0("split-support/alignments/aln",
                  formatC(i, width = 4, flag = "0"),
                  ".nex")
  )
}
