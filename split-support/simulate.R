set.seed(0)
referenceTree <- RootTree(ape::rtree(48), "t1")
plot(referenceTree)
write.tree("

for (i in 1:1000) {
  phangorn::simSeq(referenceTree, 192)
}
