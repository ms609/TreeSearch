library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
library("TreeSearch")
library("TreeDist")
library("Quartet", exclude = 'RobinsonFoulds')
nTip <- 8


generative <- RandomTree(nTip, TRUE)
generative$edge.length <- rgamma(dim(generative$edge)[1], shape = 1)
plot(generative)
dataBits <- lapply(phangorn::discrete.gamma(1, 4), function (rate)
  phangorn::simSeq(generative, l = 100, rate = rate))
dataset <- c(dataBits[[1]], dataBits[[2]], dataBits[[3]], dataBits[[4]]) # I can't remember how to do this right


morphyObj <- PhyDat2Morphy(dataset)
on.exit(morphyObj <- UnloadMorphy(morphyObj))

nTree <- NUnrooted(nTip)
trees <- lapply(seq_len(nTree) - 1L, as.phylo, nTip = nTip)

scores <- vapply(trees, function (tr) {
  c(ew = MorphyTreeLength(tr, morphyObj)
  )
}, c(ew = 0))

cid <- ClusteringInfoDistance(trees, generative)
pid <- PhylogeneticInfoDistance(trees, generative)
rf <- RobinsonFoulds(trees, generative)
qd <- QuartetDivergence(SingleTreeQuartetAgreement(trees, generative), similarity = FALSE)
Plot <- function (x, lab) {
  plot(-scores ~ x, xlab = lab, ylab = 'Parsimony score', pch = '.')
  legend('topright', bty = 'n',
         legend = signif(summary(lm(scores ~ x))$adj.r.squared))
}
par(mfrow = c(2, 2), mar = rep(2, 4))
Plot(cid, 'CID')
Plot(pid, 'PID')
Plot(qd, 'QD')
Plot(rf, 'RF')
