source("split-support/config.R")

ConcFile <- function(aln) {
  paste0("split-support/concordance/", aln, ".txt")
}

DataFile <- function(aln) {
  paste0("split-support/alignments/", aln, ".nex")
}

MBFile <- function(aln, suffix = NULL) {
  paste0("split-support/MrBayes/", aln, if(!is.null(suffix)) ".", suffix)
}

referenceTree <- read.tree("split-support/reference.tre")
refSplits <- as.Splits(referenceTree)

# Eugh, I don't like growing vectors like this!
partCorrect <- logical(0)
postProb <- numeric(0)
concord <- numeric(0)

for (aln in alns[1:21]) {
  
  parts <- read.table(MBFile(aln, "parts"), skip = 2 + nTip)
  partitions <- setNames(as.Splits(parts[, 2], tips), parts[, 1])
  truth <- partitions %in% refSplits
  
  pp <- read.table(MBFile(aln, "tstat"), skip = 1,
                   header = TRUE, comment.char = "")
  pp <- setNames(pp[, "Probability..s."], pp[, "ID"])
  
  dataset <- MatrixToPhyDat(matrix(unlist(read.nexus.data(DataFile(aln))), 
                                   nrow = nTip, byrow = TRUE,
                                   dimnames = list(tips, NULL)))
  
  if (file.exists(ConcFile(aln))) {
    conc <- as.matrix(read.table(ConcFile(aln)))
  } else {
    conc <- cbind(
      quartet = QuartetConcordance(partitions, dataset),
      cluster = ClusteringConcordance(partitions, dataset),
      phylo = PhylogeneticConcordance(partitions, dataset),
      mutual = MutualClusteringConcordance(partitions, dataset),
      shared = SharedPhylogeneticConcordance(partitions, dataset)
    )
    write.table(conc, ConcFile(aln))
  }
  
  partCorrect <- c(partCorrect, truth)
  postProb <- c(postProb, pp)
  concord <- rbind(concord, conc)
}

model <- glm(partCorrect ~ postProb + concord, family = "binomial")

model <- glm(family = "binomial",
             partCorrect ~ 
               postProb +
               concord[, "quartet"] +
               concord[, "cluster"] +
               concord[, "phylo"] +
               concord[, "mutual"] +
               concord[, "shared"]
             )
step(model) # AIC
# BIC: https://stackoverflow.com/questions/19400494
step(model, criterion = "BIC", k = log(length(partCorrect)))

# The lower the Brier score is for a set of predictions,
# the better the predictions are calibrated.
mclust::BrierScore(cbind(1 - postProb, postProb), partCorrect)
