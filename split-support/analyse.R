# Load required libraries
library("TreeTools")

# Load configuration settings
source("split-support/config.R")


referenceTree <- read.tree("split-support/reference.tre")
refSplits <- as.Splits(referenceTree)

# Eugh, I don't like growing vectors like this!
partCorrect <- logical(0)
postProb <- numeric(0)
concord <- numeric(0)
tntStat <- matrix(0, 0, 2, dimnames = list(NULL, c("sym", "freq")))

for (aln in alns[1:242]) {
  
  parts <- read.table(MBFile(aln, "parts"), skip = 2 + nTip)
  partitions <- setNames(as.Splits(parts[, 2], tips), paste0("mb", parts[, 1]))
  
  tntFile <- TNTFile(aln, "ew")
  tntTree <- ReadTntTree(tntFile, tipLabels = tips)
  tntParts <- as.Splits(tntTree)
  tntOnly <- !tntParts %in% partitions
  if (any(tntOnly)) {
    message("Unique splits in TNT tree ", aln)
    partitions <- c(partitions, tntParts[[tntOnly]])
  }
  tags <- strsplit(
    read.table(tntFile, skip = 6 + nTip, comment.char = ";")[, 3],
    "/"
  )
  tntTags <- t(vapply(tags, function(tag) {
    as.numeric(gsub("[", "-", fixed = TRUE,
                    gsub("]", "", fixed = TRUE, tag)))
    }, c(sym = 0, freq = 0)))
  
  
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
  
  partCorrect <- c(partCorrect, partitions %in% refSplits)
  postProb <- c(postProb, pp, rep(NA_real_, sum(tntOnly)))
  concord <- rbind(concord, conc)
  tntTmp <- matrix(NA_real_, length(partitions), dim(tntTags)[2])
  tntTmp[match(tntParts, partitions), ] <- tntTags
  tntStat <- rbind(tntStat, tntTmp)
}

model <- glm(partCorrect ~ postProb + concord + tntStat, family = "binomial")

Peek <- function(var) {
  m <- glm(partCorrect ~ var, family = "binomial")
  smry <- summary(m)
  print(smry)
  print(paste("R2:", signif(1 - smry$deviance / smry$null.deviance)))
  AIC(m)
}
Peek(postProb)
Peek(concord[, "quartet"])
Peek(concord[, "cluster"])
Peek(concord[, "mutual"])
  
m <- glm(partCorrect ~ concord[, "quartet"] + postProb, family = "binomial")
AIC(m)
m <- glm(partCorrect ~ concord[, "quartet"], family = "binomial")

common <- rowSums(is.na(concord)) == 0 & rowSums(is.na(tntStat)) == 0
model <- glm(family = "binomial",
             partCorrect[common] ~ 
               postProb[common] +
               concord[common, "quartet"] + # Dropped
               concord[common, "cluster"] +
               concord[common, "phylo"] +
               concord[common, "mutual"] +
               concord[common, "shared"] +
               tntStat[common, "sym"] + # Dropped
               tntStat[common, "freq"] # Dropped
             )
step(model) # AIC
# BIC: https://stackoverflow.com/questions/19400494
step(model, criterion = "BIC", k = log(length(partCorrect)))



# The lower the Brier score is for a set of predictions,
# the better the predictions are calibrated.
mclust::BrierScore(cbind(1 - postProb, postProb), partCorrect)
