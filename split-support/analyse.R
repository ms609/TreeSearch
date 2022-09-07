# Load required libraries
#library("TreeTools")
#library("TreeSearch")
devtools::load_all("../TreeTools")
if (packageVersion("TreeTools") < "1.7.3.9004") {
  stop("Please upgrade to TreeTools v1.8.0 or above")
}
devtools::load_all("../TreeSearch")

# Load configuration settings
source("split-support/config.R")

referenceTree <- read.tree("split-support/reference.tre")
refSplits <- as.Splits(referenceTree)
tips <- names(read.nexus.data(DataFile("aln0001")))

# Eugh, I don't like growing vectors like this!
partCorrect <- logical(0)
postProb <- numeric(0)
concord <- numeric(0)
bremer <- numeric(0)
tntStats <- c("symFq", "symGC", "boot", "jak", "pois")
tntStat <- matrix(0, 0, length(tntStats), dimnames = list(NULL, tntStats))
ufb <- numeric(0)
iqStats <- c("alrt", "lbp", "abayes", "ufb") # .iqtree output file gives order
iqStat <- matrix(0, 0, length(iqStats), dimnames = list(NULL, iqStats))

for (i in cli::cli_progress_along(seq_len(nAln), "Analysing")) {
  aln <- alns[i]
  
  # Load MrBayes partitions
  parts <- read.table(MBFile(aln, "parts"), skip = 2 + nTip)
  partitions <- setNames(as.Splits(parts[, 2], tips), paste0("mb", parts[, 1]))
  
  # Load TNT partitions
  tntFile <- TNTFile(aln, "ew")
  tntTree <- suppressMessages(ReadTntTree(tntFile, tipLabels = tips))
  if (!inherits(tntTree, "multiPhylo")) {
    warning("Only one tree found in file ", aln, "; missing TNT output.")
    next;
  }
  if (!all.equal(tntTree[[1]], tntTree[[2]])) {
    # Trees may differ in resolution: partitions with 0 Bremer support will
    # be collapsed in tntTree[[1]].
    if (any(!as.Splits(tntTree[[1]]) %in% as.Splits(tntTree[[2]]))) {
      warning("Trees don't match in ", aln, ". Check TNT output.")
      next;
    }
  }
  tntParts <- as.Splits(tntTree[[2]])
  nTntNode <- length(tntParts)
  tntOnly <- !tntParts %in% partitions
  if (any(tntOnly)) {
    partitions <- c(
      partitions,
      setNames(tntParts[[tntOnly]], paste0("tnt", seq_len(sum(tntOnly))))
    )
  }
  
  # Load IQ-tree partitions
  iqTree <- read.tree(IQFile(aln, ".treefile"))
  iqParts <- as.Splits(iqTree)
  iqOnly <- !iqParts %in% partitions
  if (any(iqOnly)) {
    partitions <- c(
      partitions,
      setNames(iqParts[[iqOnly]], paste0("iq", seq_len(sum(iqOnly))))
    )
  }
  ufbLines <- readLines(IQFile(aln, ".splits.nex"))[-seq_len(nTip + 13)]
  ufbLines <- ufbLines[seq_len(which.max(ufbLines == ";") - 1)]
  ufbLines <- do.call(rbind, strsplit(trimws(ufbLines), "\t"))
  ufbParts <- as.Splits(t(vapply(
    strsplit(trimws(gsub(",", "", fixed = TRUE, ufbLines[, 2])), " "),
    function(a) {
      tabulate(as.numeric(a), nTip) == 1
    }, logical(nTip))), tipLabels = tips)
  trivial <- TrivialSplits(ufbParts)
  ufbVals <- as.numeric(ufbLines[!trivial, 1])
  ufbParts <- ufbParts[[!trivial]]
  ufbOnly <- !ufbParts %in% partitions
  if (any(ufbOnly)) {
    partitions <- c(
      partitions,
      setNames(ufbParts[[ufbOnly]], paste0("ufb", seq_len(sum(ufbOnly))))
    )
  }
  
  # Once all partitions are loaded, label where possible
  
  # Populate TNT supports
  brem <- rep(NA_real_, length(partitions))
  partId2 <- match(as.Splits(tntTree[[2]]), partitions)
  brem[partId2] <- 0
  partBrem <- tntTree[[1]]$node.label[as.numeric(names(as.Splits(tntTree[[1]]))) - NTip(tntTree[[1]])]
  brem[match(as.Splits(tntTree[[1]]), partitions)] <- as.numeric(partBrem)
  
  tags <- strsplit(tntTree[[2]]$node.label, "/")
  partTags <- tags[as.numeric(names(as.Splits(tntTree[[2]]))) - NTip(tntTree[[2]])]
  tntTags <- matrix(NA_real_, length(partitions), length(tntStats),
                    dimnames = list(NULL, tntStats))
  tntTags[partId2, ] <- t(vapply(partTags, function(tag) {
    x <- gsub("[", "-", fixed = TRUE, gsub("]", "", fixed = TRUE, tag))
    x[x == "?"] <- NA_real_
    x[x == "???"] <- Inf
    as.numeric(x)
  }, numeric(length(tntStats))))
  
  
  # Populate IQ-tree supports
  iqMatch <- match(iqParts, partitions)
  tags <- strsplit(iqTree[["node.label"]], "/")
  partTags <- tags[as.numeric(names(iqParts)) - NTip(iqTree)]
  iqTags <- matrix(NA_real_, length(partitions), length(iqStats),
                   dimnames = list(NULL, iqStats))
  iqTags[iqMatch, ] <- t(vapply(partTags, as.numeric, numeric(length(iqStats))))
  
  
  # Populate Ultra-Fast bootstrap supports for partitions not in consensus
  ufbMatch <- match(ufbParts, partitions)
  iqTags[ufbMatch, "ufb"] <- ufbVals
  
  
  # Populate posterior probabilities
  pp <- read.table(MBFile(aln, "tstat"), skip = 1,
                   header = TRUE, comment.char = "")
  pp <- setNames(pp[, "Probability..s."], pp[, "ID"])
  
  
  # Calculate concordances
  dataset <- MatrixToPhyDat(matrix(unlist(read.nexus.data(DataFile(aln))), 
                                   nrow = nTip, byrow = TRUE,
                                   dimnames = list(tips, NULL)))
  
  if (file.exists(ConcFile(aln))) {
    conc <- as.matrix(read.table(ConcFile(aln)))
    if (dim(conc)[1] != dim(tntTags)[1]) {
      file.remove(ConcFile(aln))
      stop("Dimension mismatch; is concordance cache ", aln, " out of date?")
    }
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
  postProb <- c(postProb, pp, rep(0, sum(tntOnly, iqOnly, ufbOnly)))
  concord <- rbind(concord, conc)
  bremer <- c(bremer, brem)
  tntStat <- rbind(tntStat, tntTags)
  iqStat <- rbind(iqStat, iqTags)
  
  stopifnot(dim(concord)[1] == dim(tntStat)[1])
  stopifnot(dim(concord)[1] == length(postProb))
}

model <- glm(partCorrect ~ postProb + concord + bremer + tntStat + iqStat,
             family = "binomial")

Histy <- function(var, breaks = 12, even = TRUE) { # "Mosaic plot"
  outcomes <- partCorrect[!is.na(var)]
  var <- var[!is.na(var)]
  if (even) {
    breaks <- quantile(var, seq(0, 1, length.out = breaks))
  }
  bins <- cut(var, breaks = unique(breaks))
  plot(table(bins, outcomes),
       main = as.character(match.call()[-1]),
       col = c("FALSE" = 2, "TRUE" = 3),
       xlab = "",
       ylab = ""
       )
  axis(1, signif(breaks), at = seq_along(breaks) / 12)
}

par(mfrow = c(4, 4), mar = rep(2, 4))
Histy(postProb)
Histy(concord[, "quartet"])
Histy(concord[, "mutual"])
Histy(concord[, "shared"])
Histy(concord[, "phylo"])
Histy(concord[, "cluster"])
Histy(bremer)
Histy(tntStat[, "symFq"])
Histy(tntStat[, "symGC"])
Histy(tntStat[, "boot"])
Histy(tntStat[, "jak"])
Histy(tntStat[, "pois"])
Histy(iqStat[, "ufb"])
Histy(iqStat[, "lbp"])
Histy(iqStat[, "alrt"])
Histy(iqStat[, "abayes"])

Peek <- function(var) {
  m <- glm(partCorrect ~ var, family = "binomial")
  smry <- summary(m)
  print(smry)
  print(paste("R2:", signif(1 - smry$deviance / smry$null.deviance)))
  AIC(m)
}
Peek(postProb)
Peek(bremer)
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
               bremer[common] +
               concord[common, "quartet"] +
               concord[common, "cluster"] +
               concord[common, "phylo"] +
               concord[common, "mutual"] +
               concord[common, "shared"] +
               tntStat[common, "symFq"] +
               tntStat[common, "symGC"] +
               tntStat[common, "boot"] +
               tntStat[common, "jak"] + 
               tntStat[common, "pois"]
             )
step(model) # AIC
# BIC: https://stackoverflow.com/questions/19400494
step(model, criterion = "BIC", k = log(length(partCorrect)))



# The lower the Brier score is for a set of predictions,
# the better the predictions are calibrated.
mclust::BrierScore(cbind(1 - postProb, postProb), partCorrect)
