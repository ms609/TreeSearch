nTip <- 9
nRep <- 10

message(nTip, " tips; ", nRep, " replications.")
source('data-raw/small-tree.R')

results <- vapply(seq_len(nRep), CompareMethods, nTip = nTip, matrix(0, 5, 5))
write.csv(results, file = paste0('results-', nTip, '-', nRep, '.csv'))

source('data-raw/plot-small-tree.R')