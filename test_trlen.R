cat("1. Loading...\n")
devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

cat("2. Creating simple tree...\n")
tree <- ape::read.tree(text = "((1,2),3);")
cat("3. Creating dataset from StringToPhyDat...\n")
dataset <- StringToPhyDat("110 111", 1:3, byTaxon = FALSE)

cat("4. Trying to calculate TreeLength...\n")
flush.console()
result <- TreeLength(tree, dataset)
cat("Result:", result, "\n")
