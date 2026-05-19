cat("Loading...\n")
devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

cat("Creating test data...\n")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)

cat("Dataset attributes:\n")
at <- attributes(dataset)
cat("- levels:", at$levels, "\n")
cat("- nr (num chars):", at$nr, "\n")
cat("- weight:", at$weight, "\n")

cat("\nContrast matrix:\n")
contrast <- at$contrast
print(contrast)
cat("Contrast dimensions:", nrow(contrast), "x", ncol(contrast), "\n")

cat("\nDataset:\n")
print(dataset)

cat("\nCreating tree...\n")
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("Tree edge matrix:\n")
print(start_tree$edge)
cat("Tree tip labels:", start_tree$tip.label, "\n")

cat("\nNow attempting TreeLength...\n")
flush.console()

# Try to get the internal data structures
cat("Preparing data for C++ call...\n")
edge <- start_tree$edge
tip_data <- matrix(unlist(dataset, use.names = FALSE), nrow = length(dataset), byrow = TRUE)
cat("tip_data dimensions:", nrow(tip_data), "x", ncol(tip_data), "\n")
cat("tip_data:\n")
print(tip_data)

cat("\nCalling TreeLength...\n")
flush.console()
result <- TreeLength(start_tree, dataset)
cat("Result:", result, "\n")
