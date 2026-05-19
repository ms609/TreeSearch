cat("Replicating TreeLength data construction path...\n")
library(TreeTools)

# Create the same dataset and tree as the failing test
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)

cat("Dataset:\n")
print(dataset)

cat("\nUnlisting dataset...\n")
unlisted <- unlist(dataset, use.names = FALSE)
cat("Unlisted values:", unlisted, "\n")
cat("Length of unlisted:", length(unlisted), "\n")
cat("nrow for matrix:", length(dataset), "\n")

cat("\nCreating matrix like TreeLength does...\n")
tip_data <- matrix(unlisted, nrow = length(dataset), byrow = TRUE)
cat("tip_data:\n")
print(tip_data)

cat("\nExtracting attributes like TreeLength does...\n")
at <- attributes(dataset)
contrast <- at$contrast
edge <- start_tree$edge

cat("Contrast:\n")
print(contrast)
cat("Edge:\n")
print(edge)

cat("\nCalling ts_fitch_score like TreeLength does...\n")
flush.console()
result <- TreeSearch:::ts_fitch_score(edge, contrast, tip_data, at$weight, at$levels)
cat("Result:", result, "\n")
