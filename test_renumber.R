cat("Testing Renumber function...\n")
library(TreeTools)

start_tree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), c("1", "2", "3", "4", "5", "6"))

cat("Original tree:\n")
print(start_tree)
cat("Original edge:\n")
print(start_tree$edge)

cat("\nCalling Renumber...\n")
renumbered <- Renumber(start_tree)

cat("Renumbered tree:\n")
print(renumbered)
cat("Renumbered edge:\n")
print(renumbered$edge)

cat("\nCalling RenumberTips...\n")
renumbered_tips <- RenumberTips(renumbered, c("1", "2", "3", "4", "5", "6"))

cat("After RenumberTips edge:\n")
print(renumbered_tips$edge)
