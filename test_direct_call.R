cat("Loading package...\n")
library(TreeSearch)

cat("Creating minimal test data...\n")
# Minimal 2-taxon tree
edge <- matrix(c(3, 3, 1, 2), nrow = 2, ncol = 2)
cat("Edge matrix:\n")
print(edge)

# Minimal contrast matrix (2 states)
contrast <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
cat("Contrast matrix:\n")
print(contrast)

# Minimal tip data (1 character, 2 taxa)
tip_data <- matrix(c(1, 2), nrow = 2, ncol = 1)
cat("Tip data:\n")
print(tip_data)

weight <- c(1)
levels <- c("0", "1")

cat("Now calling ts_fitch_score...\n")
flush.console()

result <- TreeSearch:::ts_fitch_score(edge, contrast, tip_data, weight, levels)
cat("Result:", result, "\n")
