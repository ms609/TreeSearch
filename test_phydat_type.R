cat("Loading...\n")
devtools::load_all(quiet = TRUE, compile = FALSE)
library(TreeTools)

tree <- ape::read.tree(text = "(((1,2),3),4);")

cat("Test with DNA dataset:\n")
# Create a simple DNA dataset
dna <- phangorn::phyDat(
  list(
    "1" = c("a", "c", "g"),
    "2" = c("c", "g", "t"),
    "3" = c("a", "a", "g"),
    "4" = c("t", "c", "t")
  ),
  type = "DNA"
)

cat("DNA dataset:\n")
print(dna)
cat("Calculating TreeLength for DNA...\n")
flush.console()
result_dna <- TreeLength(tree, dna)
cat("Result:", result_dna, "\n")

cat("\nTest with USER dataset created via phangorn:\n")
user_data <- phangorn::phyDat(
  list(
    "1" = c(FALSE, TRUE),
    "2" = c(TRUE, FALSE),
    "3" = c(FALSE, TRUE),
    "4" = c(TRUE, FALSE)
  ),
  type = "USER",
  levels = c(FALSE, TRUE)
)

cat("USER dataset:\n")
print(user_data)
cat("Calculating TreeLength for USER...\n")
flush.console()
result_user <- TreeLength(tree, user_data)
cat("Result:", result_user, "\n")
