# Extracted from test-Morphy.R:53

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "TreeSearch", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library("TreeTools", quietly = TRUE)

# test -------------------------------------------------------------------------
constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
characters <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
set.seed(0)
ewResults <- Morphy(characters,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint)
expect_equal(PectinateTree(letters[1:6]), ewResults[[1]])
expect_equal(c(seed = 0, start = 1, final = 0),
               attr(ewResults, "firstHit"))
expect_equal(names(ewResults), "start_1")
expect_equal(PectinateTree(letters[1:6]),
               Morphy(characters, concavity = "p",
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint)[[1]])
expect_equal(PectinateTree(letters[1:6]),
               Morphy(characters, concavity = 10,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint)[[1]])
dataset <- characters
tree <- PectinateTree(c("a", "c", "f", "d", "e", "b"))
expect_equal(PectinateTree(letters[1:6]),
               Morphy(characters,
                      PectinateTree(c("a", "c", "f", "d", "e", "b")),
                                 ratchIter = 0, constraint = constraint)[[1]])
dataset <- MatrixToPhyDat(matrix(c(0, 0, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 0, 0, 0), ncol = 2,
                                   dimnames = list(letters[1:7], NULL)))
constraint <- MatrixToPhyDat(matrix(c(0, 0, 1, "?", 1, 1,
                                        1, 1, 1,   1, 0, 0), ncol = 2,
                                      dimnames = list(letters[1:6], NULL)))
cons <- Consensus(Morphy(dataset, constraint = constraint),
                    rooted = TRUE)
expect_true(as.Splits(as.logical(c(0, 0, 1, 1, 1)), letters[c(1:3, 5:6)]) %in% 
                as.Splits(DropTip(cons, c("d", "g"))))
