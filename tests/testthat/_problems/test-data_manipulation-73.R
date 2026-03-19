# Extracted from test-data_manipulation.R:73

# prequel ----------------------------------------------------------------------
Dehash <- function (x) {
  lapply(x, function (xi) {
    attr(xi, ".match.hash") <- NULL
    if (!is.null(dimnames(xi))) {
      dimnames(xi) <- Dehash(dimnames(xi))
    }
    xi
  })
}

# test -------------------------------------------------------------------------
mtx <- cbind(c("0", "0", 1,1,1,1),
               c(0,0,1,1,1,1),# again
               c(0,0,0,1,1,"?"))
rownames(mtx) <- letters[seq_len(nrow(mtx))]
phy1 <- TreeTools::MatrixToPhyDat(mtx)
expect_equivalent(phy1, PrepareDataProfile(phy1))
pp1 <- PrepareDataProfile(phy1)
expect_equal(attr(pp1, "weight"), attr(phy1, "weight"))
expect_equal(attr(pp1, "nr"), attr(phy1, "nr"))
expect_equal(attr(pp1, "nc"), attr(phy1, "nc"))
expect_equal(attr(pp1, "index"), attr(phy1, "index"))
mtx <- cbind(c("0", "0", 1,1,1,1),
               c(1,1,0,0,0,0),# flipped
               c(0,0,0,1,1,"{012}"))
rownames(mtx) <- letters[seq_len(nrow(mtx))]
phy2 <- TreeTools::MatrixToPhyDat(mtx)
pp2 <- PrepareDataProfile(phy2)
expect_equal(attr(pp2, "nr"), 3L)
expect_equal(attr(pp2, "nc"), attr(pp1, "nc"))
expect_equal(attr(pp2, "info.amounts")[, 1], attr(pp1, "info.amounts")[, 1])
mtx <- cbind(c("0", "0", 1,1,1, "2", "2", 3,3,3,3),
               c("?", "?", 1,1,1, "?", "?", 0,0,0,0),
               c(0,0,1,1,1,2,2,3,3,3,3),# again
               c(rep("?", 5), "2", "2", 0,0,0,0),
               c("?", "?", 1,1,1, 1,1, 0,0,0,0),
               c("0", "1", rep("?", 9))
               )
rownames(mtx) <- letters[seq_len(nrow(mtx))]
dataset <- TreeTools::MatrixToPhyDat(mtx)
q <- "?"
reduced <- matrix(c(q,q,q,q,q,q,q,q,q,q,q,  # char 1 col (binary-reduced)
                      q,q,1,1,1,q,q,q,q,q,q,    # char 2 (was "a" col)
                      1,1,1,1,1,1,1,2,2,2,2,    # remainder
                      q,q,q,q,q,q,q,q,q,q,q,
                      q,q,q,q,q,q,q,q,q,q,q,
                      q,q,q,q,q,q,q,q,q,q,q),
                    ncol = 6, dimnames = list(letters[1:11], NULL))
expect_warning(pd <- PrepareDataProfile(dataset),
                 "Multi-state characters reduced")
