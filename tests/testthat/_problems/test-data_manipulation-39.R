# Extracted from test-data_manipulation.R:39

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
expect_equal(Dehash(attributes(PrepareDataProfile(phy1))[1:10]),
               Dehash(attributes(phy1)))
mtx <- cbind(c("0", "0", 1,1,1,1),
               c(1,1,0,0,0,0),# flipped
               c(0,0,0,1,1,"{012}"))
rownames(mtx) <- letters[seq_len(nrow(mtx))]
phy2 <- TreeTools::MatrixToPhyDat(mtx)
expect_equivalent(phy1, PrepareDataProfile(phy2))
expect_equal(attributes(PrepareDataProfile(phy1)),
               attributes(PrepareDataProfile(phy2)))
