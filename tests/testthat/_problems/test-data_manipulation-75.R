# Extracted from test-data_manipulation.R:75

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
decomposed <- matrix(c(0,0,q,q,q,q,q,1,1,1,1,
                         q,q,0,0,0,q,q,1,1,1,1,
                         q,q,q,q,q,0,0,1,1,1,1,
                         
                         q,q,0,0,0,q,q,1,1,1,1,
                         
                         0,0,q,q,q,q,q,1,1,1,1,
                         q,q,0,0,0,q,q,1,1,1,1,
                         q,q,q,q,q,0,0,1,1,1,1,
                         
                         q,q,q,q,q,0,0,1,1,1,1,
                         q,q,0,0,0,0,0,1,1,1,1),
                       ncol = 9, dimnames = list(letters[1:11], NULL))
expect_warning(pd <- PrepareDataProfile(dataset))
expect_equal(decomposed, PhyDatToMatrix(pd))
expect_equal(c(1, 2, 3, 2, 1, 2, 3, 3, 4), attr(pd, "index"))
expect_equal(c(2, 3, 3, 1), attr(pd, "weight"))
dataset2 <- TreeTools::MatrixToPhyDat(mtx[!mtx[, 1] %in% c(0, 2), ])
expect_equal(attr(PrepareDataProfile(dataset2), "info.amounts"),
               attr(pd, "info.amounts")[1:3, 2, drop = FALSE])
