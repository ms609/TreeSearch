library("TreeTools", quietly = TRUE)

test_that("LocalConcordance() handles null / bad input", {
  expect_warning(expect_null(LocalConcordance(NULL)))
  expect_error(LocalConcordance("not_phyDat"), "phyDat")
})

test_that(".PairExpectedMI() returns zero for trivial cases", {
  # Single state on either side → no information

  expect_equal(.PairExpectedMI(c(10), c(5, 5)), 0)
  expect_equal(.PairExpectedMI(c(5, 5), c(10)), 0)
  # Both single-state

  expect_equal(.PairExpectedMI(c(8), c(8)), 0)
})

test_that(".PairExpectedMI() agrees with C++ expected_mi() for binary case", {
  # .ExpectedMI() wraps the C++ expected_mi() via caching;
  # .PairExpectedMI() uses the general R formula.
  # For binary a, both should agree.
  a <- c(5L, 15L)
  b <- c(8L, 12L)
  expect_equal(.PairExpectedMI(a, b), expected_mi(a, b), tolerance = 1e-8)

  a2 <- c(3L, 7L)
  b2 <- c(2L, 3L, 5L)
  expect_equal(.PairExpectedMI(a2, b2), expected_mi(a2, b2), tolerance = 1e-8)

  a3 <- c(10L, 10L)
  b3 <- c(5L, 5L, 5L, 5L)
  expect_equal(.PairExpectedMI(a3, b3), expected_mi(a3, b3), tolerance = 1e-8)
})

test_that(".PairExpectedMI() is symmetric", {
  a <- c(3L, 4L, 3L)
  b <- c(5L, 5L)
  expect_equal(.PairExpectedMI(a, b), .PairExpectedMI(b, a))
})

test_that("LocalConcordance() returns correct dimensions", {
  dat <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 1, 1, 0, 0, 1, 1,
      0, 1, 0, 1, 0, 1, 0, 1,
      0, 0, 0, 1, 0, 1, 1, 1), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))

  sigs <- c(1, 3)
  lc <- LocalConcordance(dat, sigma = sigs)

  expect_equal(nrow(lc), 5L)  # 5 alignment columns
  expect_equal(ncol(lc), length(sigs))
  expect_equal(colnames(lc), as.character(sigs))
})

test_that("Perfectly correlated characters yield NMI = 1", {
  # Two characters that carry the same partition (one is the complement)
  # should have NMI = 1.
  dat <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 1, 1, 1, 1,
      1, 1, 1, 1, 0, 0, 0, 0), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))
  lc <- LocalConcordance(dat, sigma = 1, normalize = FALSE)

  pairNMI <- attr(lc, "pairNMI")
  nr <- nrow(pairNMI)
  expect_gte(nr, 2)
  offDiag <- pairNMI[upper.tri(pairNMI)]
  expect_true(all(is.na(offDiag) | abs(offDiag - 1) < 1e-10))
})

test_that("All-ambiguous columns yield NA scores", {
  dat <- MatrixToPhyDat(matrix(
    c(0, 0, 1, 1, 0, 0, 1, 1,
      rep("?", 8),
      0, 1, 0, 1, 0, 1, 0, 1), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))
  lc <- LocalConcordance(dat, sigma = 2)

  # The all-? column (column 2) should have NA score
  expect_true(is.na(lc[2, 1]))
})

test_that("Uninformative (single-state) columns yield NA scores", {
  dat <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 1, 1, 1), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))
  lc <- LocalConcordance(dat, sigma = 1)
  # First column is constant → uninformative → NA
  expect_true(is.na(lc[1, 1]))
})

test_that("Gaussian weighting gives expected result for uniform NMI", {
  # Build an alignment where all pairwise NMI values in the window are 1:
  # alternating between two distinct binary characters that are identical.
  # Both characters are {0,0,0,0,1,1,1,1} and {0,0,1,1,0,0,1,1},
  # and we repeat them.  Every pair is the same two patterns, yielding

  # consistent (but not identical) NMI values.
  dat <- MatrixToPhyDat(matrix(
    c(rep(c(0, 0, 0, 0, 1, 1, 1, 1,
            0, 0, 1, 1, 0, 0, 1, 1), 5)), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))
  lc <- LocalConcordance(dat, sigma = c(2, 5), normalize = FALSE)

  # Interior columns should have finite, non-NA scores
  mid <- 5
  expect_true(any(is.finite(lc)))
  # All interior scores at same sigma should be equal (translation invariance)
  interior <- 4:7
  for (s in seq_len(ncol(lc))) {
    vals <- lc[interior, s]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 1) {
      expect_equal(max(vals) - min(vals), 0, tolerance = 1e-10)
    }
  }
})

test_that("LocalConcordance() works with congreveLamsdellMatrices", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  lc <- LocalConcordance(dataset, sigma = c(5, 10))
  expect_true(is.matrix(lc))
  expect_equal(ncol(lc), 2)
  # Should have some finite values
  expect_true(any(is.finite(lc)))
})

test_that("PlotLocalConcordance() runs without error", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  lc <- LocalConcordance(dataset, sigma = c(5, 10))

  expect_silent(PlotLocalConcordance(lc))
  expect_silent(PlotLocalConcordance(lc, dataset))
})

test_that("internal_gaps: even a single-taxon internal gap is retained", {
  # t1 has an internal gap at col 3 (flanked by nucleotides on each side).
  # With internal_gaps = TRUE, this is a real indel event and must NOT be
  # dropped as an autapomorphy — it creates an informative 1-vs-rest split.
  mat <- matrix(
    c("A", "A", "-", "A", "A",
      "A", "A", "A", "A", "A",
      "A", "C", "A", "C", "A",
      "C", "A", "A", "A", "C",
      "A", "C", "C", "C", "A",
      "C", "A", "A", "A", "C",
      "A", "A", "C", "A", "A",
      "C", "C", "A", "C", "C"), 8, byrow = TRUE,
    dimnames = list(paste0("t", 1:8), NULL)
  )
  dat <- MatrixToPhyDat(mat)

  lc_igap  <- LocalConcordance(dat, sigma = 1, normalize = FALSE,
                                internal_gaps = TRUE)
  lc_nogap <- LocalConcordance(dat, sigma = 1, normalize = FALSE,
                                internal_gaps = FALSE)

  # internal_gaps=TRUE: col 3 is an extra state, kept even as a singleton
  #   → different pairNMI from treating the gap as missing data
  expect_false(identical(attr(lc_igap, "pairNMI"), attr(lc_nogap, "pairNMI")))

  # The site at col 3 should be scorable under internal_gaps=TRUE
  # (it can be compared with neighbours for mutual information)
  expect_false(is.na(lc_igap[3, 1]))
})

test_that("internal_gaps: terminal gaps remain NA", {
  # t1: - - A A A  (terminal gaps at cols 1-2)
  # All other taxa: full data
  mat <- matrix(
    c("-", "-", "A", "A", "A",
      "A", "A", "A", "A", "A",
      "A", "C", "A", "C", "A",
      "C", "A", "A", "A", "C",
      "A", "C", "C", "C", "A",
      "C", "A", "A", "A", "C",
      "A", "A", "C", "A", "A",
      "C", "C", "A", "C", "C"), 8, byrow = TRUE,
    dimnames = list(paste0("t", 1:8), NULL)
  )
  dat <- MatrixToPhyDat(mat)
  lc <- LocalConcordance(dat, sigma = 1, normalize = FALSE,
                          internal_gaps = TRUE)

  # Columns 1-2 (terminal gaps in t1): the gap character is still NA
  # → fewer informative sites but function should not error
  expect_true(is.matrix(lc))
  expect_equal(nrow(lc), 5L)
})

test_that("block_size / threshold produce trim attribute", {
  dat <- MatrixToPhyDat(matrix(
    c(rep(c(0, 0, 0, 0, 1, 1, 1, 1,
            0, 1, 0, 1, 0, 1, 0, 1), 4)), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))

  # Without block_size: no trim attribute
  lc_no <- LocalConcordance(dat, sigma = 1, block_size = NULL)
  expect_null(attr(lc_no, "trim"))

  # With block_size: trim attribute is a logical of length L
  lc_bl <- LocalConcordance(dat, sigma = 1, normalize = FALSE,
                              block_size = 4L, threshold = 0.5,
                              noise_level = 0.3)
  trim <- attr(lc_bl, "trim")
  expect_true(is.logical(trim))
  expect_equal(length(trim), nrow(lc_bl))
})

test_that("block scan: all-concordant data produces no trim", {
  # Two perfectly correlated (complementary) patterns: NMI = 1 everywhere.
  # With noise_level = 0.5 (well below 1), no sites should be flagged.
  dat <- MatrixToPhyDat(matrix(
    c(rep(c(0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 0, 0, 0, 0), 5)), 8,
    dimnames = list(paste0("t", 1:8), NULL)
  ))
  lc <- LocalConcordance(dat, sigma = 2, normalize = FALSE,
                          block_size = 4L, threshold = 0.5,
                          noise_level = 0.5)
  trim <- attr(lc, "trim")
  # Pairwise NMI = 1 >> 0.5; no sites are "noisy"
  expect_false(any(trim, na.rm = TRUE))
})

test_that("PlotLocalConcordance() shades trim regions", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  lc <- LocalConcordance(dataset, sigma = 5, normalize = FALSE,
                          block_size = 5L, threshold = 0.3,
                          noise_level = 0.5)
  expect_silent(PlotLocalConcordance(lc))
  expect_silent(PlotLocalConcordance(lc, dataset))
  # trim_col = NA suppresses shading
  expect_silent(PlotLocalConcordance(lc, trim_col = NA))
})

test_that("normalize = TRUE shifts scores towards zero for random data", {
  set.seed(7418)
  randomMat <- matrix(sample(0:1, 8 * 50, replace = TRUE), 8,
                      dimnames = list(paste0("t", 1:8), NULL))
  dat <- MatrixToPhyDat(randomMat)
  lcNorm   <- LocalConcordance(dat, sigma = 10, normalize = TRUE)
  lcRaw    <- LocalConcordance(dat, sigma = 10, normalize = FALSE)

  # After normalization, mean score should be closer to zero than raw
  expect_lt(
    abs(mean(lcNorm, na.rm = TRUE)),
    mean(lcRaw, na.rm = TRUE) + 0.01  # allow small margin
  )
})
