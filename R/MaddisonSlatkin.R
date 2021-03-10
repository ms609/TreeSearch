
#Remember <- function (X) X
Remember <- addMemoization

# D: The probability that, in a randomly selected tree of _n_ taxa, _i_ of
# which have state 1, the smaller subclade of _m_ taxa will receive _j_ taxa
# with state 1
.D <- function (j, i, m, n) {
  choose(i, j) * choose(n - i, m - j) / choose(n, m)
}
# R: The probability that in a randomly selected tree of _n_ leaves, the
# smaller of the two basal subclades will have _m_ leaves.
# My function gives half the value given in Maddison & Slatkin 1991 [!]
.R <- Remember(function (m, n) {
  prod(ifelse(n == m + m, 1L, 2L) / 2L, 
       choose(n, m),
       NRooted(m),
       NRooted(n - m)
  ) / NRooted(n)
})

# B{b}: The probability that state _b_ is reconstructed at the base of a clade
# with _n_ leaves and _i_ leaves with state 1.
# My function runs from j = 0 to i, not j = 1 to i as printed in M&S91 [!]
.B0 <- Remember(function (n, i) {
  if (n == 2L && i == 1L) {
    0
  } else if (n == i) {
    0
  } else if (i == 0L) {
    1
  } else {
    sum(vapply(seq_len(n / 2L), function (m) {
      .R(m, n) * sum(vapply(0:(min(i, m)), function (j) {
        .D(j, i, m, n) * sum(
          .B0(m, j) * .B0(n - m, i - j),
          .B01(m, j) * .B0(n - m, i - j),
          .B0(m, j) * .B01(n - m, i - j)
        )
      }, double(1)))
    }, double(1)))
  }
})

.B1 <- Remember(function (n, i) {
  if (n == 2L && i == 1L) {
    0
  } else if (n == i) {
    1
  } else if (i == 0L) {
    0
  } else {
    sum(vapply(seq_len(n / 2L), function (m) {
      .R(m, n) * sum(vapply(0:min(i, m), function (j) {
        .D(j, i, m, n) * sum(
          .B1(m, j) * .B1(n - m, i - j),
          .B01(m, j) * .B1(n - m, i - j),
          .B1(m, j) * .B01(n - m, i - j)
        )
      }, double(1)))
    }, double(1)))
  }
})

.B01 <- Remember(function (n, i) {
  if (n == 2L && i == 1L) {
    1
  } else if (n == i) {
    0
  } else if (i == 0L) {
    0
  } else {
    sum(vapply(seq_len(n / 2L), function (m) {
      .R(m, n) * sum(vapply(0:min(i, m), function (j) {
        .D(j, i, m, n) * sum(
          .B01(m, j) * .B01(n - m, i - j),
          .B0(m, j) * .B1(n - m, i - j),
          .B1(m, j) * .B0(n - m, i - j)
        )
      }, double(1)))
    }, double(1)))
  }
})

.PRow <- function (PSmall, BSmall, PLarge, BLarge, i, j, m, n, s) {
  sum(vapply(0:s, function (r) prod(
      PSmall(r, m, j),
      BSmall(m, j),
      PLarge(s - r, n - m, i - j),
      BLarge(n - m, i - j)),
    double(1)))
}


# P: The probability of having _s_ steps in a clade of _n_ leaves given that
# _i_ of the leaves have state 1.
# P{b}: P, given state _b_ is reconstructed at the base of the clade.
# My function runs from j = 0 to i, not j = 1 to i as printed in M&S91 [!]
.P0 <- Remember(function (s, n, i) {
  if (n == 2L && i == 1L) {
    if (s == 1L) 1 else 0
  } else if (n == i) {
    0 # Undefined!
  } else if (i == 0L) {
    if (s == 0L) 1 else 0
  } else {
    b0 <- .B0(n, i)
    if (b0 == 0) 0 else {
      sum(vapply(seq_len(n / 2), function (m) {
        .R(m, n) / b0 * sum(vapply(0:min(i, m), function (j) {
          .D(j, i, m, n) * sum(
            .PRow(.P0, .B0, .P0, .B0, i, j, m, n, s),
            .PRow(.P0, .B0, .P01, .B01, i, j, m, n, s),
            .PRow(.P01, .B01, .P0, .B0, i, j, m, n, s)
          )
        }, double(1)))
      }, double(1)))
    }
  }
})

.P1 <- Remember(function (s, n, i) {
  if (n == 2L && i == 1L) {
    if (s == 1L) 1 else 0
  } else if (n == i) {
    if (s == 0L) 1 else 0
  } else if (i == 0L) {
    0 # Undefined!
  } else {
    b1 <- .B1(n, i)
    if (b1 == 0) 0 else {
      sum(vapply(seq_len(n / 2), function (m) {
        .R(m, n) / b1 * sum(vapply(0:min(i, m), function (j) {
          .D(j, i, m, n) * sum(
            .PRow(.P1, .B1, .P1, .B1, i, j, m, n, s),
            .PRow(.P1, .B1, .P01, .B01, i, j, m, n, s),
            .PRow(.P01, .B01, .P1, .B1, i, j, m, n, s)
          )
        }, double(1)))
      }, double(1)))
    }
  }
})

.P01 <- Remember(function (s, n, i) {
  if (n == 2L && i == 1L) {
    if (s == 1L) 1 else 0
  } else if (n == i) {
    0 # Undefined?
  } else if (i == 0L) {
    0 # Undefined?
  } else {
    b01 <- .B01(n, i)
    if (b01 == 0) 0 else {
      sum(vapply(seq_len(n / 2), function (m) {
        .R(m, n) / b01 * sum(vapply(0:min(i, m), function (j) {
          .D(j, i, m, n) * sum(
            .PRow(.P01, .B01, .P01, .B01, i, j, m, n, s),
            .PRow(.P1, .B1, .P0, .B0, i, j, m, n, s - 1),
            .PRow(.P0, .B0, .P1, .B1, i, j, m, n, s - 1)
          )
        }, double(1)))
      }, double(1)))
    }
  }
})

#' Caluclate number of trees with a given score using Maddison & Slatkin's
#' recursive approach 
#' 
#' @param a,b Number of leaves in tree with state 1 / 2
#' 
#' @references
#' \insertRef{Maddison1991}{TreeSearch}
#' @examples
#' \dontrun{MaddisonSlatkin(2, 4)}
#' @template MRS
#' @importFrom R.cache addMemoization
#' @importFrom TreeTools NRooted
#' @family profile parsimony functions
#' @export
MaddisonSlatkin <- function (a, b) {
  n <- sum(a, b)
  i <- min(a, b)
  
  
  # Return:
  vapply(seq_len(i), function (s)
    sum(
      .P0(s, n, i) * .B0(n, i),
      .P1(s, n, i) * .B1(n, i),
      .P01(s, n, i) * .B01(n, i)), double(1))
}
