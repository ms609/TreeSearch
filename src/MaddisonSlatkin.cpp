#include <Rcpp.h>
#include <memory> /* for unique_ptr */
using namespace std;
using namespace Rcpp;

typedef int_fast16_t int16;
typedef int_fast32_t int32;
const int16 MAX_TIPS = 128;

double lg2[int32(MAX_TIPS - 1) * (MAX_TIPS - 1) + 1];
double lg2_double_factorial[MAX_TIPS + MAX_TIPS - 2];
double lg2_rooted[MAX_TIPS + 1];
double lg2_unrooted[MAX_TIPS + 1];
__attribute__((constructor))
  void initialize_ldf() {
    lg2[0] = 0;
    for (int32 i = 1; i != int32(MAX_TIPS - 1) * (MAX_TIPS - 1) + 1; i++) { 
      lg2[i] = log2(i);
    }
    for (int16 i = 0; i != 3; i++) {
      lg2_double_factorial[i] = 0;
      lg2_rooted[i] = 0;
      lg2_unrooted[i] = 0;
    }
    for (int16 i = 2; i != MAX_TIPS + MAX_TIPS - 2; i++) {
      lg2_double_factorial[i] = lg2_double_factorial[i - 2] + lg2[i];
    }
    for (int16 i = 3; i != MAX_TIPS + 1; i++) {
      lg2_unrooted[i] = lg2_double_factorial[i + i - 5];
      lg2_rooted[i] = lg2_double_factorial[i + i - 3];
    }
  }


/* https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c/53983114 */
double choose (const int n, int k) {
  if (k > n) return 0;
  if (k + k > n) k = n - k;
  if (k == 0) return 1;
  
  double result = n;
  for (int16 i = 2; i <= k; i++) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

/* D: The probability that, in a randomly selected tree of _n_ taxa, _i_ of */
/* which have state 1, the smaller subclade of _m_ taxa will receive _j_ taxa */
/* with state 1 */
double D (const int16 j, const int16 i, const int16 m, const int16 n) {
  return choose(i, j) * choose(n - i, m - j) / choose(n, m);
}


/* R: The probability that in a randomly selected tree of _n_ leaves, the 
 * smaller of the two basal subclades will have _m_ leaves. 
 * My function gives half the value given in Maddison & Slatkin 1991 [!] 
double R (const int16 m, const int16 n) {
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
#' @param n Number of leaves in tree
#' @param x Vector specifying number of leaves with state 1, 2, ...
#' @param b State reconstructed at base of tree
#' States are converted to binary, so 
#' `1` denotes 'state 1',
#' `2` denotes 'state 2', and
#' `3` denotes 'state 1 OR state 2'.
#' 
#' @references
#' \insertRef{Maddison1991}{TreeSearch}
#' @examples
#' MaddisonSlatkin(2, 4)
#' @template MRS
#' @importFrom R.cache addMemoization
#' @importFrom TreeTools NRooted
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
      */