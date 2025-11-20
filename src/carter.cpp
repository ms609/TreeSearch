#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

const double LN2 = M_LN2;
const int CACHE_SIZE = 50000;

// A static vector to hold the values. 
// Being 'static', it persists between function calls within the same R session.
static std::vector<double> memoized_log2_df;

// Internal helper to compute the table once
void build_cache() {
  memoized_log2_df.reserve(CACHE_SIZE + 1);
  
  // Handle 0 and 1 (both map to result 0.0)
  memoized_log2_df.push_back(0.0); // Index 0 (unused usually, but safe padding)
  memoized_log2_df.push_back(0.0); // Index 1 (n=1)
  
  // Loop from 2 to CACHE_SIZE
  for (int n = 2; n <= CACHE_SIZE; ++n) {
    double val;
    if (n % 2 != 0) { 
      // Odd Case
      double oddN = (double)n;
      double nPlusOneOverTwo = (oddN + 1.0) / 2.0;
      val = ((std::lgamma(oddN + 1.0) - std::lgamma(nPlusOneOverTwo)) / LN2) - nPlusOneOverTwo + 1.0;
    } else {
      // Even Case
      double evenN = (double)n;
      double nOverTwo = evenN / 2.0;
      val = std::log2(evenN) + (std::lgamma(nOverTwo) / LN2) + (nOverTwo - 1.0);
    }
    memoized_log2_df.push_back(val);
  }
}

// Helper: Log2DoubleFactorial
inline double get_log2_df(int n) {
  // Check bounds for safety check: n[n < 2L] <- 1 (which is 0 in log space)
  if (n < 2) return 0.0; 
  
  // If within cache range, return cached value
  if (n <= CACHE_SIZE) {
    // Initialize cache on first run
    if (memoized_log2_df.empty()) {
      build_cache();
    }
    // R uses 1-based indexing, but our n aligns with vector index 
    // because we pushed 0 and 1 at the start.
    return memoized_log2_df[n];
  }
  
  // Fallback for n > 50000 (Slow path, rarely hit)
  if (n % 2 != 0) {
    double oddN = (double)n;
    double nPlusOneOverTwo = (oddN + 1.0) / 2.0;
    return ((std::lgamma(oddN + 1.0) - std::lgamma(nPlusOneOverTwo)) / LN2) - nPlusOneOverTwo + 1.0;
  } else {
    double evenN = (double)n;
    double nOverTwo = evenN / 2.0;
    return std::log2(evenN) + (std::lgamma(nOverTwo) / LN2) + (nOverTwo - 1.0);
  }
}

// Helper: Log2N
inline double log2_n_calc(int n, int m) {
  if (n < m) return -std::numeric_limits<double>::infinity();
  
  double n_minus_m = n - m;
  return (std::lgamma(n + n_minus_m) - 
          std::lgamma(n_minus_m + 1.0) - 
          std::lgamma(m)) / LN2 - n_minus_m;
}

// [[Rcpp::export]]
double Log2Carter1_cpp(int m, int a, int b) {
  int n = a + b;
  int twoN = n + n;
  int twoM = m + m;
  
  double term1 = std::log2(twoN - twoM - m);
  double term2 = std::lgamma(m) / LN2; 
  double term3 = get_log2_df(twoN - 5);
  double term4 = log2_n_calc(a, m);
  double term5 = log2_n_calc(b, m);
  double sub   = get_log2_df(twoN - twoM - 1);
  
  double ret = (term1 + term2 + term3 + term4 + term5) - sub;
  
  if (std::abs(ret) < std::sqrt(std::numeric_limits<double>::epsilon())) {
    return 0.0;
  }
  
  return ret;
}
