#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <cstdint> // For uint64_t

using namespace Rcpp;

const double LN2 = M_LN2;
const double INV_LN2 = 1.0 / std::log(2.0);
const int CACHE_SIZE = 50000;

// --- STATIC STORAGE ---
// 1. The Factorial Cache (Vector)
static std::vector<double> memoized_log2_df;

// 2. The Result Memoization (Hash Map)
// Key: uint64_t packed representation of (m, min(a,b), max(a,b))
// Value: The calculated result
static std::unordered_map<uint64_t, double> memoized_results;

// --- HELPERS ---

inline int ilog2_u64(uint64_t x) {
  return 63 - __builtin_clzll(x);
}

// Internal helper to compute the factorial table once
void build_factorial_cache() {
  memoized_log2_df.reserve(CACHE_SIZE + 1);
  memoized_log2_df.push_back(0.0); 
  memoized_log2_df.push_back(0.0); 
  
  for (int n = 2; n <= CACHE_SIZE; ++n) {
    double val;
    if (n % 2 != 0) { 
      double oddN = (double)n;
      double nPlusOneOverTwo = (oddN + 1.0) / 2.0;
      val = ((std::lgamma(oddN + 1.0) - std::lgamma(nPlusOneOverTwo)) * INV_LN2) - nPlusOneOverTwo + 1.0;
    } else {
      double evenN = (double)n;
      double nOverTwo = evenN / 2.0;
      val = ilog2_u64(evenN) + (std::lgamma(nOverTwo) * INV_LN2) + (nOverTwo - 1.0);
    }
    memoized_log2_df.push_back(val);
  }
}

inline double get_log2_df(int n) {
  if (n < 2) return 0.0; 
  if (n <= CACHE_SIZE) {
    if (memoized_log2_df.empty()) build_factorial_cache();
    return memoized_log2_df[n];
  }
  // Fallback logic for huge n (omitted for brevity, same as previous step)
  if (n & 1) {
    double oddN = (double)n;
    double nPlusOneOverTwo = (oddN + 1.0) / 2.0;
    return ((std::lgamma(oddN + 1.0) - std::lgamma(nPlusOneOverTwo)) * INV_LN2) - nPlusOneOverTwo + 1.0;
  } else {
    double evenN = (double)n;
    double nOverTwo = evenN / 2.0;
    return ilog2_u64(evenN) + (std::lgamma(nOverTwo) * INV_LN2) + (nOverTwo - 1.0);
  }
}

inline double log2_n_calc(int n, int m) {
  if (n < m) return -std::numeric_limits<double>::infinity();
  double n_minus_m = n - m;
  return (std::lgamma(n + n_minus_m) - std::lgamma(n_minus_m + 1.0) - std::lgamma(m)) * INV_LN2 - n_minus_m;
}

// --- EXPORTED FUNCTIONS ---
// [[Rcpp::export]]
void ClearCarterCache() {
  memoized_results.clear();
  // We generally don't clear the factorial cache as it's small and constant
}

// [[Rcpp::export]]
int CarterCacheSize() {
  return memoized_results.size();
}

// [[Rcpp::export]]
double Log2Carter1_cpp(int m, int a, int b) {
  // 1. Symmetry Optimization: a and b are interchangeable in the formula
  int min_ab, max_ab;
  if (a < b) { min_ab = a; max_ab = b; } else { min_ab = b; max_ab = a; }
  
  // 2. Generate Key
  // We assume inputs < 2,097,151 (21 bits). 
  // Layout: [m (22 bits)] [min_ab (21 bits)] [max_ab (21 bits)]
  // This fits in a 64-bit integer.
  uint64_t key = ((uint64_t)m << 42) | ((uint64_t)min_ab << 21) | (uint64_t)max_ab;
  
  // 3. Check Memoization Map
  auto it = memoized_results.find(key);
  if (it != memoized_results.end()) {
    return it->second;
  }
  
  // 4. Perform Calculation (if not found)
  int n = a + b;
  int twoN = n + n;
  int twoM = m + m;
  
  double term1 = std::log2(twoN - twoM - m);
  double term2 = std::lgamma(m) * INV_LN2; 
  double term3 = get_log2_df(twoN - 5);
  double term4 = log2_n_calc(a, m);
  double term5 = log2_n_calc(b, m);
  double sub   = get_log2_df(twoN - twoM - 1);
  
  double ret = (term1 + term2 + term3 + term4 + term5) - sub;
  
  if (std::abs(ret) < std::sqrt(std::numeric_limits<double>::epsilon())) {
    ret = 0.0;
  }
  
  // 5. Store Result
  memoized_results[key] = ret;
  
  return ret;
}
