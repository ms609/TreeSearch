#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <cstdint> 

using namespace Rcpp;

// --- CONSTANTS ---
const double INV_LN2 = 1.0 / M_LN2; 
const int CACHE_SIZE = 50000;

// --- STATIC STORAGE ---
static std::vector<double> memoized_log2_df;
static std::unordered_map<uint64_t, double> memoized_results;

// --- HELPERS ---

// Portable integer log2 (Handles MSVC vs GCC/Clang if needed, though Rtools uses GCC)
inline int ilog2_u64(uint64_t x) {
#if defined(_MSC_VER)
  unsigned long index;
  _BitScanReverse64(&index, x);
  return index;
#else
  return 63 - __builtin_clzll(x);
#endif
}

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
  // Fallback
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

// --- CORE LOGIC (Scalar, Cached) ---
// This function does ONLY the calculation. No Rcpp loops.
inline double calculate_core(int m, int a, int b) {
  int min_ab = (a < b) ? a : b;
  int max_ab = (a < b) ? b : a;
  
  uint64_t key = ((uint64_t)m << 42) | ((uint64_t)min_ab << 21) | (uint64_t)max_ab;
  
  auto it = memoized_results.find(key);
  if (it != memoized_results.end()) return it->second;
  
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
  
  if (std::abs(ret) < 1e-14) ret = 0.0;
  
  memoized_results[key] = ret;
  return ret;
}

// --- TEMPLATE DISPATCHER ---
// 0 = Log2, 1 = Natural Log, 2 = Raw (2^x)
template <int MODE>
NumericVector CarterRunner(IntegerVector m, IntegerVector a, IntegerVector b) {
  int n_size = m.size();
  NumericVector results(n_size);
  
  for(int i = 0; i < n_size; ++i) {
    double log2_val = calculate_core(m[i], a[i], b[i]);
    
    if (MODE == 0) {
      results[i] = log2_val;
    } else if (MODE == 1) {
      results[i] = log2_val * M_LN2; 
    } else if (MODE == 2) {
      results[i] = std::exp2(log2_val); // More precise than pow(2, val)
    }
  }
  return results;
}

// --- EXPORTED FUNCTIONS ---

// [[Rcpp::export]]
void ClearCarterCache() {
  memoized_results.clear();
}

// [[Rcpp::export]]
int CarterCacheSize() {
  return memoized_results.size();
}

// [[Rcpp::export]]
NumericVector Log2Carter1(IntegerVector m, IntegerVector a, IntegerVector b) {
  return CarterRunner<0>(m, a, b);
}

// [[Rcpp::export]]
NumericVector LogCarter1(IntegerVector m, IntegerVector a, IntegerVector b) {
  return CarterRunner<1>(m, a, b);
}


//' Number of trees with _m_ steps
//' 
//' Functions to compute the number of trees with a given parsimony score.
//' 
//' `Carter1()` calculates the number of trees in which Fitch parsimony will
//' reconstruct  _m_ steps, where _a_ leaves are labelled with one state,
//' and _b_ leaves are labelled with a second state, using theorem 1 of
//' \insertCite{Carter1990;textual}{TreeTools}
//' 
//' `MaddisonSlatkin()` generalises this result to characters with multiple
//' steps using the recursive approach of
//' \insertCite{Maddison1991;textual}{TreeSearch}.
//' 
//' @param m,steps Number of steps.
//' @param a,b Number of leaves labelled `0` and `1`.
//' @param states Number of leaves labelled with each possible state.
//' States are presented in binary fashion.  The first entry of the vector
//' corresponds to state `1` (binary `001`),
//' the second to state `2` (binary `010`),
//' and the third to the ambiguous state `01` (binary `011`).
//' 
//' 
//' @seealso [TreeTools::NUnrooted()]: number of unrooted trees with _n_ leaves.
//' @references 
//' \insertAllCited{}
//' 
//' See also:
//' 
//' \insertRef{Steel1993}{TreeSearch}
//' 
//' \insertRef{Steel1995}{TreeSearch}
//' 
//' (\insertRef{Steel1996}{TreeSearch})
//' @importFrom TreeTools LogDoubleFactorial
//' @examples 
//' # The character `0 0 0 1 1 1`
//' Carter1(1, 3, 3) # Exactly one step
//' Carter1(2, 3, 3) # Two steps (one extra step)
//' 
//' # Number of trees that the character can map onto with exactly _m_ steps
//' # if non-parsimonious reconstructions are permitted:
//' cumsum(sapply(1:3, Carter1, 3, 3))
//' 
//' # Three steps allow the character to map onto any of the 105 six-leaf trees.
//' Carter1(3, 3, 3)
//' NUnrooted(3 + 3)
//' 
//' @template MRS
//' @family profile parsimony functions
//' @export
// [[Rcpp::export]]
NumericVector Carter1(IntegerVector m, IntegerVector a, IntegerVector b) {
  return CarterRunner<2>(m, a, b);
}
