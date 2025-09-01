#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <Rcpp.h>
using namespace Rcpp;

#define MAX_FACTORIAL_LOOKUP 8192
static double log2_factorial_table[MAX_FACTORIAL_LOOKUP + 1];
static const double LOG2_E = 1.4426950408889634;

__attribute__((constructor))
void initialize_factorial_cache() {
  log2_factorial_table[0] = 0.0;
  for (int i = 1; i <= MAX_FACTORIAL_LOOKUP; i++) {
    log2_factorial_table[i] = log2_factorial_table[i - 1] + std::log2(i);
  }
}

// Fast lookup with bounds checking
inline double l2factorial(int n) {
  if (n <= MAX_FACTORIAL_LOOKUP) {
    return log2_factorial_table[n];
  } else {
    return lgamma(n + 1) * LOG2_E;
  }
}

// ni and nj are vectors listing the number of entitites in each cluster
// [[Rcpp::export]]
double expected_mi(const IntegerVector &ni, const IntegerVector &nj) {
  // ni = {a, N-a}; nj = counts of character states
  const int a = ni[0];
  const int N = ni[0] + ni[1];
  if (a <= 0 || a >= N) return 0.0; // trivial split
  
  const double invN = 1.0 / static_cast<double>(N);
  const double log2N = std::log2(static_cast<double>(N));
  const double log2a  = std::log2(static_cast<double>(a));
  const double log2Na = std::log2(static_cast<double>(N - a));
  const double log2_denom = l2factorial(N) - l2factorial(a) - l2factorial(N - a);
  
  double emi = 0.0;
  
  for (int j = 0; j < nj.size(); ++j) {
    int mj = nj[j];
    if (mj <= 0) continue;
    
    int kmin = std::max(0, a + mj - N);
    int kmax = std::min(a, mj);
    if (kmin > kmax) continue;
    
    const double log2mj = std::log2(static_cast<double>(mj));
    
    // compute P(K=kmin)
    double log2P = (l2factorial(mj) - l2factorial(kmin) - l2factorial(mj - kmin))
      + (l2factorial(N - mj) - l2factorial(a - kmin) - l2factorial(N - mj - (a - kmin)))
      - log2_denom;
      double Pk = std::pow(2.0, log2P);
      
      for (int k = kmin; k <= kmax; ++k) {
        if (Pk > 0.0) {
          // contribution from inside the split
          if (k > 0) {
            double mi_in = std::log2(static_cast<double>(k)) + log2N - (log2a + log2mj);
            emi += (static_cast<double>(k) * invN) * mi_in * Pk;
          }
          // contribution from outside the split
          int kout = mj - k;
          if (kout > 0) {
            double mi_out = std::log2(static_cast<double>(kout)) + log2N - (log2Na + log2mj);
            emi += (static_cast<double>(kout) * invN) * mi_out * Pk;
          }
        }
        // Update P(k) → P(k+1)
        if (k < kmax) {
          double numer = static_cast<double>((mj - k) * (a - k));
          double denom = static_cast<double>((k + 1) * (N - mj - a + k + 1));
          Pk *= numer / denom;
        }
      }
  }
  
  return emi;
}

// [[Rcpp::export]]
RawVector mi_key(IntegerVector ni, IntegerVector nj) {
  if (ni.size() != 2) {
    Rcpp::stop("ni must be a vector of length 2.");
  }
  
  std::vector<uint16_t> ni_vals = {static_cast<uint16_t>(ni[0]),
                                   static_cast<uint16_t>(ni[1])};
  if (ni_vals[0] > 65535 || ni_vals[1] > 65535) {
    Rcpp::stop("ni values must be ≤ 65535.");
  }
  std::sort(ni_vals.begin(), ni_vals.end());
  
  std::vector<uint16_t> nj_vals;
  for (int val : nj) {
    if (val > 65535) {
      Rcpp::stop("nj values must be ≤ 65535.");
    }
    nj_vals.push_back(static_cast<uint16_t>(val));
  }
  std::sort(nj_vals.begin(), nj_vals.end());
  
  
  // Total number of 16-bit ints
  size_t n = 2 + nj_vals.size();
  RawVector key_raw(n * 2);
  
  // Write ni
  key_raw[0] = ni_vals[0] >> 8;
  key_raw[1] = ni_vals[0] & 0xFF;
  key_raw[2] = ni_vals[1] >> 8;
  key_raw[3] = ni_vals[1] & 0xFF;
  
  // Write nj
  for (size_t i = 0; i < nj_vals.size(); ++i) {
    key_raw[4 + i + i]     = nj_vals[i] >> 8;
    key_raw[4 + i + i + 1] = nj_vals[i] & 0xFF;
  }
  
  return key_raw;
}
