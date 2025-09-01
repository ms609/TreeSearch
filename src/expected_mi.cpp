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


// Originally sourced from expected_MI in aricode v1.0.3 under the GPL3 license.
// doi:10.32614/CRAN.package.aricode
// Optimized and converted to base 2 by Martin R. Smith 2025-07-03
// ni and nj are vectors listing the number of entitites in each cluster
// [[Rcpp::export]]
double expected_mi(IntegerVector ni, IntegerVector nj) {
  
  int N = sum(ni);
  double emi = 0.0;
  
  // Basic assertions that should always hold
  assert(ni.size() >= 2 && "ni_ must have at least 2 elements");
  assert(nj.size() >= 2 && "n_j must have at least 2 elements");
  assert(sum(nj) == sum(ni));
  
  // Since we know sizes are >= 2, we can use more aggressive optimizations
  int ni_size = ni.size();
  int nj_size = nj.size();
  
  // Pre-allocate and compute all factorials at once
  std::vector<double> ni_f(ni_size);
  std::vector<double> nj_f(nj_size);
  std::vector<double> Nmni_f(ni_size);
  std::vector<double> Nmnj_f(nj_size);
  std::vector<double> log2_ni(ni_size);
  std::vector<double> log2_nj(nj_size);
  
  for (int i = 0; i < ni_size; ++i) {
    ni_f[i] = l2factorial(ni[i]);
    Nmni_f[i] = l2factorial(N - ni[i]);
    log2_ni[i] = (ni[i] > 0) ? std::log2(ni[i]) : -INFINITY;
  }
  for (int j = 0; j < nj_size; ++j) {
    nj_f[j] = l2factorial(nj[j]);
    Nmnj_f[j] = l2factorial(N - nj[j]);
    log2_nj[j] = (nj[j] > 0) ? std::log2(nj[j]) : -INFINITY;
  }
  
  double N_f = l2factorial(N);
  double log2_N = std::log2(N);
  double inv_N = 1.0 / N;
  
  for (int i = 0; i < ni_size; ++i) {
    if (ni[i] == 0) continue;
    
    for (int j = 0; j < nj_size; ++j) {
      if (nj[j] == 0) continue;
      
      int start_nij = std::max(1, ni[i] + nj[j] - N);
      int end_nij = std::min(ni[i], nj[j]);
      
      if (start_nij > end_nij) continue;
      
      double log2_ni_nj = log2_ni[i] + log2_nj[j];
      double base_log_prob = ni_f[i] + nj_f[j] + Nmni_f[i] + Nmnj_f[j] - N_f;
      
      for (int nij = start_nij; nij <= end_nij; ++nij) {
        double log2_nij = std::log2(nij);
        double mi_term = (log2_nij + log2_N - log2_ni_nj);
        
        double log_prob = base_log_prob -
          l2factorial(nij) - 
          l2factorial(ni[i] - nij) -
          l2factorial(nj[j] - nij) -
          l2factorial(N - ni[i] - nj[j] + nij);
        
        double prob = std::pow(2.0, log_prob);
        
        emi += nij * inv_N * mi_term * prob;
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
