#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <Rcpp.h>
using namespace Rcpp;

namespace {

constexpr int MAX_FACTORIAL_LOOKUP = 8192;
constexpr double LOG2_E = 1.4426950408889634;

// Block-scope static: C++17 guarantees the initialization runs exactly once
// even if several threads reach it together.
const std::vector<double>& log2_factorial_table() {
  static const std::vector<double> table = []() {
    std::vector<double> t(MAX_FACTORIAL_LOOKUP + 1);
    t[0] = 0.0;
    for (int i = 1; i <= MAX_FACTORIAL_LOOKUP; ++i) {
      t[i] = t[i - 1] + std::log2(i);
    }
    return t;
  }();
  return table;
}

// Fast lookup with bounds checking
inline double l2factorial(int n) {
  if (n < 0) {
    Rcpp::stop("Factorial undefined for negative arguments.");
  }
  if (n <= MAX_FACTORIAL_LOOKUP) {
    return log2_factorial_table()[n];
  } else {
    return lgamma(n + 1) * LOG2_E;
  }
}

} // namespace

//' Expected mutual information between two partitions
//'
//' Computes the mutual information expected purely by chance between two
//' partitions of the same `N` items, under the hypergeometric null in which the
//' block sizes (marginals) of each partition are fixed but the items are
//' associated at random.  Subtracting this baseline from an observed mutual
//' information yields a chance-corrected ("adjusted") mutual information, as
//' applied by [`SiteConcordance`]`(normalize = TRUE)`; it is most material for
//' small or unbalanced partitions, where raw mutual information is appreciably
//' inflated by chance agreement.
//'
//' The value is computed analytically \insertCite{@Vinh2010}{TreeDist}, summing over
//' the hypergeometric distribution of cell overlaps, and is returned in bits
//' (logarithms to base two).
//'
//' @param ni Integer vector of length two giving the sizes of the two blocks of
//'   the first (bi-)partition; these sum to the total item count `N`.
//' @param nj Integer vector giving the block sizes of the second partition
//'   (also summing to `N`).
//' @return The expected mutual information, in bits.
//' @seealso [`SiteConcordance`]
//' @examples
//' # Expected MI between a 3|4 split and a 2|5 split of 7 items:
//' expected_mi(c(3L, 4L), c(2L, 5L))
//' @export
// [[Rcpp::export]]
double expected_mi(const IntegerVector &ni, const IntegerVector &nj) {
  if (ni.size() != 2) {
    Rcpp::stop("ni must be a vector of length 2.");
  }
  // ni and nj are vectors listing the number of entitites in each cluster
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

    // Mutual information contributed by an overlap of k, per unit probability
    const auto cell_mi = [&](int k) {
      double contribution = 0.0;
      // contribution from inside the split
      if (k > 0) {
        double mi_in = std::log2(static_cast<double>(k)) + log2N - (log2a + log2mj);
        contribution += (static_cast<double>(k) * invN) * mi_in;
      }
      // contribution from outside the split
      const int kout = mj - k;
      if (kout > 0) {
        double mi_out = std::log2(static_cast<double>(kout)) + log2N - (log2Na + log2mj);
        contribution += (static_cast<double>(kout) * invN) * mi_out;
      }
      return contribution;
    };

    // Anchor the recurrence at the mode of the hypergeometric.  P(K = kmode)
    // is the largest of at most N + 1 probabilities summing to one, so it is
    // always representable; P(K = kmin) is not — at N = 1200 it is around
    // 2^-1197, and a recurrence seeded with the zero it underflows to stays
    // zero for every remaining k.
    const int kmode = std::min(kmax, std::max(kmin, static_cast<int>(
      (static_cast<double>(mj) + 1.0) * (static_cast<double>(a) + 1.0) /
        (static_cast<double>(N) + 2.0))));

    const double log2Pmode =
      (l2factorial(mj) - l2factorial(kmode) - l2factorial(mj - kmode))
      + (l2factorial(N - mj) - l2factorial(a - kmode)
           - l2factorial(N - mj - (a - kmode)))
      - log2_denom;
    const double Pmode = std::exp2(log2Pmode);

    emi += cell_mi(kmode) * Pmode;

    // Walk down: P(k - 1) = P(k) * k(N - mj - a + k) / ((mj - k + 1)(a - k + 1))
    double Pk = Pmode;
    for (int k = kmode; k > kmin; --k) {
      Pk *= (static_cast<double>(k) * (N - mj - a + k)) /
        (static_cast<double>(mj - k + 1) * (a - k + 1));
      // The distribution is unimodal, so once the tail underflows every
      // remaining term is likewise negligible.
      if (!(Pk > 0.0)) break;
      emi += cell_mi(k - 1) * Pk;
    }

    // Walk up: P(k + 1) = P(k) * (mj - k)(a - k) / ((k + 1)(N - mj - a + k + 1))
    Pk = Pmode;
    for (int k = kmode; k < kmax; ++k) {
      Pk *= (static_cast<double>(mj - k) * (a - k)) /
        (static_cast<double>(k + 1) * (N - mj - a + k + 1));
      if (!(Pk > 0.0)) break;
      emi += cell_mi(k + 1) * Pk;
    }
  }

  return emi;
}

// [[Rcpp::export]]
std::string mi_key(IntegerVector ni, IntegerVector nj) {
  if (ni.size() != 2) {
    Rcpp::stop("ni must be a vector of length 2.");
  }
  
  std::vector<uint16_t> ni_vals = {static_cast<uint16_t>(ni[0]),
                                   static_cast<uint16_t>(ni[1])};
  std::sort(ni_vals.begin(), ni_vals.end());
  
  std::vector<uint16_t> nj_vals;
  nj_vals.reserve(nj.size());
  for (int val : nj) {
    nj_vals.push_back(static_cast<uint16_t>(val));
  }
  std::sort(nj_vals.begin(), nj_vals.end());
  
  // Encode each uint16_t as 4 hex characters — no R allocation needed
  static const char hex[] = "0123456789abcdef";
  std::string key;
  key.reserve((2 + nj_vals.size()) * 4);
  
  for (uint16_t v : ni_vals) {
    key += hex[(v >> 12) & 0xF];
    key += hex[(v >> 8)  & 0xF];
    key += hex[(v >> 4)  & 0xF];
    key += hex[(v)       & 0xF];
  }
  for (uint16_t v : nj_vals) {
    key += hex[(v >> 12) & 0xF];
    key += hex[(v >> 8)  & 0xF];
    key += hex[(v >> 4)  & 0xF];
    key += hex[(v)       & 0xF];
  }
  
  return key;
}
