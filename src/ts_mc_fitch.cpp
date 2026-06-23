// Fast Monte Carlo Fitch scoring for single characters.
// Uses random_tree() from build_postorder (compiled as C).

#include <Rcpp.h>
#include <vector>
#include <cstdint>

extern "C" {
  void random_tree(int *parent_of, int *left, int *right, const int *n_tip);
  void seed_random_tree(unsigned long zs, unsigned long ws);
}

using namespace Rcpp;

//' Monte Carlo Fitch scores for a single character
//'
//' Generates `n_mc` random trees and scores each with a Fitch parsimony
//' downpass for a single character defined by `state_counts`.
//' Tree generation and scoring are done entirely in C with no R object
//' allocation per tree, making this very fast (~0.01 ms per tree).
//'
//' @param state_counts Integer vector giving the number of tips in each
//'   state.  Length determines the number of states (k); sum determines
//'   the number of tips (n).  For example, `c(13, 13, 12)` defines a
//'   3-state character with 38 tips.
//' @param n_mc Number of random trees to generate and score.
//' @return Integer vector of length `n_mc` containing the Fitch parsimony
//'   score (number of state changes) for each random tree.
//' @keywords internal
//' @export
// [[Rcpp::export]]
IntegerVector mc_fitch_scores(IntegerVector state_counts, int n_mc) {
  int k = state_counts.size();
  int n = 0;
  for (int i = 0; i < k; i++) n += state_counts[i];
  if (n < 2) return IntegerVector(n_mc, 0);

  // Seed the random-tree MWC generator from R's RNG so that profile parsimony
  // (whose Monte Carlo information estimate scores random trees via this path)
  // is reproducible under set.seed().  The Rcpp wrapper has already established
  // an RNGScope, so R's RNG state is loaded; draw two non-zero 32-bit seeds.
  unsigned long zs = 1UL + static_cast<unsigned long>(R::unif_rand() * 4294967294.0);
  unsigned long ws = 1UL + static_cast<unsigned long>(R::unif_rand() * 4294967294.0);
  seed_random_tree(zs, ws);

  // Build tip state bitmasks: one bit per state
  std::vector<uint32_t> tip_state(n);
  int idx = 0;
  for (int s = 0; s < k; s++) {
    for (int j = 0; j < state_counts[s]; j++) {
      tip_state[idx++] = (1u << s);
    }
  }

  int n_internal = n - 1;
  int n_node = 2 * n - 1;
  std::vector<int> parent(n_node), left(n_internal), right(n_internal);
  std::vector<uint32_t> state(n_node);
  std::vector<int> preorder;
  preorder.reserve(n_internal);
  std::vector<int> stack;
  stack.reserve(n_internal);

  IntegerVector scores(n_mc);
  for (int rep = 0; rep < n_mc; rep++) {
    random_tree(parent.data(), left.data(), right.data(), &n);

    // Build postorder: collect preorder via DFS, then process in reverse.
    preorder.clear();
    stack.clear();
    stack.push_back(n);  // root = n_tip
    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();
      if (node < n) continue;
      preorder.push_back(node);
      int ni = node - n;
      stack.push_back(right[ni]);
      stack.push_back(left[ni]);
    }

    // Fitch downpass in postorder (reverse of preorder)
    int score = 0;
    for (int t = 0; t < n; t++) state[t] = tip_state[t];

    for (int i = static_cast<int>(preorder.size()) - 1; i >= 0; --i) {
      int node = preorder[i];
      int ni = node - n;
      uint32_t ls = state[left[ni]];
      uint32_t rs = state[right[ni]];
      uint32_t inter = ls & rs;
      if (inter) {
        state[node] = inter;
      } else {
        state[node] = ls | rs;
        ++score;
      }
    }

    scores[rep] = score;
  }
  return scores;
}
