#ifndef TS_TBR_H
#define TS_TBR_H

// TBR (Tree Bisection and Reconnection) search.
//
// Extends SPR by trying all rerootings of the clipped subtree before
// regrafting. Uses indirect length calculation for fast candidate
// evaluation, with full rescore verification on the best candidate.
//
// The search loop is parameterized via TBRParams to support future
// ratchet and drifting extensions without refactoring.

#include "ts_data.h"
#include "ts_tree.h"

namespace ts {

struct TBRParams {
  bool accept_equal = false;     // accept Δ=0 moves?
  int max_accepted_changes = 0;  // 0 = no limit (run to convergence)
  int max_hits = 1;              // equal-score hits before stopping
};

struct TBRResult {
  double best_score;  // double for forward-compatibility with implied weights
  int n_accepted;
  int n_evaluated;
  bool converged;  // true if stopped due to no improvement
};

// Run TBR hill-climbing search on `tree` with dataset `ds`.
// Modifies `tree` in place to the best tree found.
TBRResult tbr_search(TreeState& tree, const DataSet& ds,
                     const TBRParams& params);

} // namespace ts

#endif // TS_TBR_H
