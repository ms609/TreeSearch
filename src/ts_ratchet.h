#ifndef TS_RATCHET_H
#define TS_RATCHET_H

// Parsimony ratchet: escape local optima by temporarily zeroing random
// characters (via active_mask perturbation), then searching the perturbed
// landscape with TBR before reverting to original weights.

#include "ts_data.h"
#include "ts_tree.h"

namespace ts {

struct RatchetParams {
  int n_cycles = 10;           // number of perturbation+search cycles
  double perturb_prob = 0.04;  // probability of zeroing each character
  int max_hits = 1;            // for the search-phase TBR
};

struct RatchetResult {
  double best_score;
  int n_cycles_completed;
  int total_tbr_moves;         // across all cycles
};

RatchetResult ratchet_search(TreeState& tree, DataSet& ds,
                             const RatchetParams& params);

} // namespace ts

#endif // TS_RATCHET_H
