#ifndef TS_SEARCH_H
#define TS_SEARCH_H

// NNI hill-climbing search.
//
// First-improvement strategy: iterate over NNI edges in random order,
// try both swap variants, accept immediately if score improves (or equals,
// up to maxHits). Repeat passes until no improvement found.

#include "ts_data.h"
#include "ts_tree.h"

namespace ts {

struct SearchResult {
  int score;           // best score found
  int n_moves;         // total improving moves accepted
  int n_iterations;    // total NNI candidates evaluated
};

// Run NNI hill-climbing search on `tree` with dataset `ds`.
// Modifies `tree` in place to the best tree found.
// `maxHits`: number of times the best score must be hit without improvement
//            before stopping (0 = stop on first pass with no improvement).
SearchResult nni_search(TreeState& tree, const DataSet& ds, int maxHits);

// Run SPR hill-climbing search using indirect calculation (Goloboff 1996).
// For each candidate clip node, uses incremental two-pass (Shortcut C) to
// compute divided-tree final states, then indirect length calculation for
// each destination edge. First-improvement with random clip order.
SearchResult spr_search(TreeState& tree, const DataSet& ds, int maxHits);

} // namespace ts

#endif // TS_SEARCH_H
