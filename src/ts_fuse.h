#ifndef TS_FUSE_H
#define TS_FUSE_H

// Tree fusing: combine the best parts of multiple suboptimal trees
// by exchanging shared clades (bipartitions) between a "recipient"
// and "donor" trees from a pool.
//
// Algorithm (Goloboff 1999):
//   1. Start with best tree from pool as recipient.
//   2. For each donor in pool, find shared splits (bipartitions).
//   3. Try exchanging shared clades bottom-up (smallest first).
//   4. If any exchange improves (or equals, if accept_equal) the score,
//      apply it and skip ancestor splits of the exchanged clade.
//   5. After any improvement round, run TBR to clean up.
//   6. Repeat until no improvement found or max_rounds reached.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_pool.h"

namespace ts {

struct FuseParams {
  bool accept_equal = false;   // accept equal-score exchanges?
  int max_rounds = 10;         // max improvement rounds
};

struct FuseResult {
  double best_score;
  int n_exchanges;             // number of exchanges applied
  int n_rounds;                // number of improvement rounds
};

// Fuse trees from the pool. Modifies `recipient` in place.
// Runs TBR after each round of improvements.
FuseResult tree_fuse(TreeState& recipient, const DataSet& ds,
                     const TreePool& pool, const FuseParams& params);

} // namespace ts

#endif // TS_FUSE_H
