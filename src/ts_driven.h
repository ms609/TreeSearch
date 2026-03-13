#ifndef TS_DRIVEN_H
#define TS_DRIVEN_H

// Driven search: orchestrates Wagner start, TBR, ratchet, sectorial
// search (XSS), and tree fusing into a single search strategy.
//
// Based on TNT's xmult / combosearch model (Goloboff 1999; Goloboff
// & Pol 2007; Goloboff Ch. 5).

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_pool.h"

namespace ts {

struct DrivenParams {
  int max_replicates = 20;       // max RAS+search replicates
  int target_hits = 5;           // stop after N independent hits to best

  // TBR: equal-score plateau exploration (requires accept_equal = true)
  int tbr_max_hits = 20;

  // Ratchet
  int ratchet_cycles = 10;
  double ratchet_perturb_prob = 0.04;

  // Sectorial search (XSS only; RSS not currently wired)
  int xss_rounds = 3;
  int xss_partitions = 4;
  int sector_min_size = 6;
  int sector_max_size = 50;

  // Tree fusing
  int fuse_interval = 3;         // fuse every N replicates
  bool fuse_accept_equal = false;

  // Pool
  int pool_max_size = 100;
  double pool_suboptimal = 0.0;  // 0 = keep only optimal
};

struct DrivenResult {
  double best_score;
  int replicates_completed;
  int hits_to_best;
  int pool_size;
};

// Run the full driven search. Returns the best tree found (written
// into `best_tree`) and search statistics.
DrivenResult driven_search(TreeState& best_tree, DataSet& ds,
                           const DrivenParams& params);

} // namespace ts

#endif // TS_DRIVEN_H
