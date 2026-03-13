#include "ts_driven.h"
#include "ts_fitch.h"
#include "ts_tbr.h"
#include "ts_ratchet.h"
#include "ts_drift.h"
#include "ts_sector.h"
#include "ts_fuse.h"
#include "ts_wagner.h"
#include "ts_splits.h"

#include <R.h>

namespace ts {

DrivenResult driven_search(TreeState& best_tree, DataSet& ds,
                           const DrivenParams& params) {
  TreePool pool(params.pool_max_size, params.pool_suboptimal);

  DrivenResult result;
  result.best_score = 1e18;
  result.replicates_completed = 0;
  result.hits_to_best = 0;
  result.pool_size = 0;

  bool tree_large_enough_for_sectors =
      ds.n_tips >= 2 * params.sector_min_size;

  for (int rep = 0; rep < params.max_replicates; ++rep) {

    // 1. Random addition sequence Wagner tree
    TreeState tree;
    random_wagner_tree(tree, ds);

    // 2. TBR to local optimum
    {
      TBRParams tp;
      tbr_search(tree, ds, tp);
    }

    // 3. Sectorial search (XSS) if tree is large enough
    if (tree_large_enough_for_sectors) {
      SectorParams sp;
      sp.n_partitions = params.xss_partitions;
      sp.xss_rounds = params.xss_rounds;
      sp.min_sector_size = params.sector_min_size;
      sp.max_sector_size = params.sector_max_size;
      sp.internal_ratchet_cycles = 0;
      sp.internal_max_hits = 1;
      xss_search(tree, ds, sp);
    }

    // 4. Ratchet perturbation to escape local optima
    {
      RatchetParams rp;
      rp.n_cycles = params.ratchet_cycles;
      rp.perturb_prob = params.ratchet_perturb_prob;
      rp.max_hits = params.tbr_max_hits;
      ratchet_search(tree, ds, rp);
    }

    // 5. Drifting (suboptimal + equal-score exploration)
    if (params.drift_cycles > 0) {
      DriftParams dp;
      dp.n_cycles = params.drift_cycles;
      dp.afd_limit = params.drift_afd_limit;
      dp.rfd_limit = params.drift_rfd_limit;
      dp.max_hits = params.tbr_max_hits;
      drift_search(tree, ds, dp);
    }

    // 6. Final TBR polish
    {
      TBRParams tp;
      tbr_search(tree, ds, tp);
    }

    // 7. Score and add to pool
    double score = score_tree(tree, ds);
    pool.add(tree, score);

    ++result.replicates_completed;

    // 8. Periodic tree fusing
    // Fused trees are derived from pool entries, not independent discoveries.
    // Don't let them inflate the convergence counter (hits_to_best).
    if ((rep + 1) % params.fuse_interval == 0 && pool.size() >= 2) {
      int hits_before = pool.hits_to_best();
      double best_before = pool.best_score();

      TreeState fused = pool.best().tree;
      FuseParams fp;
      fp.accept_equal = params.fuse_accept_equal;
      fp.max_rounds = 10;
      tree_fuse(fused, ds, pool, fp);

      double fused_score = score_tree(fused, ds);
      pool.add(fused, fused_score);

      if (fused_score < best_before) {
        // Fusing found a new best — no independent hits to this new score yet
        pool.set_hits_to_best(0);
      } else {
        // Fusing matched or didn't improve — restore counter
        pool.set_hits_to_best(hits_before);
      }
    }

    // 9. Convergence check
    if (pool.hits_to_best() >= params.target_hits) {
      break;
    }

    R_CheckUserInterrupt();
  }

  // Return the best tree
  result.hits_to_best = pool.hits_to_best();
  result.pool_size = pool.size();

  if (pool.size() > 0) {
    result.best_score = pool.best_score();
    best_tree = pool.best().tree;
  } else {
    // No replicates completed (e.g. max_replicates=0)
    result.best_score = -1.0;
  }

  return result;
}

} // namespace ts
