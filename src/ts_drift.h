#ifndef TS_DRIFT_H
#define TS_DRIFT_H

// Tree drifting: escape local optima by accepting suboptimal TBR moves.
//
// Uses Absolute Fit Difference (AFD) and Relative Fit Difference (RFD)
// criteria to control drift distance. Alternates suboptimal drift cycles
// with equal-score drift cycles, each followed by a standard TBR search
// phase. Based on the TNT drift algorithm (Goloboff 1999).

#include "ts_data.h"
#include "ts_tree.h"

namespace ts {

struct DriftParams {
  int n_cycles = 10;
  int afd_limit = 3;       // max absolute fit difference (steps)
  double rfd_limit = 0.1;  // max relative fit difference
  int max_hits = 1;        // for search-phase TBR
};

struct DriftResult {
  double best_score;
  int n_cycles_completed;
  int total_drift_moves;   // suboptimal moves accepted during drift phases
  int total_tbr_moves;     // moves accepted during search phases
};

// Run drift search on `tree` with dataset `ds`.
// Modifies `tree` in place to the best tree found across all cycles.
DriftResult drift_search(TreeState& tree, const DataSet& ds,
                         const DriftParams& params);

} // namespace ts

#endif // TS_DRIFT_H
