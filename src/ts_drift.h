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
#include "ts_constraint.h"
#include <functional>
#include <vector>

namespace ts {

struct DriftParams {
  int n_cycles = 10;
  int afd_limit = 3;       // max absolute fit difference (steps)
  double rfd_limit = 0.1;  // max relative fit difference
  int max_hits = 1;        // for search-phase TBR
  int tabu_size = 0;       // Tabu list size for TBR calls (0 = disabled)
};

struct DriftResult {
  double best_score;
  int n_cycles_completed;
  int total_drift_moves;   // suboptimal moves accepted during drift phases
  int total_tbr_moves;     // moves accepted during search phases
};

// Run drift search on `tree` with dataset `ds`.
// Modifies `tree` in place to the best tree found across all cycles.
//
// `sector_mask` (optional): when non-null, restricts every internal TBR clip and
// regraft to nodes flagged true, exactly as tbr_search's sector_mask does. Used
// by the sectorial search to PIN an HTU pseudo-tip: masking a sector's content
// clade keeps drift's suboptimal + search phases from re-rooting the reduced tree
// against the rest-of-tree summary, so the result stays reinsertable. nullptr
// (default) = drift the whole tree, prior behaviour.
DriftResult drift_search(TreeState& tree, const DataSet& ds,
                         const DriftParams& params,
                         ConstraintData* cd = nullptr,
                         std::function<bool()> check_timeout = nullptr,
                         const std::vector<bool>* sector_mask = nullptr);

} // namespace ts

#endif // TS_DRIFT_H
