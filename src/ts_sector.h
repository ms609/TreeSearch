#ifndef TS_SECTOR_H
#define TS_SECTOR_H

// Sectorial search: divide the tree into sectors (clades), solve each
// sub-problem independently on a reduced dataset, and reinsert improved
// resolutions. Supports Random Sectorial Search (RSS) and Exclusive
// Sectorial Search (XSS).
//
// Based on Goloboff 1999; Goloboff & Pol 2007.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_constraint.h"
#include <functional>
#include <vector>

namespace ts {

struct SplitFrequencyTable;  // forward declaration (defined in ts_pool.h)

struct SectorParams {
  int min_sector_size = 6;       // minimum tips in a sector (RSS)
  int max_sector_size = 50;      // maximum tips in a sector (RSS)
  int n_partitions = 4;          // number of partitions (XSS)
  bool accept_equal = false;     // accept equal-score sector resolutions
  int rss_picks_per_round = 0;   // 0 = auto (2 * n_tip / avg_sector_size)
  int xss_rounds = 3;            // number of XSS rounds
  int internal_ratchet_cycles = 6;
  int internal_max_hits = 1;     // max_hits for internal TBR search
  int ras_starts = 1;            // RAS+TBR restarts per sector (start 0 = TBR on
                                 // the existing subtree; 1 = prior behaviour;
                                 // TNT uses 3 random-addition restarts)
  int collapse_target = 0;       // >0: collapse a selected clade's deep sub-clades
                                 // into ~this many composite first-pass terminals
                                 // (Goloboff 1999 coarse-grained sector). 0 = off.

  // Conflict-guided sector selection.
  // When non-null, RSS uses weighted random selection that biases toward
  // sectors containing splits absent from the pool consensus.
  // Owned externally (by driven_search); never freed by sector code.
  const SplitFrequencyTable* split_freq = nullptr;

  int clip_order = 0;  // ClipOrder cast to int (RANDOM = 0)
};

struct SectorResult {
  double best_score;
  int n_sectors_searched;
  int n_sectors_improved;
  int total_steps_saved;
};

// Reduced dataset for a single sector.
struct ReducedDataset {
  DataSet data;                    // bit-packed, same block structure
  TreeState subtree;               // topology for the sector
  int sector_root;                 // node index in full tree
  int n_real_tips;                 // leaf descendants (real OTUs)
  int n_htus;                      // HTU pseudo-tips (typically 1)
  std::vector<int> sector_to_full; // sector node → full tree node
  std::vector<int> full_to_sector; // full tree node → sector node (-1 if absent)
};

// Build a reduced dataset for the clade rooted at `sector_root`.
// Requires that `tree` has current `final_` states (run fitch_score first).
ReducedDataset build_reduced_dataset(const TreeState& tree,
                                     const DataSet& ds,
                                     int sector_root);

// Count the number of leaf descendants of `node`.
int count_clade_tips(const TreeState& tree, int node);

// Random Sectorial Search: pick random sectors, search, reinsert.
// Modifies `tree` in place.
SectorResult rss_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params,
                        ConstraintData* cd = nullptr,
                        std::function<bool()> check_timeout = nullptr);

// Exclusive Sectorial Search: partition tree into non-overlapping sectors.
// Modifies `tree` in place.
SectorResult xss_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params,
                        ConstraintData* cd = nullptr,
                        std::function<bool()> check_timeout = nullptr);

// Constrained Sectorial Search: sector-restricted TBR on the full tree.
// Unlike RSS/XSS, does not build a reduced dataset with an HTU pseudo-tip.
// Instead, restricts TBR clips and regrafts to within each sector, scoring
// against the full tree for exact evaluation. Eliminates HTU approximation
// errors at the cost of higher per-candidate evaluation.
SectorResult css_search(TreeState& tree, DataSet& ds,
                        const SectorParams& params,
                        ConstraintData* cd = nullptr,
                        std::function<bool()> check_timeout = nullptr);

} // namespace ts

#endif // TS_SECTOR_H
