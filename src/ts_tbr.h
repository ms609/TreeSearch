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
#include "ts_constraint.h"
#include "ts_pool.h"
#include <functional>
#include <vector>

namespace ts {

// Clip ordering strategy for TBR search.
enum class ClipOrder {
  RANDOM = 0,     // Current default: uniform random shuffle
  INV_WEIGHT = 1, // Weighted random, w = 1/(1+s) where s = subtree size
  TIPS_FIRST = 2, // All tip clips first (shuffled), then rest (shuffled)
  BUCKET = 3,     // Three size buckets: tips/small/large, random within each
  ANTI_TIP = 4,   // Non-tip clips (shuffled) first, tip clips (shuffled) last
  LARGE_FIRST = 5 // Large (>√n) first, then small (2..√n), then tips; random within each
};

struct TBRParams {
  bool accept_equal = false;     // accept Δ=0 moves?
  int max_accepted_changes = 0;  // 0 = no limit (run to convergence)
  int max_hits = 1;              // equal-score hits before stopping
  int tabu_size = 0;             // tabu list capacity (0 = disabled)
  ClipOrder clip_order = ClipOrder::RANDOM;
  bool diagnostics = false;      // collect per-pass diagnostic records
  // True unrooted TBR (default on).  At apparent convergence the kernel checks
  // the one root edge and (for NA) performs an exact full-neighbourhood sweep,
  // ensuring the result is a genuine unrooted-TBR optimum.  Gated out when
  // sector_mask / cd / tabu / pool are active (state would be invalidated).
  // See dev/plans/2026-06-18-tbr-shared-start.md.
  bool unrooted = true;
};

// Per-pass diagnostic record (populated only when TBRParams::diagnostics
// is true). One record per pass of the outer while loop.
struct TBRPassRecord {
  int pass_index;
  bool productive;              // true if a move was accepted
  int accepted_clip_size;       // subtree size of accepted clip (0 if null)
  int n_clips_tried;            // clips evaluated before acceptance (or total)
  int n_candidates_evaluated;   // total regraft×reroot evaluations this pass
};

struct TBRResult {
  double best_score;  // double for forward-compatibility with implied weights
  int n_accepted;
  int n_evaluated;
  int n_zero_skipped; // clips skipped due to zero-length edge (opt #7)
  bool converged;  // true if stopped due to no improvement
  std::vector<TBRPassRecord> pass_records; // populated when diagnostics=true
};

// Run TBR hill-climbing search on `tree` with dataset `ds`.
// Modifies `tree` in place to the best tree found.
// If `cd` is non-null and active, constraint-violating moves are skipped.
// If `sector_mask` is non-null, only clips and regrafts within the sector
// are considered (CSS = Constrained Sectorial Search).
// If `check_timeout` is non-null, it is polled periodically (every n_tip
// clips) and the search returns early if it returns true.
TBRResult tbr_search(TreeState& tree, const DataSet& ds,
                     const TBRParams& params,
                     ConstraintData* cd = nullptr,
                     const std::vector<bool>* sector_mask = nullptr,
                     TreePool* collect_pool = nullptr,
                     std::function<bool()> check_timeout = nullptr);

// Cache key used by exact_verify_sweep's optimum memoization (the NA
// convergence certifier).  A FALSE ("genuine optimum") verdict for a topology
// is valid only under the current (topology, dataset, weighting-regime) triple,
// so all three are mixed in.  The ratchet mutates the weighting regime
// (active_mask / upweight_mask / pattern_freq) in place mid-search, so the
// regime MUST be in the key or a base-regime verdict leaks into a perturbed
// pass and silently skips improvers.  Exposed so the regression test
// (test-ts-na-evcache.R) can assert the key is sensitive to each regime field
// AND to topology, using the exact code the cache uses.
uint64_t exact_verify_cache_key(const TreeState& tree, const DataSet& ds);

// Re-root `tree` so tip `t` becomes a direct child of the root pseudo-node.
// Parsimony length is root-invariant, so this only changes the representation
// (which edges sit adjacent to the root).  Rebuilds postorder; does NOT refresh
// the Fitch state arrays, so the caller must full_rescore() / score_tree()
// before reading scores or state-dependent quantities (e.g. collapse flags).
// Rooting on a TIP makes the root-adjacent edges trivial splits, which is what
// makes a subsequent collapse rooting-invariant (the criterion otherwise skips
// root-adjacent edges).
void reroot_at_tip(TreeState& tree, int t);

} // namespace ts

#endif // TS_TBR_H
