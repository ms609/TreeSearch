#ifndef TS_DRIVEN_H
#define TS_DRIVEN_H

// Driven search: orchestrates Wagner start, TBR, ratchet, sectorial
// search (XSS), and tree fusing into a single search strategy.
//
// Based on TNT's xmult / combosearch model (Goloboff 1999; Goloboff
// & Pol 2007; Goloboff Ch. 5).
//
// Agent C additions: all-pool return, timeout, interrupt checks,
// verbosity.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_pool.h"
#include "ts_constraint.h"
#include "ts_strategy.h"
#include <functional>

namespace ts {

// Progress information passed to the callback at each reporting point.
struct ProgressInfo {
  int replicate;          // 1-based current replicate
  int max_replicates;     // configured maximum
  double best_score;      // pool's current best (1e18 if pool empty)
  int hits_to_best;       // independent discoveries of best
  int target_hits;        // convergence target
  int pool_size;          // trees in pool
  const char* phase;      // "wagner", "tbr", "xss", "rss", "ratchet",
                          // "drift", "final_tbr", "fuse", "replicate", "done"
  double elapsed_seconds; // wall-clock since search start
  double phase_score;     // score after this phase (-1 if N/A)
};

struct DrivenParams {
  int max_replicates = 20;       // max RAS+search replicates
  int target_hits = 5;           // stop after N independent hits to best

  // TBR
  int tbr_max_hits = 1;

  // Ratchet
  int ratchet_cycles = 5;
  double ratchet_perturb_prob = 0.04;
  int ratchet_perturb_mode = 0;       // 0=zero, 1=upweight, 2=mixed
  int ratchet_perturb_max_moves = 0;  // 0=auto
  bool ratchet_adaptive = false;

  // NNI perturbation: topology-space escape mechanism (IQ-TREE-style).
  // Randomly NNI-swap a fraction of internal branches, then TBR to a
  // new local optimum.  Complementary to weight-perturbation ratchet.
  int nni_perturb_cycles = 0;           // 0 = disabled
  double nni_perturb_fraction = 0.5;    // fraction of branches to perturb

  // Drifting
  int drift_cycles = 2;
  int drift_afd_limit = 3;
  double drift_rfd_limit = 0.1;

  // Simulated annealing: single-chain linear cooling schedule using
  // stochastic TBR.  Runs between drift and final TBR polish.
  // n_phases = 0 disables annealing.
  int anneal_phases = 0;           // 0 = disabled
  double anneal_t_start = 20.0;    // initial Boltzmann temperature
  double anneal_t_end = 0.0;       // final temperature
  int anneal_moves_per_phase = 0;  // 0 = n_tip

  // Sectorial search
  int xss_rounds = 3;
  int xss_partitions = 4;
  int rss_rounds = 1;           // RSS rounds after XSS; 0 = skip
  int css_rounds = 0;           // CSS rounds after RSS; 0 = skip
  int css_partitions = 4;       // partitions for CSS
  int sector_min_size = 6;
  int sector_max_size = 50;

  // Tree fusing
  int fuse_interval = 3;         // fuse every N replicates
  bool fuse_accept_equal = false;

  // Pool
  int pool_max_size = 100;
  double pool_suboptimal = 0.0;  // 0 = keep only optimal

  // Timeout (seconds). 0 or negative = no timeout.
  double max_seconds = 0.0;

  // Verbosity: 0 = silent, 1 = per-replicate summary, 2 = per-phase detail
  int verbosity = 0;

  // Progress callback. When set, invoked instead of Rprintf for progress
  // reporting. When empty, falls back to Rprintf.
  std::function<void(const ProgressInfo&)> progress_callback;

  // Tabu list size for TBR plateau exploration (0 = disabled)
  int tabu_size = 100;

  // NNI warmup: run NNI hill-climbing before SPR/TBR.
  // At ≤88 tips, overhead is negligible (~1.5s at 180 tips, <0.1s at ≤88).
  // At ≥180 tips, NNI saves ~50% of initial descent time and leads TBR to
  // better basins of attraction (empirically ~100 steps better at 180 tips).
  bool nni_first = true;

  // SPR→TBR escalation: run SPR first (cheaper per move), then TBR.
  // When true, initial hill-climbing is SPR followed by TBR to escape
  // moves that SPR cannot find. When false, goes straight to TBR.
  // Default false: TBR-only often gives better final scores because
  // the shallower initial optimum gives ratchet/drift more room to explore.
  bool spr_first = false;

  // Number of random Wagner trees per replicate (keep best-scoring)
  int wagner_starts = 1;

  // Biased taxon-addition order for Wagner tree construction.
  // 0 = RANDOM (default), 1 = GOLOBOFF (non-ambiguous chars), 2 = ENTROPY.
  // Applied only to the first Wagner start; remaining starts use random order
  // to preserve basin diversity.  Goloboff 2014 §3.3.
  int wagner_bias = 0;
  double wagner_bias_temp = 0.3;   // softmax temperature; 0 = greedy argmax

  // Outer search cycle count: number of times the [XSS → Ratchet →
  // NNI-perturb → Drift → TBR] block is repeated per replicate.
  // Default 1 = current linear pipeline.  Values > 1 interleave fresh
  // XSS passes after each ratchet/drift escape, matching TNT's xmult
  // pattern.  Ratchet/drift/NNI-perturb cycles are divided evenly among
  // outer cycles; total budget is approximately unchanged.
  // Goloboff 1999 §2.3 (sectorial + ratchet interleaving).
  int outer_cycles = 1;

  // Optional starting tree edge matrix (R format: n_edge × 2, 1-based).
  // When non-empty, replicate 0 uses this topology instead of Wagner.
  // Subsequent replicates still use random Wagner trees.
  std::vector<int> start_edge;  // flattened column-major [parent|child]
  int start_n_edge = 0;

  // Consensus-stability stopping criterion.
  // 0 = disabled (default). When > 0, stop if the strict consensus of
  // best-score pool trees has been unchanged for this many consecutive
  // replicates. Checked after each replicate completes and the pool is
  // updated. Sits alongside targetHits — whichever fires first wins.
  int consensus_stable_reps = 0;

  // Adaptive search level.
  // When true, dynamically scale ratchet_cycles and drift_cycles based
  // on the hit rate (fraction of replicates that find the current best
  // score). High hit rates → reduce effort; low hit rates → increase.
  // The base values are the initially configured cycles; adaptation
  // applies a multiplier each replicate.
  bool adaptive_level = false;

  // Cross-replicate consensus constraint tightening.
  // When true, after a minimum number of replicates, extract the strict
  // consensus splits from the pool and enforce them as topological
  // constraints for subsequent replicates. This focuses search on
  // uncertain parts of the tree. Constraints are cleared whenever the
  // best score improves. Only active when no user-supplied constraint
  // is present.
  bool consensus_constrain = false;
  int consensus_constrain_min_reps = 5;  // minimum replicates before engaging

  // Fraction of the time budget reserved for MPT enumeration (T-202).
  // The main search loop exits at budget × (1 - enum_time_fraction),
  // leaving the remainder for the plateau walk.  Default 0.1 = 10%.
  // Set to 0 to disable (old behaviour: skip enumeration on timeout).
  double enum_time_fraction = 0.1;

  // Adaptive starting-tree strategy selection (T-190).
  // When true, each replicate draws its starting strategy from a Thompson
  // sampling bandit over {Wagner-random, Wagner-Goloboff, Wagner-entropy,
  // random-tree, pool-ratchet, pool-NNI-perturb}. The bandit learns
  // which strategies hit the best score for this dataset.
  // When false, all replicates use the fixed `wagner_bias` strategy.
  // Only affects the serial path; parallel uses round-robin.
  bool adaptive_start = false;
};

// Cumulative per-phase wall-clock timing (milliseconds).
struct PhaseTimings {
  double wagner_ms = 0.0;
  double nni_ms = 0.0;
  double tbr_ms = 0.0;
  double xss_ms = 0.0;
  double rss_ms = 0.0;
  double css_ms = 0.0;
  double ratchet_ms = 0.0;
  double nni_perturb_ms = 0.0;
  double drift_ms = 0.0;
  double anneal_ms = 0.0;
  double final_tbr_ms = 0.0;
  double fuse_ms = 0.0;

  void operator+=(const PhaseTimings& o) {
    wagner_ms    += o.wagner_ms;
    nni_ms       += o.nni_ms;
    tbr_ms       += o.tbr_ms;
    xss_ms       += o.xss_ms;
    rss_ms       += o.rss_ms;
    css_ms       += o.css_ms;
    ratchet_ms   += o.ratchet_ms;
    nni_perturb_ms += o.nni_perturb_ms;
    drift_ms     += o.drift_ms;
    anneal_ms    += o.anneal_ms;
    final_tbr_ms += o.final_tbr_ms;
    fuse_ms      += o.fuse_ms;
  }
};

struct DrivenResult {
  double best_score;
  int replicates_completed;
  int hits_to_best;
  int pool_size;
  int n_topologies_at_best;      // distinct topologies at best score
  int last_improved_rep;         // 1-based replicate that last improved score (0 = not tracked)
  bool timed_out;                // true if search ended due to timeout
  bool consensus_stable;         // true if stopped by consensus stability
  PhaseTimings timings;          // cumulative across all replicates

  // Per-strategy diagnostics (populated when adaptive_start is true)
  std::array<int, N_STRAT> strategy_attempts{};
  std::array<int, N_STRAT> strategy_successes{};
};

// Result of a single replicate (tree + score, no pool interaction).
struct ReplicateResult {
  TreeState tree;
  double score;
  bool interrupted;  // true if stopped by interrupt or timeout
  PhaseTimings timings;          // per-replicate phase timings
};

// Run one replicate: Wagner → NNI → SPR → TBR → XSS → RSS → ratchet → drift → TBR.
// Does NOT interact with the pool — caller handles that.
// `check_timeout` should return true when time limit is exceeded.
// Verbosity is the effective verbosity for this replicate (0 in parallel).
struct SplitFrequencyTable;  // forward declaration (defined in ts_pool.h)

// If `starting_tree` is non-null, use it instead of building a Wagner tree.
// If `split_freq` is non-null, RSS uses conflict-guided sector selection.
// `strategy` controls how the starting tree is built when `starting_tree`
// is null. For pool-based strategies, the caller should perturb a pool tree
// and pass it as `starting_tree`.
ReplicateResult run_single_replicate(
    DataSet& ds,
    const DrivenParams& params,
    ConstraintData* cd,
    std::function<bool()> check_timeout,
    int verbosity,
    TreeState* starting_tree = nullptr,
    const SplitFrequencyTable* split_freq = nullptr,
    StartStrategy strategy = StartStrategy::WAGNER_RANDOM);

// Run the full driven search. Returns search statistics.
// The pool contents (all retained trees) are accessible via the pool
// reference stored in `pool_out`. Caller should extract edge matrices.
DrivenResult driven_search(TreePool& pool_out, DataSet& ds,
                           const DrivenParams& params,
                           ConstraintData* cd = nullptr);

} // namespace ts

#endif // TS_DRIVEN_H
