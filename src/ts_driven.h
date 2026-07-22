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

  // Taxon pruning-reinsertion (T-266): complementary perturbation that
  // drops a fraction of leaves, TBR-optimizes the backbone, then greedily
  // re-adds the dropped taxa via Wagner insertion + TBR polish.
  int prune_reinsert_cycles = 0;          // 0 = disabled
  double prune_reinsert_drop = 0.10;      // fraction of tips to drop
  int prune_reinsert_selection = 0;       // 0 = random, 1 = instability
  int prune_reinsert_tbr_moves = 5;       // TBR moves on reduced tree (0=converge)
  int prune_reinsert_full_moves = 0;      // TBR moves on full tree (0=converge)
  int prune_reinsert_nni = 0;            // 1 = NNI polish on full tree (cheaper at large n)

  // Drifting
  int drift_cycles = 2;
  int drift_afd_limit = 3;
  double drift_rfd_limit = 0.1;

  // Simulated annealing perturbation (PCSA: post-convergence SA).
  // Multi-cycle SA with best-tree restart, inserted after drift phase.
  // Each cycle: SA cooling schedule -> TBR reconverge -> keep if improved.
  // Effective at escaping deep basins under EW at >=100 tips.
  // 0 = disabled (default). Typical: 3-5 cycles for large trees.
  int anneal_cycles = 0;
  int anneal_phases = 5;             // temperature steps per cycle
  double anneal_t_start = 20.0;      // initial Boltzmann temperature
  double anneal_t_end = 0.0;         // final temperature
  int anneal_moves_per_phase = 0;    // 0 = n_tip

  // Sectorial search
  int xss_rounds = 3;
  int xss_partitions = 4;
  int rss_rounds = 1;           // RSS rounds after XSS; 0 = skip
  int css_rounds = 0;           // CSS rounds after RSS; 0 = skip
  int css_partitions = 4;       // partitions for CSS
  int sector_min_size = 6;
  int sector_max_size = 50;
  int ras_starts = 1;           // RAS+TBR restarts per sector (1 = polish the
                                // existing subtree; >1 rebuilds it that many
                                // times and keeps best, per Goloboff 1999 RSS;
                                // TNT uses 3). Plumbs SectorParams::ras_starts.
  bool sector_accept_equal = false;  // accept equal-score sector resolutions
                                     // (Goloboff 2014 plateau traversal);
                                     // plumbs SectorParams::accept_equal.
  int sector_max_hits = 1;           // equal-length trees the internal sector TBR
                                     // holds while swapping (1 = old; TNT holds
                                     // many). Plumbs SectorParams::internal_max_hits.
  int sector_collapse_target = 0;    // >0: collapse a big selected clade's deep
                                     // sub-clades into ~this many composite
                                     // first-pass terminals (Goloboff 1999 coarse
                                     // sector). Plumbs SectorParams::collapse_target.
  int rss_picks_per_round = 0;       // 0 = auto (2*n_tip/avg_size, ~5). >0 sets
                                     // sector picks BETWEEN global-TBR rounds
                                     // (Goloboff 1999 sequential replacements;
                                     // TNT ~20-25). Plumbs
                                     // SectorParams::rss_picks_per_round.

  // TNT-style in-sector drift / combined analysis for LARGE sectors.
  // Solve big sectors (>= threshold real tips) by drift or RAS+drift+fuse,
  // reaching sector topologies plain RAS+TBR cannot escape (TNT `sectsch`
  // godrift/gocomb/drift/combstarts/fuse). 0 thresholds = off (default).
  int sector_go_drift = 0;           // Plumbs SectorParams::sector_go_drift.
  int sector_go_comb = 0;            // Plumbs SectorParams::sector_go_comb.
  int sector_drift_cycles = 5;       // Plumbs SectorParams::sector_drift_cycles.
  int sector_drift_afd = 3;          // Plumbs SectorParams::sector_drift_afd.
  double sector_drift_rfd = 0.1;     // Plumbs SectorParams::sector_drift_rfd.
  int sector_comb_starts = 3;        // Plumbs SectorParams::sector_comb_starts.
  int sector_fuse_rounds = 3;        // Plumbs SectorParams::sector_fuse_rounds.

  // Post-ratchet sectorial search (T-257).
  // When true, run XSS+RSS+CSS again after ratchet perturbation using the
  // same round counts and sector parameters.  TNT interleaves sectorial
  // search throughout each replicate; this approximates that pattern by
  // exploiting the new basin reached after ratchet before TBR polish.
  bool post_ratchet_sectorial = false;

  // Tree fusing
  int fuse_interval = 3;         // fuse every N replicates (between-replicate)
  bool fuse_accept_equal = false;

  // Intra-replicate fusing (T-258).
  // When true, fuse the current tree against pool donors after TBR polish
  // in each outer cycle.  TNT fuses within each replicate; this approximates
  // that pattern.  The pool is read-only — the fused tree replaces the
  // current replicate tree but is only added to the pool after the replicate
  // completes.  Requires pool.size() >= 2 to have meaningful donors.
  bool intra_fuse = false;

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

  // TBR clip ordering strategy (see ClipOrder enum in ts_tbr.h).
  // 0=RANDOM (default), 1=INV_WEIGHT, 2=TIPS_FIRST, 3=BUCKET,
  // 4=ANTI_TIP, 5=LARGE_FIRST
  int clip_order = 0;

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
  // Default 1 = single pass through the pipeline.  Values > 1 interleave
  // fresh XSS passes after each ratchet/drift escape, matching TNT's
  // xmult pattern.  Ratchet/drift/NNI-perturb cycles are divided evenly
  // among outer cycles; total budget is approximately unchanged.
  // Goloboff 1999 §2.3 (sectorial + ratchet interleaving).
  int outer_cycles = 1;

  // Maximum number of improvement-triggered resets of the outer cycle
  // counter.  When a cycle improves the score, the counter resets to 0
  // so the search keeps exploiting the new basin — but at most this many
  // times.  0 = no resets (outer_cycles is exact).  -1 = unlimited.
  // Default 0: outer_cycles controls the total number of cycles exactly.
  // Strategy presets may set higher values (e.g. 2–3) to allow productive
  // re-exploration after escaping local optima.
  int max_outer_resets = 0;

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


  // Perturbation-count stopping rule (IQ-TREE-style).
  // 0 = disabled. Default 2 (set in R SearchControl).
  // When > 0, stop after nTip * perturb_stop_factor consecutive
  // replicates that fail to improve the best score.  Resets on
  // every improvement.
  int perturb_stop_factor = 0;
  // Adaptive search level.
  // When true, dynamically scale ratchet_cycles and drift_cycles based
  // on the hit rate (fraction of replicates that find the current best
  // score). High hit rates → reduce effort; low hit rates → increase.
  // The base values are the initially configured cycles; adaptation
  // applies a multiplier each replicate.
  bool adaptive_level = false;

  // Adaptive ratchet perturbation probability (T-182).
  // When true, the ratchet perturbation probability is tapered across
  // replicates based on pool stability.  Early replicates (unstable pool)
  // use ratchet_perturb_prob at full strength; later replicates (high
  // hit rate, stable consensus) use a reduced probability for finer
  // local exploration.
  //
  // The taper factor is:  max(taper_floor, 1.0 - taper_strength * stability)
  // where stability = hits_to_best / replicates_completed.
  // The effective probability = ratchet_perturb_prob * taper_factor.
  bool ratchet_taper = false;
  double ratchet_taper_floor = 0.5;    // minimum taper factor (prob ≥ 50% of base)
  double ratchet_taper_strength = 0.6; // how aggressively to reduce (0..1)

  // Cross-replicate stall escalation (per-dataset adaptive perturbation).
  // When a run stalls (unsuccessful_reps >= ceil(nTip/10)), escalate ratchet
  // perturbation for subsequent replicates:
  //   ratchet_perturb_prob = min(0.5, base * stall_escalate_factor^k),
  //   k = floor((unsuccessful_reps - s0) / s0),
  // and force ratchet_adaptive = TRUE for that replicate. An improvement
  // resets unsuccessful_reps to 0, which restores the base strength. The
  // factor 1.0 (default) disables the whole rule — a true no-op.
  double stall_escalate_factor = 1.0;

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
  double prune_reinsert_ms = 0.0;
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
    prune_reinsert_ms += o.prune_reinsert_ms;
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
  bool perturb_stop;             // true if stopped by perturb_stop_factor
  PhaseTimings timings;          // cumulative across all replicates

  // Per-strategy diagnostics (populated when adaptive_start is true)
  std::array<int, N_STRAT> strategy_attempts{};
  std::array<int, N_STRAT> strategy_successes{};

  // Score from each completed replicate's local optimum, in order of
  // completion.  Used by ScoreSpectrum() for Chao1-style coverage estimation.
  std::vector<double> replicate_scores;

  // Total TBR/SPR-class candidate rearrangements evaluated across the whole
  // search (TNT "Total rearrangements examined" analogue). Serial path only;
  // 0 when run in parallel. See DataSet::n_candidates_evaluated.
  long long candidates_evaluated = 0;
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
    StartStrategy strategy = StartStrategy::WAGNER_RANDOM,
    const TreePool* pool = nullptr);

// Run the full driven search. Returns search statistics.
// The pool contents (all retained trees) are accessible via the pool
// reference stored in `pool_out`. Caller should extract edge matrices.
DrivenResult driven_search(TreePool& pool_out, DataSet& ds,
                           const DrivenParams& params,
                           ConstraintData* cd = nullptr);

} // namespace ts

#endif // TS_DRIVEN_H
