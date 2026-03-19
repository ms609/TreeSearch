#include "ts_ratchet.h"
#include "ts_tbr.h"
#include "ts_fitch.h"
#include "ts_rng.h"

#include <random>
#include <algorithm>
#include <cmath>
#include <R.h>

namespace ts {

namespace {

// --- Perturbation state save/restore ---

struct PerturbSnapshot {
  std::vector<uint64_t> active_masks;
  std::vector<uint64_t> upweight_masks;
  std::vector<int> pattern_freq;  // for IW
};

void save_perturb_state(const DataSet& ds, PerturbSnapshot& snap) {
  snap.active_masks.resize(ds.n_blocks);
  snap.upweight_masks.resize(ds.n_blocks);
  for (int b = 0; b < ds.n_blocks; ++b) {
    snap.active_masks[b] = ds.blocks[b].active_mask;
    snap.upweight_masks[b] = ds.blocks[b].upweight_mask;
  }
  snap.pattern_freq = ds.pattern_freq;
}

void restore_perturb_state(DataSet& ds, const PerturbSnapshot& snap) {
  for (int b = 0; b < ds.n_blocks; ++b) {
    ds.blocks[b].active_mask = snap.active_masks[b];
    ds.blocks[b].upweight_mask = snap.upweight_masks[b];
  }
  ds.pattern_freq = snap.pattern_freq;
}

// --- Perturbation modes ---

// ZERO_ONLY: clear random active_mask bits (original Nixon/Goloboff approach).
void perturb_zero(DataSet& ds, double prob, std::mt19937& rng) {
  std::bernoulli_distribution coin(prob);
  for (auto& blk : ds.blocks) {
    uint64_t mask = blk.active_mask;
    for (int i = 0; i < blk.n_chars; ++i) {
      if (coin(rng)) {
        mask &= ~(uint64_t(1) << i);
      }
    }
    blk.active_mask = mask;
  }
}

// UPWEIGHT_ONLY: set upweight_mask bits (double selected characters).
// For IW/profile, also increment pattern_freq for selected patterns
// (each upweighted character contributes one extra copy to the pattern
// frequency, matching the EW behaviour where upweight_mask adds one
// extra step per selected bit position).
void perturb_upweight(DataSet& ds, double prob, std::mt19937& rng,
                      bool use_iw) {
  std::bernoulli_distribution coin(prob);
  for (auto& blk : ds.blocks) {
    uint64_t up = 0;
    for (int i = 0; i < blk.n_chars; ++i) {
      if (coin(rng)) {
        up |= (uint64_t(1) << i);
        if (use_iw) {
          int pat = blk.pattern_index[i];
          if (pat >= 0) ds.pattern_freq[pat] += 1;
        }
      }
    }
    blk.upweight_mask = up & blk.active_mask;
  }
}

// MIXED: zero some characters, upweight others (disjoint sets).
// Each character is independently: zeroed with prob p, or upweighted with
// prob p, or left unchanged with prob 1-2p.  If 2p > 1, clamp.
void perturb_mixed(DataSet& ds, double prob, std::mt19937& rng,
                   bool use_iw) {
  double p_zero = prob;
  double p_up = prob;
  // Clamp so total probability ≤ 1
  if (p_zero + p_up > 1.0) {
    p_zero = 0.5;
    p_up = 0.5;
  }

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  for (auto& blk : ds.blocks) {
    uint64_t new_active = blk.active_mask;
    uint64_t up = 0;
    for (int i = 0; i < blk.n_chars; ++i) {
      if (!(blk.active_mask & (uint64_t(1) << i))) continue;
      double r = unif(rng);
      if (r < p_zero) {
        new_active &= ~(uint64_t(1) << i);
      } else if (r < p_zero + p_up) {
        up |= (uint64_t(1) << i);
        if (use_iw) {
          int pat = blk.pattern_index[i];
          if (pat >= 0) ds.pattern_freq[pat] += 1;
        }
      }
    }
    blk.active_mask = new_active;
    blk.upweight_mask = up & new_active;
  }
}

// Copy topology from src to dst (same-sized trees).
void copy_topology(TreeState& dst, const TreeState& src) {
  dst.parent = src.parent;
  dst.left   = src.left;
  dst.right  = src.right;
}

} // anonymous namespace

RatchetResult ratchet_search(TreeState& tree, DataSet& ds,
                             const RatchetParams& params,
                             ConstraintData* cd,
                             std::function<bool()> check_timeout) {
  const bool use_iw = std::isfinite(ds.concavity);

  // Initial TBR to get a baseline
  TBRParams search_params;
  search_params.accept_equal = false;
  search_params.max_accepted_changes = 0;
  search_params.max_hits = params.max_hits;
  search_params.tabu_size = params.tabu_size;

  TBRResult initial = tbr_search(tree, ds, search_params, cd);

  double best_score = initial.best_score;
  int total_moves = initial.n_accepted;
  int n_escapes = 0;

  // Save best topology
  TreeState best_tree = tree;

  // Perturbation TBR params
  int actual_max_moves = params.perturb_max_moves > 0
      ? params.perturb_max_moves
      : std::max(20, std::min(200, tree.n_tip / 8));

  TBRParams perturb_params;
  perturb_params.accept_equal = params.perturb_accept_equal;
  perturb_params.max_accepted_changes = actual_max_moves;
  perturb_params.max_hits = 1;
  perturb_params.tabu_size = params.tabu_size;

  // Seed RNG (from R in serial mode, from thread-local in parallel mode)
  std::mt19937 rng = ts::make_rng();

  // Adaptive perturbation state
  double current_prob = params.perturb_prob;
  const int adapt_batch = 3;
  int recent_escapes = 0;

  PerturbSnapshot snap;
  int cycles_completed = 0;

  for (int cycle = 0; cycle < params.n_cycles; ++cycle) {
    // 1. Perturbation phase: modify weights
    save_perturb_state(ds, snap);

    switch (params.perturb_mode) {
      case PerturbMode::ZERO_ONLY:
        perturb_zero(ds, current_prob, rng);
        break;
      case PerturbMode::UPWEIGHT_ONLY:
        perturb_upweight(ds, current_prob, rng, use_iw);
        break;
      case PerturbMode::MIXED:
        perturb_mixed(ds, current_prob, rng, use_iw);
        break;
    }

    // 2. Short TBR on perturbed landscape
    TBRResult perturb_result = tbr_search(tree, ds, perturb_params, cd);
    total_moves += perturb_result.n_accepted;

    // 3. Restore original weights, full TBR to new local optimum
    restore_perturb_state(ds, snap);
    TBRResult search_result = tbr_search(tree, ds, search_params, cd);
    total_moves += search_result.n_accepted;

    if (search_result.best_score < best_score) {
      best_score = search_result.best_score;
      best_tree = tree;
      ++n_escapes;
      ++recent_escapes;
    } else {
      // Reset to best known tree
      copy_topology(tree, best_tree);
      tree.build_postorder();
      tree.reset_states(ds);
    }

    ++cycles_completed;

    // 4. Adaptive tuning (opt-in)
    if (params.adaptive && cycle > 0 && (cycle + 1) % adapt_batch == 0) {
      double escape_rate =
          static_cast<double>(recent_escapes) / adapt_batch;
      if (escape_rate < params.target_escape_rate * 0.5) {
        // Not escaping enough — perturb harder
        current_prob = std::min(params.adapt_max_prob, current_prob * 1.5);
      } else if (escape_rate > params.target_escape_rate * 2.0) {
        // Escaping too easily — may be over-disrupting
        current_prob = std::max(params.adapt_min_prob, current_prob * 0.7);
      }
      recent_escapes = 0;
    }

    if (ts::check_interrupt()) break;
    if (check_timeout && check_timeout()) break;
  }

  // Ensure tree holds the best result
  if (cycles_completed > 0) {
    copy_topology(tree, best_tree);
    tree.build_postorder();
    tree.reset_states(ds);
  }

  return RatchetResult{
    best_score,
    cycles_completed,
    total_moves,
    n_escapes,
    current_prob
  };
}

} // namespace ts
